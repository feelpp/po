#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "eigenlap.h"

EigenLap::EigenLap( mesh_ptrtype mesh ):super()
{
    this->nev = option(_name="solvereigen.nev").template as<int>();
    this->ncv = option(_name="solvereigen.ncv").template as<int>();
    this->mesh = mesh;
    this->Vh = vSpace_type::New( mesh );
    this->Mlh = mlSpace_type::New( mesh );
    this->g = std::vector<vElement_type>(nev, Vh->element() );
    this->psi = std::vector<mlElement_type>(nev, Mlh->element() );
    this->lambda = std::vector<double>(nev, 0 );

    this->g0 = std::vector<vElement_type>(nev, Vh->element() );
    this->gradu = std::vector<vElement_type>(nev, Vh->element() );
    this->modebis = std::vector<vElement_type>(nev, Vh->element() );

}

void
EigenLap::run()
{
    if( option( _name="needEigen").as<bool>() )
        compute_eigens();
    else
        load_eigens();
    if( option( _name="needDecomp").as<bool>() )
        decomp();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "-----End Eigen-----" << std::endl;
}

void
EigenLap::compute_eigens()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "-----Eigen Problem-----" << std::endl;
        std::cout << "number of eigenvalues computed = " << nev <<std::endl;
        std::cout << "number of column vector = " << ncv <<std::endl;
    }

    int penal = option(_name="penal").as<double>();

    auto Xh = sSpace_type::New( mesh );
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto u3 = U.template element<2>();
    auto V = Xh->element();
    auto v1 = V.template element<0>();
    auto v2 = V.template element<1>();
    auto v3 = V.template element<2>();

    auto p = U.template element<3>();
    auto q = V.template element<3>();

    auto l = form1( _test=Xh );

    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ), 
                   _expr=dxt(u1)*dx(v1) + dyt(u1)*dy(v1) + dzt(u1)*dz(v1)
                   + dxt(u2)*dx(v2) + dyt(u2)*dy(v2) + dzt(u2)*dz(v2)
                   + dxt(u3)*dx(v3) + dyt(u3)*dy(v3) + dzt(u3)*dz(v3)
                   + (dxt(u1)+dyt(u2)+dzt(u3))*id(q)
                   + (dx(v1)+dy(v2)+dz(v3))*idt(p)
                   + 1e-20*idt(p)*id(q) );

    a += on(_range=markedfaces(mesh,2), _rhs=l, _element=p, _expr=cst(0.));
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u1, _expr=cst(0.));
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u2, _expr=cst(0.));
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u3, _expr=cst(0.));

    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), idt( u1 )*id( v1 )
                   + idt( u2 )*id( v2 )
                   + idt( u3 )*id( v3 )                   
                   + 0.*idt(p)*id(q) );

    SolverEigen<double>::eigenmodes_type modes;
    modes = eigs( _matrixA=a.matrixPtr(),
                  _matrixB=b.matrixPtr(),
                  _nev=nev,
                  _ncv=ncv,
                  _transform=SINVERT,
                  _spectrum=SMALLEST_MAGNITUDE,
                  _verbose = true );

    auto modeTmp = Xh->element();

    if ( !modes.empty() )
    {
        int i = 0;
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("lambda", std::fstream::out);
        for( auto const& mode : modes )
        {
            modeTmp = *mode.second.get<2>();
            g[i] = vf::project(_space=Vh, _range=elements(mesh),
                               _expr=vec(idv(modeTmp.template element<0>()),
                                      idv(modeTmp.template element<1>()),
                                      idv(modeTmp.template element<2>()) ) );
            std::string path = (boost::format("mode-%1%")%i).str();
            g[i].save(_path=path);
            lambda[i] = mode.second.get<0>();
            if ( Environment::worldComm().isMasterRank() )
                s << lambda[i] << std::endl;

            double erreurEigS = normL2( _range=elements(mesh),
                                        _expr=curlv(g[i])-sqrt(lambda[i])*idv(g[i]) );

            auto cg = Vh->element();
            auto c = form2( _trial=Vh, _test=Vh );
            c = integrate( _range=elements(mesh), _expr=inner(idt(cg),id(cg)) );
            auto f = form1( _test=Vh );
            f = integrate( _range=elements(mesh), _expr=inner(curlv(g[i]),id(cg)) );
            c.solve( _name="curl", _rhs=f, _solution=cg );
            double erreurEig = normL2( _range=elements(mesh),
                                       _expr=curlv(cg)-lambda[i]*idv(g[i]) );

            double di = normL2(elements(mesh), divv(g[i]));
            double bord = normL2(boundaryfaces(mesh), trans(curlv(g[i]))*N() );

            if ( Environment::worldComm().isMasterRank() ){
                std::cout << i << " : " << lambda[i] << std::endl;
                std::cout << "curl-sqrt(lambda) = " << erreurEigS << " curl-lambda = " << erreurEig << " div = " << di << " curl(g).n = " << bord << std::endl;
            }

            i++;
            if(i>=nev)
                break;
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
}

void
EigenLap::load_eigens()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "-----Load Eigen-----" << std::endl;
        std::cout << "number of eigenvalues = " << nev <<std::endl;
    }

    std::fstream s;
    s.open ("lambda", std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found\ntry to launch with --needEigen=true" << std::endl;
        exit(0);
    }

    int i;
    for( i=0; i<nev && s.good(); i++ ){
        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].load(_path=path);
        s >> lambda[i];
    }

    s.close();

    if ( i != nev ){
        std::cout << "Number of eigenvalues different from nev !" << std::endl;
        exit(0);
    }
}

void
EigenLap::decomp()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "-----Decomposition-----" << std::endl;

    auto a2 = form2( _test=Vh, _trial=Vh );
    auto l2 = form1( _test=Vh );

    auto a3 = form2( _test=Mlh, _trial=Mlh );
    auto l3 = form1( _test=Mlh );

    for(int i=0; i<nev; i++){
        auto vg0 = Vh->element();
        l2 = integrate( _range=elements(mesh),
                        _expr=lambda[i]*inner(idv(g[i]),id(vg0)) );
        a2 = integrate( _range=elements(mesh),
                        _expr=-divt(g0[i])*div(vg0) );
        a2+= integrate( _range=elements(mesh),
                        _expr=inner(gradt(g0[i]), grad(vg0) ) );
        a2+= on( _range=boundaryfaces(mesh),
                 _element=g0[i], _rhs=l2, _expr=cst(0.) );
        a2.solve( _name="gi0", _rhs=l2, _solution=g0[i] );


        auto psii = psi[i].template element<0>();
        auto psil = psi[i].template element<1>();
        a3 = integrate( _range=elements(mesh),
                        _expr=inner(gradt(psii),grad(psii))
                        );//+ id(psii)*idt(psil) + idt(psii)*id(psil) );
        l3 = integrate( _range=elements(mesh), _expr=divv(g0[i])*id(psii) );
        a3.solve( _name="psi", _rhs=l3, _solution=psi[i] );


        auto gv = Vh->element();
        auto c = form2( _trial=Vh, _test=Vh );
        c = integrate( _range=elements(mesh), _expr=inner(idt(gradu[i]),id(gv)) );
        auto f = form1( _test=Vh );
        f = integrate( _range=elements(mesh), _expr=gradv( psi[i].template element<0>() )*id(gv));
        c.solve( _name="gradpsi", _rhs=f, _solution=gradu[i] );

        modebis[i] = vf::project( _space=Vh, _range=elements(mesh),
                                  _expr=idv(g0[i])+idv(gradu[i]) );

        double erreurL2 = normL2( elements(mesh), idv(g[i])-idv(modebis[i]) );
        double nG = normL2( elements(mesh), idv(g[i]) );
        double nG0 = normL2( elements(mesh), idv(g0[i]) );
        double nPsi = normL2( elements(mesh), idv(gradu[i]) );
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "||g-(g0+grad(psi)||_L2 = " << erreurL2 << " ||g|| = " << nG << " ||g0|| = " << nG0 << " ||psi|| = " << nPsi << std::endl;
    }
}
