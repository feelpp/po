#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

#include "eigenlapZ.h"

EigenLapZ::EigenLapZ( mesh_ptrtype mesh ):super()
{
    this->nev = ioption(_name="solvereigen.nev");
    this->ncv = ioption(_name="solvereigen.ncv");
    this->mesh = mesh;
    this->Vh = vSpace_type::New( mesh );
    this->Sh = sSpace_type::New( mesh );
    this->g = std::vector<vElement_type>(nev, Vh->element() );
    this->psi = std::vector<sElement_type>(nev, Sh->element() );
    this->lambda = std::vector<double>(nev, 0 );

    this->g0 = std::vector<vElement_type>(nev, Vh->element() );
    this->gradu = std::vector<vElement_type>(nev, Vh->element() );
    this->modebis = std::vector<vElement_type>(nev, Vh->element() );

}

void
EigenLapZ::run()
{
    if( boption( _name="needEigen") )
        compute_eigens();
    else
        load_eigens();

    if( boption( _name="testBessel") )
        testBessel();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End Eigen -----" << std::endl;
}

void
EigenLapZ::compute_eigens()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Eigen Problem -----" << std::endl;
        std::cout << "number of eigenvalues computed = " << nev <<std::endl;
        std::cout << "number of column vector = " << ncv <<std::endl;
    }

    int bc = doption(_name="bc");


    auto MLh = mlSpace_type::New(mesh);
    auto upi = MLh->element();
    auto vxi = MLh->element();

    auto u3 = upi.template element<0>();
    auto pi = upi.template element<1>();
    auto v3 = vxi.template element<0>();
    auto xi = vxi.template element<1>();

    auto l = form1( _test=MLh );

    auto a = form2( _test=MLh, _trial=MLh);
    // [formA]
    a = integrate( _range=elements( mesh ),
                   _expr= dxt(u3)*dx(v3) + dyt(u3)*dy(v3) + dzt(u3)*dz(v3) );
    // [formA]
    if( boption( _name="needRelev" )){
        // [relev]
        a += integrate( _range=elements(mesh), _expr=1.e-20*idt(pi)*id(xi) );
        a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u3, _expr=cst(bc));
        a += on(_range=markedfaces(mesh, 2), _rhs=l, _element=u3, _expr=cst(bc));
        // [relev]
    } else {
        // [nrelev]
        a += integrate( _range=elements( mesh ),
                        _expr=dzt(u3)*id(xi) + dz(v3)*idt(pi) );
        a += on(_range=markedfaces(mesh, 3), _rhs=l, _element=u3, _expr=cst(bc));
        // [nrelev]
    }

    // [formB]
    auto b = form2( _test=MLh, _trial=MLh);
    b = integrate( _range=elements(mesh),
                   _expr=idt( u3 )*id( v3 ) + 0.*idt(pi)*id(xi) );
    // [formB]

    auto modes = veigs( _formA=a, _formB=b );

    int i = 0;
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);

    auto gd = Sh->element();

    auto a2 = form2( _test=Sh, _trial=Sh );
    auto l2 = form1( _test=Sh );

    auto a3 = form2( _test=Sh, _trial=Sh );
    auto l3 = form1( _test=Sh );

    for( auto const& mode: modes )
    {
        u3 = mode.second.element<0>();
        g[i] = vf::project(_space=Vh, _range=elements(mesh),
                           _expr=vec(cst(0.),
                                     cst(0.),
                                     idv(u3) ) );

        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].save(_path=path);
        lambda[i] = mode.first;
        if ( Environment::worldComm().isMasterRank() )
            s << lambda[i] << std::endl;

        auto normL2Div = normL2( _range=elements(mesh), _expr=divv(g[i]) );
        auto div = integrate(_range=elements(mesh), _expr=divv(g[i]) ).evaluate()(0,0);

        if ( Environment::isMasterRank() )
            std::cout << "Lambda_" << i << " = " <<  mode.first << "\tDivergence = " << div << "\t||Div|| = " << normL2Div << "\n";

        if( boption(_name="needDecomp") ){
            // [gi0]
            l2 = integrate( _range=elements(mesh),
                            _expr=mode.first*idv(u3)*id(gd) );
            a2 = integrate( _range=elements(mesh),
                            _expr=-dzt(gd)*dz(gd)
                            + dxt(gd)*dx(gd) + dyt(gd)*dy(gd) + dzt(gd)*dz(gd) );
            a2 += on(_range=boundaryfaces(mesh), _rhs=l2, _element=gd, _expr=cst(0.) );
            // [gi0]
            a2.solve(_name="gi0", _rhs=l2, _solution=gd);
            g0[i] = vf::project(_space=Vh, _range=elements(mesh),
                                _expr=vec(cst(0.),
                                          cst(0.),
                                          idv(gd) ) );

            // [psi]
            l3 = integrate(_range=elements(mesh),
                           _expr=dzv(gd)*id(gd) );
            a3 = integrate(_range=elements(mesh),
                           _expr=dxt(psi[i])*dx(psi[i]) + dyt(psi[i])*dy(psi[i]) + dzt(psi[i])*dz(psi[i]) );
            // [psi]
            a3.solve(_name="psi", _rhs=l3, _solution=psi[i]);

            auto gv = Vh->element();
            auto c = form2( _trial=Vh, _test=Vh );
            c = integrate( _range=elements(mesh), _expr=inner(idt(gradu[i]),id(gv)) );
            auto f = form1( _test=Vh );
            f = integrate( _range=elements(mesh), _expr=gradv( psi[i] )*id(gv));
            c.solve( _name="gradpsi", _rhs=f, _solution=gradu[i] );

            modebis[i] = vf::project( _space=Vh, _range=elements(mesh),
                                      _expr=idv(g0[i])+idv(gradu[i]) );

            double erreurL2 = normL2( elements(mesh), idv(g[i])-idv(modebis[i]) );
            double nG = normL2( elements(mesh), idv(g[i]) );
            double nG0 = normL2( elements(mesh), idv(g0[i]) );
            double nPsi = normL2( elements(mesh), idv(gradu[i]) );
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "||g-(g0+grad(psi)||_L2 = " << erreurL2 << " ||g|| = " << nG << " ||g0|| = " << nG0 << " ||psi|| = " << nPsi << std::endl;

            std::string pathPsi = (boost::format("psi-%1%")%i).str();
            psi[i].save(_path=pathPsi);

        }
        i++;
        if(i>=nev)
            break;
    }
    if ( Environment::worldComm().isMasterRank() )
        s.close();

}

void
EigenLapZ::load_eigens()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Load Eigen -----" << std::endl;
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
        std::string pathPsi = (boost::format("psi-%1%")%i).str();
        psi[i].load(_path=pathPsi);
        s >> lambda[i];
    }

    s.close();

    if ( i != nev ){
        std::cout << "Number of eigenvalues different from nev !" << std::endl;
        exit(0);
    }
}

void EigenLapZ::testBessel()
{
    double l=0.5;
    double R=0.05;
    double k = doption(_name="k");
    int j = ioption(_name="j");
    int m = ioption(_name="m");

    // unsigned int n_roots = 5U;
    // std::vector<double> roots;
    // for(int i=0; i<5; i++){
    //     std::cout << "k = " << i << std::endl;
    //     boost::math::cyl_bessel_j_zero(0.0+i, 1, n_roots, std::back_inserter(roots));
    //     if ( Environment::worldComm().isMasterRank() ){
    //         std::copy(roots.begin(),
    //                   roots.end(),
    //                   std::ostream_iterator<double>(std::cout, "\n"));
    //     }
    //     auto s = std::sqrt(2)/(std::sqrt(l*pi)*R*std::abs(-boost::math::cyl_bessel_j(1.0+i, roots[0])))*boost::math::cyl_bessel_j(0.0+i, 0)*std::sin(pi/l*0);
    //     if ( Environment::worldComm().isMasterRank() ){
    //         std::cout << "g3_{"<<i<<",0,0}(0,0,0) = " << s << std::endl;
    //     }
    // }

    auto zero = boost::math::cyl_bessel_j_zero(k,j);
    double lam = zero*zero/(R*R)+m*m*pi*pi/(l*l);
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "zero(" << k << "," << j << ") = " << zero << "\tLambda(" << k << "," << j << "," << m << ") = " << lam << std::endl;


    psi[0] = vf::project(_space=Sh,_range=elements(mesh),
                       _expr=std::sqrt(2)/(std::sqrt(l*pi)*R*std::abs(-boost::math::cyl_bessel_j(k+1, zero)))*sin(m*pi/l*Pz()) );

    //*boost::math::cyl_bessel_j(0.0,zero*sqrt(Px()*Px()+Py()*Py())/R )

}
