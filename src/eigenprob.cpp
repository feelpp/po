#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/projectors.hpp>

#include "eigenprob.hpp"

EigenProb::EigenProb( mesh_ptrtype mesh ):super()
{
    this->nev = ioption(_name="solvereigen.nev");

    this->mesh = mesh;

    // L2^3
    this->Vh = space_vtype::New( mesh );
    // L2
    this->Sh = space_stype::New( mesh );

    this->g = std::vector<element_vtype>(nev, Vh->element() );
    this->psi = std::vector<element_stype>(nev, Sh->element() );
    this->lambda = std::vector<double>(nev, 0 );

    // debug
    this->g0 = std::vector<element_vtype>(nev, Vh->element() );
    this->gradu = std::vector<element_vtype>(nev, Vh->element() );
    this->modebis = std::vector<element_vtype>(nev, Vh->element() );

}

void
EigenProb::run()
{
    if( boption( _name="needEigen") )
        compute_eigens();
    else
        load_eigens();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End Eigen -----" << std::endl;
}

void
EigenProb::compute_eigens()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Eigen Problem -----" << std::endl;
        std::cout << "number of eigenvalues computed = " << nev <<std::endl;
        std::cout << "number of column vector = " << ioption(_name="solvereigen.ncv") <<std::endl;
    }

    double alpha = doption(_name="parameters.alpha");
    double beta = doption(_name="parameters.beta");
    double gamma = doption(_name="parameters.gamma");

    // L2^3xL2
    auto Xh = space_ptype::New( mesh );
    auto U = Xh->element();
    auto V = Xh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    // eigen problem Ax=lambda*Bx
    auto l = form1( _test=Xh );
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ),
                   trans(curlt(u))*curl(v)
                   + grad(q)*idt(u) + gradt(p)*id(v)
                   );
    if( boption(_name="divdiv") )
        a += integrate( elements(mesh), alpha*divt(u)*div(v) );
    if ( boption(_name="bccurln" ) )
        a += integrate( boundaryfaces(mesh), beta*(trans(curlt(u))*N())*(trans(curl(v))*N())/hFace() );
    if( boption(_name="bcn" ) )
        a += integrate( boundaryfaces(mesh), gamma*(trans(idt(u))*N())*(trans(id(v))*N())/hFace() );

    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), inner(idt(u),id(v)) );

    auto modes = veigs( _formA=a, _formB=b );


    // preparation for the decomposition
    auto vg0 = Vh->element();
    auto a2 = form2( _test=Vh, _trial=Vh );
    auto l2 = form1( _test=Vh );

    // L2xR
    auto Mlh = space_mltype::New( mesh );
    auto Psi = Mlh->element();
    auto psii = Psi.element<0>();
    auto nu = Psi.element<1>();
    auto a3 = form2( _test=Mlh, _trial=Mlh );
    auto l3 = form1( _test=Mlh );

    if( boption( _name="needDecomp") ){
        // left side of gi0
        a2 = integrate( _range=elements(mesh),
                        _expr=trans(curlt(vg0))*curl(vg0) );

        // left side of psi
        a3 = integrate( _range=elements(mesh),
                        _expr=inner(gradt(psii),grad(psii) )
                        + id(psii)*idt(nu) + idt(psii)*id(nu) );
    }

    // left side for grad psi
    auto gv = Vh->element();
    auto c = form2( _trial=Vh, _test=Vh );
    auto f = form1( _test=Vh );

    if( boption( _name="needDebug") && boption(_name="needDecomp") )
        c = integrate( _range=elements(mesh), _expr=inner(idt(gv),id(gv)) );


    // opening file to store lambda
    int i = 0;
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);


    for( auto const& mode : modes )
    {
        g[i] = mode.second.element<0>();
        lambda[i] = mode.first;

        // storing g and lambda
        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].save(_path=path);
        if ( Environment::worldComm().isMasterRank() )
            s << lambda[i] << std::endl;


        if( boption(_name="needDebug")){
            double di = normL2(elements(mesh), divv(g[i]));
            double gn = normL2(boundaryfaces(mesh), trans(idv(g[i]))*N());
            double curlgn = normL2(boundaryfaces(mesh), trans(curlv(g[i]))*N());
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "lambda_" << i << " = " << lambda[i] << "  //  ||div|| = " << di << "  //  ||g.n|| = " << gn << "  //  ||curlg.n|| = " << curlgn;
     }

        if ( boption( _name="needDecomp") ){
            // right side of gi0
            l2 = integrate( _range=elements(mesh),
                            _expr=lambda[i]*inner(idv(g[i]),id(vg0)) );
            a2+= on( _range=boundaryfaces(mesh),
                     _element=vg0, _rhs=l2, _expr=zero<3>() );
            a2.solve( _name="gi0", _rhs=l2, _solution=g0[i] );

            // right side of psi
            l3 = integrate( _range=elements(mesh), _expr=divv(g0[i])*id(psii) );
            a3.solve( _name="psi", _rhs=l3, _solution=Psi );
            psi[i] = Psi.element<0>();

            // storing psi
            std::string pathPsi = (boost::format("psi-%1%")%i).str();
            psi[i].save(_path=pathPsi);


            if( boption( _name="needDebug") ){
                // right side of grad psi
                f = integrate( _range=elements(mesh), _expr=gradv( psi[i] )*id(gv));
                c.solve( _name="gradpsi", _rhs=f, _solution=gradu[i] );

                // g0 + grad psi
                modebis[i] = vf::project( _space=Vh, _range=elements(mesh),
                                          _expr=idv(g0[i])+idv(gradu[i]) );

                double erreurL2 = normL2( elements(mesh), idv(g[i])-idv(modebis[i]) );
                double nG = normL2( elements(mesh), idv(g[i]) );
                double nG0 = normL2( elements(mesh), idv(g0[i]) );
                double nPsi = normL2( elements(mesh), idv(gradu[i]) );
                if ( Environment::worldComm().isMasterRank() )
                    std::cout << "    ||g-(g0+grad(psi)|| = " << erreurL2 << "    ||g|| = " << nG << "    ||g0|| = " << nG0 << "    ||psi|| = " << nPsi;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
             std::cout << std::endl;

        i++;
        if(i>=nev)
            break;
    }
    if ( Environment::worldComm().isMasterRank() )
        s.close();
}

void
EigenProb::load_eigens()
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
