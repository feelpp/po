/** \file eigenprob.cpp
    \brief Source file for the class EigenProb
*/

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/mpi/timer.hpp>

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
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/ginac.hpp>

#include "eigenprob.hpp"


EigenProb::EigenProb( mesh_ptrtype mesh ):super()
{
    this->nev = ioption(_name="solvereigen.nev");

    this->mesh = mesh;


    // H1^3
    this->Vh = space_vtype::New( mesh );
    // H1
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
    boost::mpi::timer t;

    if( boption( _name="computeEigen") ){

        compute_eigens();

        LOG(INFO) << "eigen = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "eigen = " << t.elapsed() << " sec" << std::endl;
        t.restart();

        if( boption( "needDecomp" ) ){

            compute_decomp();

            LOG(INFO) << "decomp = " << t.elapsed() << " sec" << std::endl;
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "decomp = " << t.elapsed() << " sec" << std::endl;
        }
    }
    else {
        load_eigens();
        if( boption( "needDecomp") ){
            if( boption( "computeDecomp") ){

                t.restart();
                compute_decomp();

                LOG(INFO) << "decomp = " << t.elapsed() << " sec" << std::endl;
                if ( Environment::worldComm().isMasterRank() )
                    std::cout << "decomp = " << t.elapsed() << " sec" << std::endl;
            }
            else {
                load_decomp();
            }
        }
    }

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End Eigen -----" << std::endl;
}


void
EigenProb::compute_eigens()
{
    LOG(INFO) << "----- Eigen Problem -----" << std::endl;
    LOG(INFO) << "number of eigenvalues computed = " << nev <<std::endl;
    LOG(INFO) << "number of column vector = " << ioption(_name="solvereigen.ncv") <<std::endl;

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Eigen Problem -----" << std::endl;
        std::cout << "number of eigenvalues computed = " << nev <<std::endl;
        std::cout << "number of column vector = " << ioption(_name="solvereigen.ncv") <<std::endl;
    }


    // [options]
    double alpha = doption(_name="parameters.alpha");
    double beta = doption(_name="parameters.beta");
    double gamma = doption(_name="parameters.gamma");
    // [options]


    // H1^3xH1
    auto Xh = space_ptype::New( mesh );

    LOG(INFO) << "----- Xh -----\n";
    LOG(INFO) << "[dof] number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(U): " << Xh->template functionSpace<0>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(P): " << Xh->template functionSpace<1>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Xh -----\n";
        std::cout << "[dof] number of dof: " << Xh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Xh->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(U): " << Xh->template functionSpace<0>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(P): " << Xh->template functionSpace<1>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof() << "\n";
    }


    // eigen problem Ax=lambda*Bx
    auto U = Xh->element();
    auto V = Xh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();


    auto l = form1( _test=Xh );
    auto a = form2( _test=Xh, _trial=Xh);
    auto b = form2( _test=Xh, _trial=Xh);

    // g[i] = mode.second.element<0>();


    // LOG(INFO) << "----- Vh -----\n";
    // LOG(INFO) << "[dof] number of dof: " << Vh->nDof() << "\n";
    // LOG(INFO) << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";

    // if ( Environment::isMasterRank() ){
    //     std::cout << "----- Vh -----\n";
    //     std::cout << "[dof] number of dof: " << Vh->nDof() << "\n";
    //     std::cout << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
    // }


    // auto u = Vh->element();
    // auto v = Vh->element();


    // auto l = form1( _test=Vh );
    // auto a = form2( _test=Vh, _trial=Vh);
    // auto b = form2( _test=Vh, _trial=Vh);


    if( boption( "useCurl" ) )
        // [acurl]
        a = integrate( elements( mesh ),
                       trans(curlt(u))*curl(v) );
        // [acurl]
    else
        // [agrad]
        a = integrate( elements( mesh ),
                       inner(gradt(u), grad(v)) );
        // [agrad]

    if( boption( "usePresDiv" ) )
        // [presdiv]
        a += integrate( elements(mesh),
                        divt(u)*id(q) + div(v)*idt(p) );
        // [presdiv]

    if( boption( "usePresGrad" ) )
        // [presgrad]
        a += integrate( elements(mesh),
                        grad(q)*idt(u) + gradt(p)*id(v) );
        // [presgrad]


    if( boption("useDiric" ) )
        a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=zero<3>() );

    if( boption(_name="divdiv") )
        // [divdiv]
        a += integrate( elements(mesh), alpha*divt(u)*div(v) );
        // [divdiv]

    if ( boption(_name="bccurln" ) )
        // [bccurln]
        a += integrate( boundaryfaces(mesh), beta*(trans(curlt(u))*N())*(trans(curl(v))*N()) );
        // [bccurln]

    if( boption(_name="bcn" ) )
        // [bcn]
        a += integrate( boundaryfaces(mesh), gamma*(trans(idt(u))*N())*(trans(id(v))*N()) );
        // [bcn]


    // [rhsB]
    b = integrate( elements(mesh), inner(idt(u),id(v)) );
    // [rhsB]


    Environment::logMemoryUsage("memory usage before:");
    // [eigen] resolution of the eigen problem
    auto modes = veigs( _formA=a, _formB=b );
    // [eigen]
    Environment::logMemoryUsage("memory usage after:");


    auto cu = Vh->element();
    auto ccu = Vh->element();
    auto acu = form2( _test=Vh, _trial=Vh );
    auto bcu = form1( _test=Vh );

    if( boption("needDebug") )
        acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );


    // opening file to store lambda
    int i = 0;
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);


    for( auto const& mode : modes )
    {
        // bcu = integrate(elements(mesh),
        //                 inner(idv(mode.second.element<0>()), id(cu)) );
        // acu.solve(_rhs=bcu, _solution=g[i], _name="gi0" );

        g[i] = mode.second.element<0>();
        lambda[i] = mode.first;


        // [store] storing g and lambda
        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].save(_path=path);
        if ( Environment::worldComm().isMasterRank() )
            s << lambda[i] << std::endl;
        // [store]


        if( boption(_name="needDebug")){

            bcu = integrate(elements(mesh),
                            trans(curlv(g[i]))*id(cu) );
            acu.solve(_rhs=bcu, _solution=cu, _name="curl" );

            bcu = integrate(elements(mesh),
                            trans(curlv(cu))*id(cu) );
            acu.solve(_rhs=bcu, _solution=ccu, _name="curl2" );


            double di = normL2(elements(mesh), divv(g[i]) );
            double gn = normL2(boundaryfaces(mesh), trans(idv(g[i]))*N() );
            double curlgn = normL2(boundaryfaces(mesh), trans(idv(cu))*N() );
            double curl2gn = normL2(boundaryfaces(mesh), trans(idv(ccu))*N() );

            double erreur = normL2(elements(mesh), idv(cu) - sqrt(lambda[i])*idv(g[i]));
            double erreur2 = normL2(elements(mesh), idv(ccu) - lambda[i]*idv(g[i]) );
            // double erreur = normL2(elements(mesh), curlv(g[i]) - lambda[i]*idv(g[i]));
            // double erreur2 = normL2(elements(mesh), curlv(cu) - lambda[i]*lambda[i]*idv(g[i]) );

            // double avX = integrate(elements(mesh), trans(vec(cst(1),cst(0.),cst(0.)))*curlv(g[i]) ).evaluate()(0,0);
            // double avY = integrate(elements(mesh), trans(vec(cst(0.),cst(1),cst(0.)))*curlv(g[i]) ).evaluate()(0,0);
            // double avZ = integrate(elements(mesh), trans(vec(cst(0.),cst(0.),cst(1)))*curlv(g[i]) ).evaluate()(0,0);

            double normG = normL2( elements(mesh), idv(g[i]) );
            double normCurlG = normL2( elements(mesh), idv(cu) );
            double normCurl2G = normL2( elements(mesh), idv(ccu) );

            if ( Environment::worldComm().isMasterRank() ){
                std::cout << "lambda_" << i << " = " << lambda[i] << "    ||div|| = " << di << "    ||g.n|| = " << gn << "    ||curlg.n|| = " << curlgn/* << "    avX = " << avX << "    avY = " << avY << "    avZ = " << avZ*/ << "    ||curl2g.n|| = " << curl2gn << "     ||curl(g)-sqrt(l)*g|| = " << erreur << "    ||curl2(g)-l*g|| = " << erreur2/* << "     ||l*g||+||curlG||-2<curlG,l*g> = " << erreurb*/ << "    ||g|| = " << normG << "    ||curl(g)|| = " << normCurlG << "    ||curl2(g)|| = " << normCurl2G << std::endl;
            }
        }


        // gradu[i] = cu;
        // g0[i] = ccu;


        if ( Environment::worldComm().isMasterRank() && boption("needDebug") )
            std::cout << std::endl;


        i++;
        if(i>=nev)
            break;
    }


    if ( Environment::worldComm().isMasterRank() )
        s.close();
}


void
EigenProb::compute_decomp()
{
    // preparation for the decomposition

    auto vg0 = Vh->element();
    auto a2 = form2( _test=Vh, _trial=Vh );
    auto l2 = form1( _test=Vh );


    // H1xR
    auto Mlh = space_mltype::New( mesh );

    auto Psi = Mlh->element();
    auto psii = Psi.element<0>();
    auto nu = Psi.element<1>();
    auto a3 = form2( _test=Mlh, _trial=Mlh );
    auto l3 = form1( _test=Mlh );


    LOG(INFO) << "----- Mlh -----\n";
    LOG(INFO) << "[dof] number of dof: " << Mlh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Mlh->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(U): " << Mlh->template functionSpace<0>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Mlh->template functionSpace<0>()->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(P): " << Mlh->template functionSpace<1>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Mlh->template functionSpace<1>()->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Mlh -----\n";
        std::cout << "[dof] number of dof: " << Mlh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Mlh->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(U): " << Mlh->template functionSpace<0>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(U): " << Mlh->template functionSpace<0>()->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(P): " << Mlh->template functionSpace<1>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(P): " << Mlh->template functionSpace<1>()->nLocalDof() << "\n";
    }


    // left side of gi0
    a2 = integrate( _range=elements(mesh),
                    _expr=trans(curlt(vg0))*curl(vg0) );

    // left side of psi
    a3 = integrate( _range=elements(mesh),
                    _expr=inner(gradt(psii),grad(psii) )
                    + id(psii)*idt(nu) + idt(psii)*id(nu) );


    // left side for grad psi
    auto gv = Vh->element();
    auto c = form2( _trial=Vh, _test=Vh );
    auto f = form1( _test=Vh );

    if( boption( _name="needDebug") )
        c = integrate( _range=elements(mesh), _expr=inner(idt(gv),id(gv)) );


    for( int i = 0; i < nev; i++ ){
        // right side of gi0
        l2 = integrate( _range=elements(mesh),
                        _expr=lambda[i]*inner(idv(g[i]),id(vg0)) );

        a2+= on( _range=boundaryfaces(mesh),
                 _element=vg0, _rhs=l2, _expr=zero<3>() );

        a2.solve( _name="gi0", _rhs=l2, _solution=g0[i] );


        // right side of psi
        double meanPsi = doption(_name="meanPsi");

        l3 = integrate( _range=elements(mesh),
                        _expr=divv(g0[i])*id(psii)
                        + meanPsi*id(nu)
                        );

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
            double nGBis = normL2( elements(mesh), idv(g0[i])+idv(gradu[i]) );
            double nG0 = normL2( elements(mesh), idv(g0[i]) );
            double nPsi = normL2( elements(mesh), idv(gradu[i]) );

            if ( Environment::worldComm().isMasterRank() )
                std::cout << "||g-(g0+grad(psi)|| = " << erreurL2 << "    ||g0|| = " << nG0 << "    ||Grad(psi)|| = " << nPsi << "    ||g0+grad(psi)|| = " << nGBis << std::endl;
        }
    }

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
        std::cout << "Eigen values not found\ntry to launch with --computeEigen=true" << std::endl;
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
EigenProb::load_decomp()
{
    for( int i = 0; i < nev; i ++ ){
        std::string pathPsi = (boost::format("psi-%1%")%i).str();
        psi[i].load(_path=pathPsi);
    }
}

