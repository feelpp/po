#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/mpi/timer.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/cross.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/ginac.hpp>

#include "spectralproblem.hpp"


using namespace Eigen;


SpectralProblem::SpectralProblem( mesh_ptrtype mesh ):super()
{
    this->mesh = mesh;
    Vh = space_vtype::New( mesh );
    Sh = space_stype::New( mesh );

    M = ioption( _name="solvereigen.nev" );
    j = MatrixXd(M,M);
    f = VectorXd(M);
    c = VectorXd(M);

    alpha2 = soption( _name="alpha2" );
    double r = doption( _name="radius" );
    double s = doption( _name="speed" );
    double n = doption( _name="nu" );
    Re = 2*r*s/n;
    u = Vh->element();
    v = Vh->element();
}

void SpectralProblem::init( vector_vtype G, vector_stype P, std::vector<double> L, element_vtype A )
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Initialization Spectral Problem -----" << std::endl;
        std::cout << "----- Re = " << Re << " -----" << std::endl;
        std::cout << "----- alpha2 = " << alpha2 << " -----" << std::endl;
    }


    // [lambda]
    lambda = VectorXd(M);

    for(int i = 0; i < M; i++)
        lambda(i) = L[i];
    // [lambda]

    g = G;
    psi = P;
    a = A;

    boost::mpi::timer t;

    initRiak();

    LOG(INFO) << "Riak = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Riak = " << t.elapsed() << " sec" << std::endl;
    t.restart();

    // initRijk();

    LOG(INFO) << "Rijk = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Rijk = " << t.elapsed() << " sec" << std::endl;
    t.restart();

    initRfk();

    LOG(INFO) << "Rfk = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Rfk = " << t.elapsed() << " sec" << std::endl;
    t.restart();

    initRpk();

    LOG(INFO) << "Rpk = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Rpk = " << t.elapsed() << " sec" << std::endl;

    c = VectorXd::Ones(M);
}

void SpectralProblem::initRiak()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Riak -----" << std::endl;

    // [riakInit]
    Riak = MatrixXd(M,M);
    // [riakInit]
    if( boption("computeRiak") ){
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("riak", std::fstream::out);
        // [riakComp]
        for(int i = 0; i< M ; i++){
            for(int k = 0; k < M; k++){
                Riak(k,i) = integrate( _range=elements( mesh ),
                                       _expr=inner( cross( idv(g[i]),idv(a) ), idv(g[k])) ).evaluate()(0,0) * sqrt(lambda(i));
                // [riakComp]
                if ( Environment::worldComm().isMasterRank() ){
                    std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
                    s << Riak(k,i) << std::endl;
                }
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else {
        std::fstream s;
        s.open ("riak", std::fstream::in);
        if( !s.is_open() ){
            std::cout << "Riak not found\ntry to launch with --computeRiak=true" << std::endl;
            exit(0);
        }
        for(int i = 0; i< M ; i++)
            for(int k = 0; k < M; k++)
                s >> Riak(k,i);
        s.close();
    }
}

void SpectralProblem::initRijk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;
    // [rijkInit]
    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);
    // [rijkInit]
    if( boption("computeRijk") ){
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rijk", std::fstream::out);

        // [rijkCompute]
        for(int k = 0; k < M; k++){
            for(int i = 0; i < M; i++){
                for(int j = 0; j < M; j++){
                    Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                                              _expr=inner( cross( idv(g[i]),idv(g[j]) ), idv(g[k]) ) ).evaluate()(0,0) * sqrt(lambda[i]);
                    // [rijkCompute]
                    if ( Environment::worldComm().isMasterRank() ){
                        std::cout << "Rijk(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
                        s << Rijk(k,i) << std::endl;
                    }
                }
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else {
        std::fstream s;
        s.open ("rijk", std::fstream::in);
        if( !s.is_open() ){
            std::cout << "Rijk not found\ntry to launch with --computeRijk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
            for(int i = 0; i < M; i++)
                for(int j = 0; j < M; j++)
                    s >> Rijk(k)(i,j);
        s.close();
    }
}

void SpectralProblem::initRfk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rfk -----" << std::endl;

    // [rfkInit]
    Rfk = VectorXd(M);
    // [rfkInit]
    if( boption("computeRfk") ){
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rfk", std::fstream::out);

        // [rfkComp]
        auto ff = expr<3,1>(soption(_name="f"));
        for(int k = 0; k < M; k++){
            Rfk(k) = integrate( _range=elements( mesh ),
                                _expr=trans(ff)*idv(g[k]) ).evaluate()(0,0);
            // [rfkComp]
            if ( Environment::worldComm().isMasterRank() ){
                std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
                s << Rfk(k) << std::endl;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else {
        std::fstream s;
        s.open ("rfk", std::fstream::in);
        if( !s.is_open() ){
            std::cout << "Rfk not found\ntry to launch with --computeRfk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
            s >> Rfk(k);
        s.close();
    }
}

void SpectralProblem::initRpk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rpk -----" << std::endl;

    // auto vars = Symbols{ "radius", "speed" };
    // auto a2_e = parse( this->alpha2, vars );
    // auto a2 = expr( a2_e, vars );
    // a2.setParameterValues( {
    //         { "speed", doption( _name="speed" ) },
    //             { "radius", doption( _name="radius" ) } } );

    // [rpkInit]
    Rpk = VectorXd(M);
    // [rpkInit]
    if( boption("computeRpk") ){
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rpk", std::fstream::out);

        // [rpkComp]
        auto a2 = expr( alpha2 );
        for(int k = 0; k < M; k++){
            Rpk(k) = integrate( _range=markedfaces( mesh, 1 ),
                                _expr=a2*idv(psi[k]) ).evaluate()(0,0);
            Rpk(k) -= integrate( _range=markedfaces( mesh, 2 ),
                                 _expr=a2*idv(psi[k]) ).evaluate()(0,0);
            // [rpkComp]

            //Rpk(k) = integrate( _range=boundaryfaces(mesh),
            //                    _expr=a2*idv(psi[k]) ).evaluate()(0,0);

            if ( Environment::worldComm().isMasterRank() ){
                std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
                s << Rpk(k) << std::endl;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else {
        std::fstream s;
        s.open ("rpk", std::fstream::in);
        if( !s.is_open() ){
            std::cout << "Rpk not found\ntry to launch with --computeRpk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
            s >> Rpk(k);
        s.close();
    }
}

void SpectralProblem::run()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Start Spectral Problem -----" << std::endl;

    boost::mpi::timer t;

    // [StokesA]
    MatrixXd A = MatrixXd(M,M);
    A = Riak;
    A += (lambda/Re).asDiagonal();
    // [StokesA]

    // [StokesB]
    VectorXd b = VectorXd(M);
    b = Rfk + Rpk;
    // [StokesB]

    // [StokesSolve]
    HouseholderQR<MatrixXd> qr(M,M);
    qr.compute(A);
    c = qr.solve(b);
    // [StokesSolve]

    LOG(INFO) << "stokes = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "stokes = " << t.elapsed() << " sec" << std::endl;


    for(int i=0; i<M; i++){
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "c(" << i << ") = " << c(i) << std::endl;
        u += vf::project( _space=Vh, _range=elements(mesh),
                          _expr = c(i)*idv(g[i]) );
    }
    v = vf::project( _space=Vh, _range=elements(mesh), _expr = idv(a)+idv(u) );


    // // [NSInit]
    // VectorXd dc = VectorXd::Matrix(M);
    // double tol = 1.e-6;
    // // [NSInit]
    // // [NSSys1]
    // HouseholderQR<MatrixXd> qr(M,M);
    // // [NSSys1]

    // int i=0;
    // do{
    //     // [NSMatF]
    //     f = c.cwiseProduct(lambda)/Re + Riak*c - Rfk;
    //     for (int k = 0; k < M; k++)
    //         f(k) += c.transpose()*Rijk(k)*c;
    //     // [NSMatF]

    //     // [NSMatJ]
    //     for (int k = 0; k < M; k++)
    //         j.row(k) = c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
    //     j += Riak;
    //     j += lambda.asDiagonal();
    //     // [NSMatJ]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "j = " << j << std::endl << "f = " << f << std::endl;

    //     // [NSSys2]
    //     qr.compute(j);
    //     dc = qr.solve(-f);
    //     // [NSSys2]
    //     // [NSAdd]
    //     c += dc;
    //     // [NSAdd]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << c << std::endl;

    //     i++;
    // } while(i < 10 && dc.norm() > tol);

    // if(i==10)
    //     std::cout << "Newton does not converge\n";

    // else{
    //     for( i = 0; i < M; i++){
    //         u += vf::project( _space=Vh, _range=elements(mesh),
    //                           _expr = c(i)*idv(g[i]) );
    //     }
    //     u += a;
    // }
}
