#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
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
    Re = r*s/n;
    u = Vh->element();
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

    initRiak();
    // initRijk();
    initRfk();
    initRpk();

    c = VectorXd::Ones(M);
}

void SpectralProblem::initRiak()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Riak -----" << std::endl;

    // [riak]
    Riak = MatrixXd(M,M);
    for(int i = 0; i< M ; i++)
        for(int k = 0; k < M; k++){
            Riak(k,i) = integrate( _range=elements( mesh ),
                                   _expr=inner(cross(idv(g[i]), idv(a)),idv(g[k])) ).evaluate()(0,0)*sqrt(lambda(i));
            // [riak]
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "Rki(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
        }
}

void SpectralProblem::initRijk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;

    // [rijk]
    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++){
        Rijk(k) = MatrixXd(M,M);
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                                          _expr=inner(cross( idv(g[i]), idv(g[j]) ), idv(g[k]) ) ).evaluate()(0,0)*sqrt(lambda[i]);
                // [rijk]
                if ( Environment::worldComm().isMasterRank() )
                    std::cout << "Rkij(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
            }
    }
}

void SpectralProblem::initRfk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rfk -----" << std::endl;

    auto ff = expr<3,1>(soption(_name="f"));
    // [rfk]
    Rfk = VectorXd(M);
    for(int k = 0; k < M; k++){
        Rfk(k) = integrate( _range=elements( mesh ),
                            _expr=trans(ff)*idv(g[k]) ).evaluate()(0,0);
    // [rfk]
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
    }
}

void SpectralProblem::initRpk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rpk -----" << std::endl;

    auto vars = Symbols{ "radius", "speed" };
    auto a2_e = parse( this->alpha2, vars );
    auto a2 = expr( a2_e, vars );
    a2.setParameterValues( {
            { "speed", doption( _name="speed" ) },
                { "radius", doption( _name="radius" ) } } );

    // [rpk]
    Rpk = VectorXd(M);
    for(int k = 0; k < M; k++){
        Rpk(k) = integrate( _range=markedfaces( mesh, 1 ),
                            _expr=a2*idv(psi[k]) ).evaluate()(0,0);
        Rpk(k) += integrate( _range=markedfaces( mesh, 2 ),
                             _expr=a2*idv(psi[k]) ).evaluate()(0,0);
        // [rpk]
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
    }
}

void SpectralProblem::run()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Start Spectral Problem -----" << std::endl;

    // [StokesA]
    MatrixXd A = MatrixXd(M,M);
    A = Riak;
    A += (lambda/Re).asDiagonal();
    // [StokesA]

    // [StokesB]
    VectorXd b = VectorXd(M);
    b= Rfk;
    // [StokesB]

    // [StokesSolve]
    HouseholderQR<MatrixXd> qr(M,M);
    qr.compute(A);
    c = qr.solve(b);
    // [StokesSolve]

    for(int i=0; i<M; i++){
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "c(" << i << ") = " << c(i) << std::endl;
        u += vf::project( _space=Vh, _range=elements(mesh),
                          _expr = c(i)*idv(g[i]) );
    }
    u += vf::project( _space=Vh, _range=elements(mesh), _expr = idv(a) );

    // [NSInit]
    // VectorXd dc = VectorXd::Matrix(M);
    // double tol = 1.e-6;
    // HouseholderQR<MatrixXd> qr(M,M);
    // [NSInit]

    // [NSNewton]
    // int i=0;
    // do{
    //     f = c.cwiseProduct(lambda)/Re + Riak*c - Rfk - Rpk/Re;
    //     for (int k = 0; k < M; k++)
    //         f(k) += c.transpose()*Rijk(k)*c;

    //     for (int k = 0; k < M; k++)
    //         j.row(k) = Riak.row(k) + c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
    //     j += lambda.asDiagonal();
    // [NSNewton]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "j = " << j << std::endl << "f = " << f << std::endl;

    // [NSSolve]
    //     qr.compute(j);
    //     dc = qr.solve(-f);
    //     c += dc;
    // [NSSolve]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << c << std::endl;

    // [NSNewtonEnd]
    //     i++;
    // } while(i < 10 && dc.norm() > tol);

    // if(i==10)
    //     std::cout << "Newton does not converge\n";
    // [NSNewtonEnd]
    // else{
    //     for( i = 0; i < M; i++){
    //         u += vf::project( _space=Vh, _range=elements(mesh),
    //                           _expr = c(i)*idv(g[i]) );
    //     }
    //     u += a;
    // }
}
