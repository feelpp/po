#include <feel/feelvf/vf.hpp>

#include "spectralproblem.h"

using namespace Eigen;

SpectralProblem::SpectralProblem( mesh_ptrtype mesh ):super()
{
    this->mesh = mesh;
    Vh = space_vtype::New( mesh );
    Sh = space_stype::New( mesh );

    M = option( _name="solvereigen.nev" ).as<int>();
    j = MatrixXd::Matrix(M,M);
    f = VectorXd::Matrix(M);
    c = VectorXd::Matrix(M);

    fa = Vh->element();
    alpha2 = option( _name="alpha2" ).as<std::string>();
    Re = option( _name="Re" ).as<double>(); // 2.*radius*speed/nu

}

void SpectralProblem::init( vector_vtype G, vector_stype P, std::vector<double> L, element_vtype A )
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Initialization Spectral Problem -----" << std::endl;

    lambda = VectorXd::Matrix(M);
    for(int i = 0; i < M; i++)
        lambda(i) = L[i];
    g = G;
    psi = P;
    a = A;

    initRiak();
    initRijk();
    initRfk();
    initRpk();

    c = VectorXd::Ones(M);
}

void SpectralProblem::initRiak()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Riak -----" << std::endl;
    Riak = MatrixXd::Matrix(M,M);
    for(int i = 0; i< M ; i++)
        for(int k = 0; k < M; k++){
            Riak(k,i) = integrate( _range=elements( mesh ),
                                   _expr=inner(cross(idv(g[i]), idv(a)),idv(g[k])) ).evaluate()(0,0)*sqrt(lambda[i]);
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "Rki(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
        }

}

void SpectralProblem::initRijk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;
    Rijk = Matrix<MatrixXd, Dynamic, 1>::Matrix(M,1);
    for(int k = 0; k < M; k++){
        Rijk(k) = MatrixXd::Matrix(M,M);
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                                          _expr=inner(cross( idv(g[i]), idv(g[j]) ), idv(g[k]) ) ).evaluate()(0,0)*sqrt(lambda[i]);
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "Rkij(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
        }
    }
}

void SpectralProblem::initRfk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rfk -----" << std::endl;
    // f as an expression, h = project(Vh, f-idv(a)) ?
    Rfk = VectorXd::Matrix(M);
    for(int k = 0; k < M; k++)
        Rfk(k) = integrate( _range=elements( mesh ),
                            _expr=inner( idv(fa), idv(g[k]) ) ).evaluate()(0,0);
}

void SpectralProblem::initRpk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rpk -----" << std::endl;
    auto vars = Symbols{ "x", "y", "z" };
    auto alpha_e = parse( this->alpha2, vars );
    auto alpha = expr( alpha_e, vars );

    Rpk = VectorXd::Matrix(M);
    for(int k = 0; k < M; k++)
        Rpk(k) = 0;
        // Rpk(k) = integrate( _range=elements( mesh ),
        //                     _expr=alpha*idv(psi[k].template element<0>()) ).evaluate()(0,0);
}

void SpectralProblem::run()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Start Spectral Problem -----" << std::endl;    
    
    VectorXd dc = VectorXd::Matrix(M);
    double tol = 1.e-6;
    HouseholderQR<MatrixXd> qr(M,M);

    int i=0;
    do{
        f = c.cwiseProduct(lambda)/Re + Riak*c - Rfk - Rpk/Re;
        for (int k = 0; k < M; k++)
            f(k) += c.transpose()*Rijk(k)*c;

        for (int k = 0; k < M; k++)
            j.row(k) = Riak.row(k) + c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
        j += lambda.asDiagonal();

        if ( Environment::worldComm().isMasterRank() )
            std::cout << "j = " << j << std::endl << "f = " << f << std::endl;

        qr.compute(j);
        dc = qr.solve(-f);
        c += dc;

        if ( Environment::worldComm().isMasterRank() )
            std::cout << "iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

        if ( Environment::worldComm().isMasterRank() )
            std::cout << c << std::endl;

        i++;
    } while(i < 10 && dc.norm() > tol);

    if(i==10)
        std::cout << "Newton does not converge\n";
    else{
        u = Vh->element();
        for( i = 0; i < M; i++){
            u += vf::project( _space=Vh, _range=elements(mesh),
                              _expr = c(i)*idv(g[i]) );
        }
        u += a;
    }

}
