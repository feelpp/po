#include <feel/feelvf/vf.hpp>

#include "spectralproblem.h"

using namespace Eigen;

SpectralProblem::SpectralProblem( mesh_ptrtype mesh ):super()
{
    this->mesh = mesh;
    Vh = space_vtype::New( mesh );
    Sh = space_stype::New( mesh );

    M = option( _name="solvereigen.nev" ).as<int>();
    j = MatrixXd(M,M);
    f = VectorXd(M);
    c = VectorXd(M);

    fa = Vh->element();
    alpha2 = option( _name="alpha2" ).as<std::string>();
    Re = option( _name="Re" ).as<double>(); // 2.*radius*speed/nu

}

void SpectralProblem::init( vector_vtype G, vector_stype P, std::vector<double> L, element_vtype A )
{
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
    Riak(M,M);
    for(int i = 0; i< M ; i++)
        for(int k = 0; k < M; k++)
            Riak(i,k) = integrate( _range=elements( mesh ),
                                   _expr=inner(cross( idv(g[i]), idv(a) ), idv(g[k]) ) ).evaluate()(0,0)*lambda[i];

}

void SpectralProblem::initRijk()
{
    Rijk(M,1);
    for(int k = 0; k < M; k++){
        Rijk(k) = MatrixXd(M,M);
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++)
                Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                                           _expr=inner(cross( idv(g[i]), idv(g[j]) ), idv(g[k]) ) ).evaluate()(0,0)*lambda[i];
    }
}

void SpectralProblem::initRfk()
{
    // f as an expression, h = project(Vh, f-idv(a)) ?
    Rfk = VectorXd(M);
    for(int k = 0; k < M; k++)
        Rfk(k) = integrate( _range=elements( mesh ),
                            _expr=inner( idv(fa), idv(g[k]) ) ).evaluate()(0,0);
}

void SpectralProblem::initRpk()
{
    auto vars = Symbols{ "x", "y", "z" };
    auto alpha_e = parse( this->alpha2, vars );
    auto alpha = expr( alpha_e, vars );

    Rpk = VectorXd(M);
    for(int k = 0; k < M; k++)
        Rpk(k) = integrate( _range=elements( mesh ),
                            _expr=alpha*idv(psi[k].template element<0>()) ).evaluate()(0,0);
}

void SpectralProblem::run()
{
    VectorXd dc = VectorXd::Ones(M);
    double tol = 1.e-6;
    HouseholderQR<MatrixXd> qr(M,M);

    int i;
    for(i = 0; i < 10 && dc.norm() > tol; i++){
        std::cout << "iteration : " << i << std::endl;
        f = c.cwiseProduct(lambda)/Re + c*Riak  - Rfk - Rpk/Re;
        for (int k = 0; k < M; k++)
            f(k) += c.transpose()*Rijk(k)*c;

        for (int k = 0; k < M; k++)
            j.row(k) = Riak.row(k) + c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
        j += lambda.asDiagonal();

        qr.compute(j);
        dc = qr.solve(-f);
        c += dc;
    }

    if(i==10)
        std::cout << "Newton does not converge\n";

}
