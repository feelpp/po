#include "spectralproblem.h"

SpectralProblem::SpectralProblem(mesh_ptrtype mesh):super(),M_backend( backend_type::build( this->vm() ) )
{
    this->mesh = mesh;
    n = option(_name="nev").as<int>();
    c = vector_1_double( n, 0);
    m = vector_2_double( n, vector_1_double(n, 0));
    lambda = vector_1_double( n, 0);
    // auto g = std::vector<element_type>(n, element);
    Riak = vector_2_double(n, vector_1_double(n, 0));
    Rijk = vector_3_double(n, vector_2_double(n, vector_1_double(n, 0) ) );
    // M_backend->nlSolver()->residual =
    //     boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    // M_backend->nlSolver()->jacobian =
    //     boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );
}

void SpectralProblem::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
}
void SpectralProblem::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
}

void SpectralProblem::run()
{
    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J){
        double res;
        for(int k=0; k < n; k++){
            for(int i=0; i < n; i++){
                res = 0;
                if(i==k)
                    res += lambda[i];
                res += lambda[i]*Riak[i][k];
                for(int j=0; j < n; j++){
                    res += lambda[i]*c[j]*Rijk[i][j][k];
                    res += c[j]*lambda[j]*Rijk[j][i][k];
                }
                m[k][i]=res;
            }
        }
    };

    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R){
    };

    backend()->nlSolver()->residual =Residual;
    backend()->nlSolver()->jacobian =Jacobian;

}
