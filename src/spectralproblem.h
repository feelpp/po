/*
  lambda_k^2*c_k + Re*sum_i c_i*lambda_i int(g_ixa).g_k
                 + Re*sum_ij c_i*lambda_i*c_j*int(g_ixg_j).g_k
                 - Re/2*int h_a.g_k - int_dO alpha_2*psi_k
  for all 0 <= k <= M
*/
#ifndef __SPECTRALPROBLEM_H
#define __SPECTRALPROBLEM_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelalg/backend.hpp>

using namespace Feel;
using namespace Feel::vf;

class SpectralProblem : public Application
{
    typedef Application super;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Simplex<3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef std::vector<double> vector_1_double;
    typedef std::vector<std::vector<double> > vector_2_double;
    typedef std::vector<std::vector<std::vector<double> > > vector_3_double;

    backend_ptrtype M_backend;
    mesh_ptrtype mesh;
    vector_1_double c;
    vector_2_double m;
    vector_1_double lambda;
    vector_2_double Riak;
    vector_3_double Rijk;
    int n;

    SpectralProblem(mesh_ptrtype);
    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    void run();
};

#endif
