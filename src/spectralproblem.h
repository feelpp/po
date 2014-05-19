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
#include <Eigen/Dense>

using namespace Feel;
using namespace Feel::vf;
using namespace Eigen;

class SpectralProblem : public Application
{
    typedef Application super;

    typedef Simplex<3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<2, Vectorial> > basis_vtype;
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
    typedef boost::shared_ptr<space_vtype> space_ptrvtype;
    typedef space_vtype::element_type element_vtype;

    //typedef bases<Lagrange<2, Scalar>, Lagrange<0, Scalar> > basis_stype; // ML
    typedef bases<Lagrange<2, Scalar> > basis_stype;
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
    typedef boost::shared_ptr<space_stype> space_ptrstype;
    typedef space_stype::element_type element_stype;

    typedef std::vector<element_vtype> vector_vtype;
    typedef std::vector<element_stype> vector_stype;

    mesh_ptrtype mesh;
    space_ptrvtype Vh;
    space_ptrstype Sh;

    int M;
    MatrixXd j;
    VectorXd f;

    vector_vtype g;
    VectorXd lambda;
    vector_stype psi;
    element_vtype a;
    std::string fa_s;
    std::string alpha2;
    double Re;

    MatrixXd Riak;
    Matrix<MatrixXd, Dynamic, 1 > Rijk;
    VectorXd Rfk;
    VectorXd Rpk;

    void initRiak();
    void initRijk();
    void initRfk();
    void initRpk();

 public:
    SpectralProblem( mesh_ptrtype );
    void init( vector_vtype, vector_stype, std::vector<double>, element_vtype );
    void run();

    VectorXd c;
    element_vtype u;
    element_vtype v;
};

#endif
