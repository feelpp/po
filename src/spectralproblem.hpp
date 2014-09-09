/** \file spectralproblem.hpp
    \brief Header file for the class SpectralProblem
*/

#ifndef __SPECTRALPROBLEM_H
#define __SPECTRALPROBLEM_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <Eigen/Dense>


using namespace Feel;
using namespace Feel::vf;
using namespace Eigen;

/** \brief Resolves the spectral problem

    \f{align*}
    \frac{\partial c_k}{\partial t} + \lambda_k^2 c_k &+ Re\sum_i c_i\lambda_i \int (g_i\times a)\cdot g_k\\
                 &+ Re\sum_{ij} c_i\lambda_i c_j\int (g_i \times g_j)\cdot g_k\\
                 &- \frac{Re}{2}\int h_a \cdot g_k - \int_{\partial\Omega} \alpha_2 \psi_k\\
  \forall 0 \leq k \leq M &
  \f}
*/
class SpectralProblem : public Application
{
  /// Inherits Application
    typedef Application super;
  /// Order of the function space
    static const uint16_type Order = 2;

  /// Simplex of dimension 3
    typedef Simplex<3> convex_type;
  /// Mesh of the type convex_type
    typedef Mesh<convex_type> mesh_type;
  /// Pointer on the mesh
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

  /// Basis of \f$ [P_{Order}]^3 \f$
    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
  /// Function space \f$ [P_{Order}]^3 \f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
  /// Pointer on the function space \f$ [P_{Order}]^3 \f$
    typedef boost::shared_ptr<space_vtype> space_ptrvtype;
  /// Element of \f$ [P_{Order}]^3 \f$
    typedef space_vtype::element_type element_vtype;

    //typedef bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > basis_stype; // ML
  /// Basis of \f$ P_{Order} \f$
    typedef bases<Lagrange<Order, Scalar> > basis_stype;
  /// Function space \f$ P_{Order} \f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
  /// Pointer on the function space \f$ P_{Order} \f$
    typedef boost::shared_ptr<space_stype> space_ptrstype;
  /// Element of \f$ P_{Order} \f$
    typedef space_stype::element_type element_stype;

  /// Array of element of \f$ [P_{Order}]^3 \f$
    typedef std::vector<element_vtype> vector_vtype;
  /// Array of element of \f$ P_{Order} \f$
    typedef std::vector<element_stype> vector_stype;

  /// The mesh used for the application.
    mesh_ptrtype mesh;
  /// The vectorial function space used.
    space_ptrvtype Vh;
  /// The scalar function space used.
    space_ptrstype Sh;

  /// The number of eigen modes
    int M;
  /// The matrix containing the jacobian
    MatrixXd j;
  /// The vecotr containing the second member
    VectorXd f;

  /// An array containing the eigen modes
    vector_vtype g;
  /// A vector containing the eigen values
    VectorXd lambda;
  /// An array containing the \f$ \psi \f$ of the decomposition
    vector_stype psi;
  /// An element of \f$ P_{Order} \f$ containing \f$\mathbf{a}\f$
    element_vtype a;
  /// A string containing the right side \f$ \mathbf{f}\f$
    std::string fa_s;
  /// A string containing the boundary condition \f$\alpha_2\f$
    std::string alpha2;
  /// The Reynolds
    double Re;

    // [ri]
  /// The matrix containing \f$ R_{iak}\f$
    MatrixXd Riak;
  /// The vector of matrix containing \f$ R_{ijk}\f$
    Matrix<MatrixXd, Dynamic, 1 > Rijk;
  /// The vector containing \f$ R_{fk}\f$
    VectorXd Rfk;
  /// The vector containing \f$ R_{pk} \f$
    VectorXd Rpk;
    // [ri]

/** \brief Initialize \f$R_{iak}\f$

    If the option computeRiak is set to true, compute \f$ \int_\Omega (\mathbf{g}_i\times\mathbf{a})\cdot \mathbf{g}_k\quad \forall i,k \f$\n
    else, load them from the disk
*/
    void initRiak();

/** \brief Initialize \f$R_{ijk}\f$

    If the option computeRijk is set to true, compute \f$ \int_\Omega (\mathbf{g}_i\times\mathbf{g}_j)\cdot \mathbf{g}_k\quad \forall i,j,k \f$\n
    else, load them from the disk
*/
    void initRijk();

/** \brief Initialize \f$R_{fk}\f$

    If the option computeRfk is set to true, compute \f$ \int_\Omega \mathbf{f}\cdot \mathbf{g}_k\quad \forall k \f$\n
    else, load them from the disk
*/
    void initRfk();

/** \brief Initialize \f$R_{pk}\f$

    If the option computeRpk is set to true, compute \f$ \int_{\partial\Omega} \alpha_2\psi_k \quad \forall k \f$\n
    else, load them from the disk
*/
    void initRpk();

 public:
  /** \brief Constructs an object of type SpectralProblem

      Just claim the memory, not fit to use
      \param mesh the pointer on the mesh
   */
    SpectralProblem( mesh_ptrtype mesh );

  /** \brief Initialize the object with the eigen modes and \f$ \mathbf{a}\f$

      Initialialize the eigen functions, the eigen values, the \f$\psi\f$ of the decomposition and \f$\mathbf{a}\f$\n
      Call initRiak(), initRijk(), initRfk() and initRpk()
      \param g the array containing the eigen functions
      \param psi the array containing the \f$\psi\f$ of the decomposition
      \param l the array containg the eigen values
      \param a the element of \f$[P_{Order}]^3\f$ containing \f$\mathbf{a}\f$
   */
    void init( vector_vtype g, vector_stype psi, std::vector<double> l, element_vtype a );

/** \brief Resolves the spectral problem

    Use a Newton method to resolves a Navier-Stokes problem or a linear equation in the case of Stokes problem
*/
    void run();

  /// The coefficients \f$c_i\f$
    VectorXd c;
  /// \f$\mathbf{u} = \sum_i c_i\mathbf{g}_i\f$
    element_vtype u;
  /// \f$\mathbf{v} = \mathbf{a}+\mathbf{u}\f$
    element_vtype v;
};

#endif
