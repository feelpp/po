/** \file eigenprob.hpp
    \brief Header file for the class EigenProb
*/

#ifndef EIGENPROB_H
#define EIGENPROB_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>


using namespace Feel;
using namespace Feel::vf;

/** \brief Resolves an eigen problem.

    \f{align*}
    \rott \mathbf{u} &= \Lambda\mathbf{u}\\
    \div\mathbf{u} &= 0\\
    \mathbf{u}\cdot\mathbf{n}\restr &= 0\\
    \rot \mathbf{u}\cdot\mathbf{n}\restr &= 0\\
    \rott \mathbf{u}\cdot\mathbf{n}\restr &= 0\\
    \f}

    which leads to the following weak form :
    \f{align*}
    \int_\Omega (\rot\mathbf{u})\cdot(\rot\bm{\varphi}) &+ \int_\Omega\bm{\varphi}\grad p + \int_\Omega \mathbf{u}\grad q \\
&+ \alpha\int_\Omega \div\mathbf{g}\div\bm{\varphi} \\
&+ \beta\int_{\partial\Omega}(\mathbf{g}\cdot\mathbf{n})(\bm{\varphi}\cdot\mathbf{n}) \\
&+ \gamma\int_{\partial\Omega}(\rot\mathbf{u}\cdot\mathbf{n})(\rot\bm{\varphi}\cdot\mathbf{n})  = \Lambda\int_\Omega \mathbf{g}\cdot\bm{\varphi}
    \f}
 */
class EigenProb : public Application
{
  /// Inherits Application
    typedef Application super;

  /// Dimension of the mesh
    static const uint16_type Dim = 3;
  /// Order of the function space
    static const uint16_type Order = 2;

  /// Simplex of dimension Dim
    typedef Simplex<Dim> convex_type;
  /// Mesh of type convex_type
    typedef Mesh<convex_type> mesh_type;
  /// Pointer on the mesh
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //typedef bases<Nedelec<0,NedelecKind::NED1>, Lagrange<1, Scalar> > basis_ptype;
    // [space]
  /// Basis of \f$[P_{Order}]^3\times P_1\f$
    typedef bases<Lagrange<Order, Vectorial>, Lagrange<1, Scalar> > basis_ptype;
    // [space]
  /// Function space \f$[P_{Order}]^3\times P_1\f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_ptype > space_ptype;

  /// Basis of \f$ P_{Order}\times P_0\f$
    typedef bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > basis_mltype;
  /// Function space \f$ P_{Order}\times P_0\f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_mltype > space_mltype;


    //typedef bases<Nedelec<0,NedelecKind::NED1> > basis_vtype;
  /// Basis of \f$[P_{Order}]^3\f$
    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
  /// Function space \f$[P_{Order}]^3\f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
  /// Element of \f$[P_{Order}]^3\f$
    typedef typename space_vtype::element_type element_vtype;
  /// Pointer on the function space \f$[P_{Order}]^3\f$
    typedef boost::shared_ptr<space_vtype> space_vptrtype;

  /// Basis of \f$ P_{Order}\f$
    typedef bases<Lagrange<Order, Scalar> > basis_stype;
  /// Function space \f$ P_{Order}\f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
  ///Element of \f$ P_{Order}\f$
    typedef typename space_stype::element_type element_stype;
  /// Pointer on the function space \f$ P_{Order}\f$
    typedef boost::shared_ptr<space_stype> space_sptrtype;

 public:
  /// Array of element of \f$[P_{Order}]^3\f$ containing the eigen functions
    std::vector<element_vtype> g;
  /// Array of double containing the eigen values
    std::vector<double> lambda;
  /// Array of element of \f$ P_{Order}\f$ containing the \f$\psi\f$ of the decomposition of the eigen functions
    std::vector<element_stype> psi;

    // debug
  /// Debug
    std::vector<element_vtype> g0;
  /// Debug
    std::vector<element_vtype> gradu;
  /// Debug
    std::vector<element_vtype> modebis;

  /** \brief Constructs an object of type EigenProb

      \param mesh the pointer on the mesh
  */
    EigenProb( mesh_ptrtype mesh);

  /** \brief Find the eigen modes

      If the option computeEigen is set to true, launch compute_eigens(), else launch load_eigens()\n
      If the option needDecomp and computeDecomp are set to true, launch compute_decomp(), else if computeDecomp is set to false, launch load_decomp()
  */
    void run();

 private:
  /** \brief Compute the eigen modes

      Resolves the eigen problem.\n
      Depend on the options useCurl, usePresDiv, usePresGrad, useDiric, bccurln, bcn, divdiv
  */
    void compute_eigens();

  /** \brief Compute the decomposition of the eigen functions

      Resolves the two system in order to write \f$\mathbf{g}=\mathbf{g_0}+\grad \psi\f$\n
      Depend on the option meanPsi
   */
    void compute_decomp();

  /** \brief Load the eigen modes.

      Load the eigen modes from the disk, the partition mesh and the number of eigen modes need to be identical
  */
    void load_eigens();

  /** \brief Load the \f$\psi\f$ of the decomposition

      Load the \f$\psi\f$ of the decomposition from the disk, the partition mesh and the number of eigen modes need to be identical
   */
    void load_decomp();

  /// The mesh used for the application.
    mesh_ptrtype mesh;
  /// The vectorial function space used.
    space_vptrtype Vh;
  /// The scalar function space used.
    space_sptrtype Sh;
  /// The number of eigen modes needed.
    int nev;
};

#endif
