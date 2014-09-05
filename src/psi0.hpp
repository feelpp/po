/** \file psi0.hpp
    \brief Header file for the class Psi0
 */

#ifndef __PSI0_H
#define __PSI0_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>


using namespace Feel;
using namespace Feel::vf;

/// Resolves a Poisson problem
/**
  \f{align*}
  -\laplace u &= 0 \\
  \grad u \cdot \mathbf{n}\restr &= g
  \f}

  which, with the use of Lagrange multipliers for the constraint \f$ \int_\Omega u = 0 \f$, leads to the variationnal formulation :
  \f{equation*}
  \int_\Omega \grad u \cdot \grad v + u \nu + v\lambda = \int_{\partial\Omega} gv
  \f}
 */
class Psi0 : public Application
{
    typedef Application super;

    static const uint16_type Order = 2;

 public:
    /// Mesh of dimension 3
    typedef Mesh<Simplex<3> > mesh_type;
    /// Pointer on the mesh
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // [space]
    /// Function space \f$P_2\times P_0\f$ on the mesh
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > > mlSpace_type;
    // [space]

    /// Function space \f$[P_2]^3\f$ on the mesh
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Vectorial> > > space_type;
    /// Pointer on the function space
    typedef boost::shared_ptr<space_type> space_ptrtype;
    /// Element of the function space
    typedef space_type::element_type element_type;
    /// Pointer on the element
    typedef boost::shared_ptr<element_type> element_ptrtype;

    /// \f$ \grad u\f$
    element_type gradu;

  /** \brief Initializer
      \param mesh the pointer on the mesh
      \param g the the boundary condition g
  */
    Psi0(mesh_ptrtype mesh, std::string g);

  /** \brief Find \f$\grad u\f$
      
      If the option computeP0 is set to true, launch compute_psi0()
      else, launch load_psi0()
  */
    void run();

  /** \brief Compute \f$ \grad u\f$

      Resolves the problem on the mesh with the boundary condition g
  */
    void compute_psi0();


  /** \brief Load \f$\grad u \f$

      Load \f$\grad u \f$ from the disk, the partition mesh needs to be identical
  */
    void load_psi0();

 private:
    mesh_ptrtype mesh;
    std::string g_s;
    space_ptrtype Xh;
};

#endif
