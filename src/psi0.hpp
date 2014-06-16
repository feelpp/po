/*
  -laplace(u) = 0
  grad(u).n = g

  int_O grad(u)*grad(v) + u*nu + v*lambda = int_pO g*v
  use of Lagrange multipliers int_0 u = 0
 */
#ifndef __PSI0_H
#define __PSI0_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

using namespace Feel;
using namespace Feel::vf;

class Psi0 : public Application
{
    typedef Application super;
 public:
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<RaviartThomas<0>, Lagrange<2, Scalar> > prod_basis_type;
    // typedef bases<RaviartThomas<0>, Lagrange<2, Scalar>, Lagrange<0, Scalar> > prod_basis_type;
    typedef FunctionSpace<mesh_type, prod_basis_type > space_type;
    typedef space_type::element_type element_type;

    Psi0(mesh_ptrtype, std::string);
    void run();

    element_type U;

 private:
    mesh_ptrtype mesh;
    std::string g_s;

};

#endif
