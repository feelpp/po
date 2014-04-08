/*
  u = grad(psi)
  div(u) = 0
  u.n = alpha

  -int_O u*v + grad(p)*v + grad(q)*u = int_pO alpha*q
 */
#ifndef __DARCY_H
#define __DARCY_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

using namespace Feel;
using namespace Feel::vf;

class Darcy : public Application
{
    typedef Application super;
 public:
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<RaviartThomas<0>, Lagrange<1, Scalar> > prod_basis_type;
    // typedef bases<RaviartThomas<0>, Lagrange<1, Scalar>, Lagrange<0, Scalar> > prod_basis_type;
    typedef FunctionSpace<mesh_type, prod_basis_type > space_type;
    typedef space_type::element_type element_type;

    Darcy(mesh_ptrtype, std::string);
    void run();

    element_type U;

 private:
    mesh_ptrtype mesh;
    std::string g_s;
};

#endif
