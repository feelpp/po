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

using namespace Feel;
using namespace Feel::vf;

class Psi0 : public Application
{
    typedef Application super;
 public:
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // [space]
    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar>, Lagrange<0, Scalar> > > mlSpace_type;
    // [space]

    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Vectorial> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    element_type gradu;

    Psi0(mesh_ptrtype, std::string);
    void run();

 private:
    mesh_ptrtype mesh;
    std::string g_s;
    space_ptrtype Xh;
};

#endif
