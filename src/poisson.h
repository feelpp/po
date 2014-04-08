/*
  -laplace(u) = 0
  grad(u).n = g

  int_O grad(u)*grad(v) + u*nu + v*lambda = int_pO g*v
  use of Lagrange multipliers int_0 u = 0
 */
#ifndef __POISSON_H
#define __POISSON_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>

using namespace Feel;
using namespace Feel::vf;

class Poisson : public Application
{
    typedef Application super;
 public:
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Scalar>, Lagrange<0, Scalar> > > mlSpace_type;
    typedef mlSpace_type::element_type mlElement_type;
    //typedef typename mlElement_type:: sub_element<0>::type element_0_type;
    typedef boost::shared_ptr<mlElement_type> mlElement_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial> > > space_type;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    mlElement_type U;
    element_type gradu;

    Poisson(mesh_ptrtype, std::string);
    void run();

 private:
    mesh_ptrtype mesh;
    std::string g_s;
};

#endif
