/*
  (1)
  curl2(g) = lambda*g
  g.n = curl(g).n = curl2(g).n = 0

  int_O curl(g)*curl(v) = lambda*int_O g*v

  (2)
  grad(div(g0))-laplace(g0)=-laplace(g)
  g0 = 0 on pO

  -int_O div(g0)*div(v) + int_0 grad(g0)*grad(v) = int_O grad(g)*grad(v)

  (3)
  -laplace(psi) = div(g0)
  grad(psi).n = 0

  int_0 grad(psi)*grad(v) = int_0 div(g0)*v
 */
#ifndef __EIGEN_CURL_H
#define __EIGEN_CURL_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>

using namespace Feel;
using namespace Feel::vf;

class Eigen_Curl : public Application
{
    typedef Application super;
    static const uint16_type Order = 2;
    static const uint16_type Dim = 3;

 public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order, Scalar>, Lagrange<Order, Scalar>, Lagrange<Order, Scalar> > sBasis_type;
    typedef FunctionSpace<mesh_type, sBasis_type> sSpace_type;
    typedef sSpace_type::element_type sElement_type;
    typedef bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > mlBasis_type;
    typedef FunctionSpace<mesh_type, mlBasis_type > mlSpace_type;
    typedef boost::shared_ptr<mlSpace_type> mlSpace_ptrtype;
    typedef mlSpace_type::element_type mlElement_type;
    typedef bases<Lagrange<Order, Vectorial> > vBasis_type;
    typedef FunctionSpace<mesh_type, vBasis_type > vSpace_type;
    typedef vSpace_type::element_type vElement_type;
    typedef boost::shared_ptr<vSpace_type> vSpace_ptrtype;

    std::vector<vElement_type> g;
    std::vector<double> lambda;
    std::vector<mlElement_type> psi;

    Eigen_Curl(bool, mesh_ptrtype);
    void run();
    void load_eigens();
    void decomp();

 private:
    mesh_ptrtype mesh;
    vSpace_ptrtype Vh;
    mlSpace_ptrtype Mlh;
    int nev;
    int ncv;
};

#endif
