#ifndef TESTCURL_H
#define TESTCURL_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>

using namespace Feel;


class TestCurl : public Application
{
    typedef Application super;

    static const uint16_type Dim = 3;
    static const uint16_type Order = 2;
 public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
    typedef typename space_vtype::element_type element_vtype;
    typedef boost::shared_ptr<space_vtype> space_vptrtype;

    std::vector<element_vtype> test;

    TestCurl( mesh_ptrtype );
    void run();
 private:

    mesh_ptrtype mesh;
    space_vptrtype Vh;
};

#endif
