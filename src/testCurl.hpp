#ifndef TESTCURL_H
#define TESTCURL_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>

using namespace Feel;


class TestCurl : public Application
{
    typedef Application super;

    static const uint16_type Dim = 3;
    static const uint16_type Order = 3;

public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

//typedef bases<Nedelec<0,NedelecKind::NED1> > basis_nedtype;
    typedef bases<Lagrange<Order, Vectorial> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    std::vector<element_type> test;

    TestCurl( mesh_ptrtype );
    void run();
 private:

    mesh_ptrtype mesh;
    space_ptrtype Vh;
};

#endif
