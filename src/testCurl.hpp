/** \file testCurl.hpp
    \brief Header file for the class TestCurl
*/

#ifndef TESTCURL_H
#define TESTCURL_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>

using namespace Feel;

/// Test the operator curl
class TestCurl : public Application
{
    /// Inherits Application
    typedef Application super;

    /// Dimension of the domain
    static const uint16_type Dim = 3;
    /// Order of the function space
    static const uint16_type Order = 3;

public:
    /// Simplex of dimension Dim
    typedef Simplex<Dim> convex_type;
    /// Mesh of type convex_type
    typedef Mesh<convex_type> mesh_type;
    /// Pointer on the mesh
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

//typedef bases<Nedelec<0,NedelecKind::NED1> > basis_nedtype;
    /// Basis of \f$ [P_{Order}]^3 \f$
    typedef bases<Lagrange<Order, Vectorial> > basis_type;
    /// Function space of \f$ [P_{Order}]^3 \f$ on the mesh
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    /// Element of \f$ [P_{Order}]^3 \f$
    typedef typename space_type::element_type element_type;
    /// Pointer on the function space of \f$ [P_{Order}]^3 \f$
    typedef boost::shared_ptr<space_type> space_ptrtype;

    /// Array of element of \f$ [P_{Order}]^3 \f$
    std::vector<element_type> test;

    /// Constructs an object of type TestCurl
    TestCurl( mesh_ptrtype );
    /// Test the operator curl()
    void run();

private:
    /// The mesh used for the application
    mesh_ptrtype mesh;
    /// The vectorial function space used
    space_ptrtype Vh;
};

#endif
