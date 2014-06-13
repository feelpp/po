#ifndef EIGENLAP_H
#define EIGENLAP_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pchv.hpp>

using namespace Feel;
using namespace Feel::vf;

class EigenLap : public Application
{
    typedef Application super;
    static const uint16_type Order = 2;
    static const uint16_type Dim = 3;

 public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<Order, Vectorial>, Lagrange<1, Scalar> > basis_ptype;
    typedef FunctionSpace<mesh_type, basis_ptype > space_ptype;

    typedef bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > basis_mltype;
    typedef FunctionSpace<mesh_type, basis_mltype > space_mltype;

    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
    typedef space_vtype::element_type element_vtype;
    typedef boost::shared_ptr<space_vtype> space_vptrtype;

    typedef bases<Lagrange<Order, Scalar> > basis_stype;
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
    typedef space_stype::element_type element_stype;
    typedef boost::shared_ptr<space_stype> space_sptrtype;

    std::vector<element_vtype> g;
    std::vector<double> lambda;
    std::vector<element_stype> psi;

    std::vector<element_vtype> g0;
    std::vector<element_vtype> gradu;
    std::vector<element_vtype> modebis;

    EigenLap( mesh_ptrtype );
    void run();
    void compute_eigens();
    void load_eigens();

 private:
    mesh_ptrtype mesh;
    space_vptrtype Vh;
    space_sptrtype Sh;
    int nev;
};

#endif
