#ifndef EIGENPROB_H
#define EIGENPROB_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>


using namespace Feel;
using namespace Feel::vf;


class EigenProb : public Application
{
    typedef Application super;

    static const uint16_type Dim = 3;
    static const uint16_type Order = 2;

 public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //typedef bases<Nedelec<0,NedelecKind::NED1>, Lagrange<1, Scalar> > basis_ptype;
    // [space]
    typedef bases<Lagrange<Order, Vectorial>, Lagrange<1, Scalar> > basis_ptype;
    // [space]
    typedef FunctionSpace<mesh_type, basis_ptype > space_ptype;


    typedef bases<Lagrange<Order, Scalar>, Lagrange<0, Scalar> > basis_mltype;
    typedef FunctionSpace<mesh_type, basis_mltype > space_mltype;


    //typedef bases<Nedelec<0,NedelecKind::NED1> > basis_vtype;
    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
    typedef typename space_vtype::element_type element_vtype;
    typedef boost::shared_ptr<space_vtype> space_vptrtype;


    typedef bases<Lagrange<Order, Scalar> > basis_stype;
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
    typedef typename space_stype::element_type element_stype;
    typedef boost::shared_ptr<space_stype> space_sptrtype;


    std::vector<element_vtype> g;
    std::vector<double> lambda;
    std::vector<element_stype> psi;

    // debug
    std::vector<element_vtype> g0;
    std::vector<element_vtype> gradu;
    std::vector<element_vtype> modebis;


    EigenProb( mesh_ptrtype );
    void run();

 private:
    void compute_eigens();
    void decomp();
    void load_eigens();


    mesh_ptrtype mesh;
    space_vptrtype Vh;
    space_sptrtype Sh;
    int nev;
};

#endif
