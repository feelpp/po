#ifndef EIGENLAP_H
#define EIGENLAP_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>

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
    typedef bases<Lagrange<Order, Scalar>, Lagrange<Order, Scalar>, Lagrange<Order, Scalar>, Lagrange<Order, Scalar> > sBasis_type;
    //typedef bases<Lagrange<Order, Scalar>, Lagrange<Order, Scalar>, Lagrange<Order, Scalar> > sBasis_type;
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

    std::vector<vElement_type> g0;
    std::vector<vElement_type> gradu;
    std::vector<vElement_type> modebis;

    EigenLap( mesh_ptrtype );
    void run();
    void compute_eigens();
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
