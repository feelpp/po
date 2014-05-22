#ifndef EIGENLAPZ_H
#define EIGENLAPZ_H

#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>

using namespace Feel;
using namespace Feel::vf;

class EigenLapZ : public Application
{
    typedef Application super;

 public:
    typedef Simplex<3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<2, Scalar> > sBasis_type;
    typedef FunctionSpace<mesh_type, sBasis_type> sSpace_type;
    typedef sSpace_type::element_type sElement_type;
    typedef boost::shared_ptr<sSpace_type> sSpace_ptrtype;

    typedef bases<Lagrange<2, Vectorial> > vBasis_type;
    typedef FunctionSpace<mesh_type, vBasis_type > vSpace_type;
    typedef vSpace_type::element_type vElement_type;
    typedef boost::shared_ptr<vSpace_type> vSpace_ptrtype;

    // [typedef]
    typedef bases<Lagrange<2, Scalar>, Lagrange<1, Scalar> > mlBasis_type;
    typedef FunctionSpace<mesh_type, mlBasis_type> mlSpace_type;
    // [typedef]

    std::vector<vElement_type> g;
    std::vector<double> lambda;
    std::vector<sElement_type> psi;

    std::vector<vElement_type> g0;
    std::vector<vElement_type> gradu;
    std::vector<vElement_type> modebis;

    EigenLapZ( mesh_ptrtype );
    void run();
    void compute_eigens();
    void load_eigens();
    void testBessel();

 private:
    mesh_ptrtype mesh;
    vSpace_ptrtype Vh;
    sSpace_ptrtype Sh;
    int nev;
    int ncv;
};

#endif
