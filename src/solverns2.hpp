#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "solvereigenns2.hpp"
#include "solvera0.hpp"
// #include "solvera1.hpp"
// #include "solvera2.hpp"
#include "solverspectralproblem.hpp"

using namespace Feel;
using namespace Feel::vf;

class SolverNS2
{
public:
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Nedelec<0, NedelecKind::NED1> ned_basis_type;
    typedef Lagrange<1,Scalar> lag1scalar_basis_type;
    typedef Lagrange<2,Vectorial> lag2vec_basis_type;

    typedef bases<ned_basis_type, lag1scalar_basis_type> eigen_basis_type;
    typedef FunctionSpace<mesh_type, eigen_basis_type> eigen_space_type;
    typedef boost::shared_ptr<eigen_space_type> eigen_space_ptrtype;

    typedef typename eigen_space_type::template sub_functionspace<0>::type ned_space_type;
    typedef boost::shared_ptr<ned_space_type> ned_space_ptrtype;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef std::pair<value_type,ned_element_type> eigenpair_type;
    typedef std::vector<eigenpair_type> eigenmodes_type;

    typedef FunctionSpace<mesh_type, bases<lag2vec_basis_type> > lag2vec_space_type;
    typedef boost::shared_ptr<lag2vec_space_type> lag2vec_space_ptrtype;
    typedef lag2vec_space_type::element_type lag2vec_element_type;

    void solve();

private:
    mesh_ptrtype mesh;

    eigen_space_ptrtype Xh;
    ned_space_ptrtype Nh;
    eigenmodes_type eigenModes;

    lag2vec_space_ptrtype Vh;
    lag2vec_element_type a0;

    lag2vec_element_type u;
    lag2vec_element_type v;

    void setEigen();
    void setA0();
    void setA1();
    void setA2();
    void solveSP();
};

void
SolverNS2::solve()
{
    mesh = loadMesh( _mesh = new mesh_type,
                     _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS
                     );

    if ( Environment::isMasterRank() && FLAGS_v > 0)
    {
        std::cout << " number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << " number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << " number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << " number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << " number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " number of local vertices : " << mesh->numLocalVertices() << std::endl;
    }

    auto e = exporter( mesh );

    if( boption("needEigen"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute eigenmodes ----------\n";

        setEigen();
        int i = 0;
        for( auto const& pair : eigenModes)
            e->add( ( boost::format( "mode-%1%" ) % i++ ).str(), pair.second );
    }

    if( boption("needA0"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute a0 ----------\n";

        setA0();
        e->add( "a0", a0);
    }

    setA1();
    setA2();

    if( boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute spectral problem ----------\n";

        solveSP();
        e->add("u", u);
    }

    e->save();
}

void
SolverNS2::setEigen()
{
    Xh = eigen_space_type::New( mesh );
    Nh = Xh->template functionSpace<0>();

    auto solverEigen = SolverEigenNS2<decltype(Xh)>::build(mesh, Xh);
    eigenModes = solverEigen->solve();
}

void
SolverNS2::setA0()
{
    Vh = lag2vec_space_type::New( mesh );

    auto solvera0 = SolverA0<decltype(Vh)>::build(mesh, Vh);
    a0 = solvera0->solve();
}

void
SolverNS2::setA1()
{
}

void
SolverNS2::setA2()
{
}

void
SolverNS2::solveSP()
{
    auto solversp = SolverSpectralProblem<decltype(Nh), decltype(Vh)>::build(mesh, Nh, Vh);
    solversp->setA0(a0);
    solversp->setEigen(eigenModes);

    if ( Environment::isMasterRank() && FLAGS_v > 0)
        std::cout << " ---------- init R coeff ----------\n";
    solversp->init();

    if ( Environment::isMasterRank() && FLAGS_v > 0)
        std::cout << " ---------- solve spectral problem ----------\n";
    u = solversp->solve();
}
