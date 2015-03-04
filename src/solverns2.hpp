#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "solvereigenns2.hpp"
// #include "solvera0.hpp"
// #include "solvera1.hpp"
// #include "solvera2.hpp"
// #include "solverspectralproblem.hpp"

using namespace Feel;
using namespace Feel::vf;

class SolverNS2
{
public:
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Nedelec<0, NedelecKind::NED1> basis_edge_type;
    typedef Lagrange<1,Scalar> basis_vertex_type;
    typedef bases<basis_edge_type, basis_vertex_type> basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::template sub_functionspace<0>::type ned_space_type;
    typedef boost::shared_ptr<ned_space_type> ned_space_ptrtype;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef std::pair<value_type,ned_element_type> eigenpair_type;
    typedef std::vector<eigenpair_type> eigenmodes_type;

    void solve();

private:
    mesh_ptrtype mesh;
    space_ptrtype Xh;
    ned_space_ptrtype Nh;
    eigenmodes_type eigenModes;

    void setEigen();
    void setA0();
    void setA1();
    void setA2();
    void initPS();
    void solvePS();
};

void
SolverNS2::solve()
{
    mesh = loadMesh( _mesh = new mesh_type,
                     _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS
                     );

    if ( Environment::isMasterRank() && FLAGS_v > 0){
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
        setEigen();
        int i = 0;
        for( auto const& pair : eigenModes)
            e->add( ( boost::format( "mode-%1%" ) % i++ ).str(), pair.second );
    }

    setA0();
    setA1();
    setA2();
    initPS();
    solvePS();

    e->save();
}

void
SolverNS2::setEigen()
{
    Xh = space_type::New( mesh );
    Nh = Xh->template functionSpace<0>();

    auto eigen = SolverEigenNS2<decltype(Xh)>::build(mesh, Xh);
    eigenModes = eigen->solve();
}

void
SolverNS2::setA0()
{

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
SolverNS2::initPS()
{

}

void
SolverNS2::solvePS()
{

}
