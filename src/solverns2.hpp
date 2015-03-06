#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "utilities.h"
#include "solvereigenns2.hpp"
#include "solvera0.hpp"
// #include "solvera1.hpp"
// #include "solvera2.hpp"
#include "solverspectralproblem.hpp"

using namespace Feel;
using namespace Feel::vf;
                                                               // move logInfo here + retrieve Zhsize
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

    void load_mesh();
    void setEigen();
    void setA0();
    void setA1();
    void setA2();
    void solveSP();
    void logInfo(std::ostream& out);
};

void
SolverNS2::solve()
{
    load_mesh();

    auto e = exporter( mesh );

    if( boption("needEigen") || boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute eigenmodes ----------\n";

        setEigen();
        int i = 0;
        for( auto const& pair : eigenModes)
            e->add( ( boost::format( "mode-%1%" ) % i++ ).str(), pair.second );
    }

    if( boption("needA0") || boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute a0 ----------\n";

        setA0();
        e->add( "a0", a0);
    }

    if( boption("needA1") || boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute a1 ----------\n";

        setA1();
    }

    if( boption("needA2") || boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute a2 ----------\n";

        setA2();
    }

    if( boption("needSP"))
    {
        if ( Environment::isMasterRank() && FLAGS_v > 0)
            std::cout << " ---------- compute spectral problem ----------\n";

        solveSP();
        //e->add("u", u);
    }

    logInfo(LOG(INFO));
    logInfo(std::cout);
    e->save();
}

void
SolverNS2::load_mesh()
{
    Feel::fs::path mypath(soption( _name="gmsh.filename" ));
    std::string mesh_name = ( boost::format( "%1%.msh" )
                              %mypath.stem().string() ).str();


    if( boption("computeEigen")
        && boption("computeA0")
        && boption("computeA1")
        && boption("computeA2")
        && boption("computeRijk")
        && boption("computeRiak")
        && boption("computeRfk") )
    {
        mesh = loadMesh( _mesh=new mesh_type,
                         _rebuild_partitions=(mypath.extension() == ".msh"),
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS
                         );
    } else {
        mesh = loadMesh( _mesh=new mesh_type,
                         _filename=mesh_name,
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS
                         );
    }

    logMesh(LOG(INFO), mesh);
    logMesh(std::cout, mesh);
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

    // if ( Environment::isMasterRank() && FLAGS_v > 0)
    //     std::cout << " ---------- solve spectral problem ----------\n";
    // u = solversp->solve();
}

void
SolverNS2::logInfo(std::ostream& out)
{
    if(Environment::isMasterRank())
    {
        out << "[info] np = " << Environment::numberOfProcessors() << std::endl
            << "[info] geo = " << soption("gmsh.filename") << std::endl
            << "[info] h = " << doption("gmsh.hsize") << std::endl
            << "[info] elt = " << mesh->numGlobalElements() << std::endl;
        if(boption("needEigen") || boption("needSP"))
            out << "[info] Nh dof = " << Nh->nDof() << std::endl;
        if(boption("needA0") || boption("needSP"))
            out << "[info] Vh dof = " << Vh->nDof() << std::endl;
    }
}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "needEigen", po::value<bool>()->default_value( true ), "need eigenmodes" )
        ( "markerList", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "computeEigen", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "nbMode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "print", po::value<bool>()->default_value( false ), "print matrices" )

        ( "needA0", po::value<bool>()->default_value( false ), "need relief a0" )
        ( "computeA0", po::value<bool>()->default_value( true ), "need to compute a0, else load it" )
        ( "radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius,speed" )

        ( "needA1", po::value<bool>()->default_value( false ), "need relief a1" )
        ( "computeA1", po::value<bool>()->default_value( true ), "need to compute a1, else load it" )
        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "needA2", po::value<bool>()->default_value( false ), "need relief a2" )
        ( "computeA2", po::value<bool>()->default_value( true ), "need to compute a2, else load it" )
        ( "alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "needSP", po::value<bool>()->default_value( true ), "need to run the spectral problem" )
        ( "nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "computeRijk", po::value<bool>()->default_value( true ), "compute or load Rijk" )
        //( "computeRaik", po::value<bool>()->default_value( false ), "compute or load Raik" )
        ( "computeRiak", po::value<bool>()->default_value( false ), "compute or load Riak" )
        ( "computeRfk", po::value<bool>()->default_value( false ), "compute or load Rfk" )
        ( "f", po::value<std::string>()->default_value( "{0,0,1}" ), "f" )
        ;
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "a0" ) ).add( backend_options( "grada0" ) );
    return libOptions.add( feel_options() );
}

