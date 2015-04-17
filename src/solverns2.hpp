#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "utilities.h"
#include "solvereigenns2.hpp"
#include "solvera.hpp"
#include "solverspectralproblem.hpp"

using namespace Feel;
using namespace Feel::vf;

class SolverNS2
{
public:
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Nedelec<0, NedelecKind::NED1> ned_fct_type;
    typedef Lagrange<0, Scalar> cst_fct_type;
    typedef Lagrange<1, Scalar> scalar1_fct_type;
    typedef Lagrange<2, Scalar> scalar2_fct_type;
    typedef Lagrange<1, Vectorial> vec_fct_type;

    typedef bases<ned_fct_type, scalar1_fct_type> eigen_basis_type;
    typedef FunctionSpace<mesh_type, eigen_basis_type> eigen_space_type;
    typedef boost::shared_ptr<eigen_space_type> eigen_space_ptrtype;

    typedef bases<scalar2_fct_type, cst_fct_type> ml_basis_type;
    typedef FunctionSpace<mesh_type, ml_basis_type> ml_space_type;
    typedef boost::shared_ptr<ml_space_type> ml_space_ptrtype;

    typedef typename eigen_space_type::template sub_functionspace<0>::type ned_space_type;
    typedef boost::shared_ptr<ned_space_type> ned_space_ptrtype;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef typename ml_space_type::template sub_functionspace<0>::type scalar_space_type;
    typedef boost::shared_ptr<scalar_space_type> scalar_space_ptrtype;
    typedef typename scalar_space_type::element_type scalar_element_type;

    typedef std::tuple<value_type, ned_element_type, scalar_element_type> eigentuple_type;
    typedef std::vector<eigentuple_type> eigenmodes_type;

    typedef FunctionSpace<mesh_type, bases<vec_fct_type> > vec_space_type;
    typedef boost::shared_ptr<vec_space_type> vec_space_ptrtype;
    typedef vec_space_type::element_type vec_element_type;

    void solve();

private:
    mesh_ptrtype mesh;

    eigen_space_ptrtype Xh;
    ned_space_ptrtype Nh;
    ml_space_ptrtype Mh;
    scalar_space_ptrtype Sh;
    vec_space_ptrtype Vh;

    eigenmodes_type eigenModes;
    vec_element_type a;

    vec_element_type u;
    vec_element_type v;
    vec_element_type vex;

    void load_mesh();
    void initSpaces();
    void setEigen();
    void setA();
    void solveSP();
    void post();
    void logInfo(std::ostream& out);
};

void
SolverNS2::solve()
{
    boost::mpi::timer t;
    boost::mpi::timer total;

    load_mesh();
    logTime(t, "mesh", ioption("solverns2.verbose") > 0);

    auto e = exporter( mesh );

    initSpaces();
    logTime(t, "spaces", ioption("solverns2.verbose") > 0);

    if( boption("solverns2.needEigen") || boption("solverns2.needSP"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute eigenmodes ----------\n";

        setEigen();
        int i = 0;
        for( auto const& tuple : eigenModes)
        {
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), std::get<1>(tuple) );
            e->add( ( boost::format( "psi-%1%" ) % i++ ).str(), std::get<2>(tuple) );
        }
        logTime(t, "eigenmodes", ioption("solverns2.verbose") > 0);
    }

    if( boption("solverns2.needA0") || boption("solverns2.needA1") || boption("solverns2.needA2") || boption("solverns2.needSP"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute a ----------\n";

        setA();
        e->add( "a", a);
        logTime(t, "a", ioption("solverns2.verbose") > 0);
    }

    if( boption("solverns2.needSP"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute spectral problem ----------\n";

        solveSP();
        e->add("u", u);
        logTime(t, "sp", ioption("solverns2.verbose") > 0);
        post();
        e->add("v", v);
        e->add("vex", vex);
        logTime(t, "post", ioption("solverns2.verbose") > 0);
    }


    logInfo(LOG(INFO));
    logInfo(std::cout);
    e->save();

    logTime(total, "total", ioption("solverns2.verbose") > 0);
}

void
SolverNS2::load_mesh()
{
    Feel::fs::path mypath(soption( _name="gmsh.filename" ));
    std::string mesh_name = ( boost::format( "%1%.msh" )
                              %mypath.stem().string() ).str();


    if( boption("solverns2.computeEigen")
        && boption("solverns2.computeA0")
        && boption("solverns2.computeA1")
        && boption("solverns2.computeA2")
        && boption("solverns2.computeRijk")
        && boption("solverns2.computeRiak")
        && boption("solverns2.computeRfk") )
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
SolverNS2::initSpaces()
{
    Xh = eigen_space_type::New( mesh );
    Nh = Xh->template functionSpace<0>();
    Mh = ml_space_type::New( mesh );
    Sh = Mh->template functionSpace<0>();
    Vh = vec_space_type::New( mesh );
}

void
SolverNS2::setEigen()
{
    auto solverEigen = SolverEigenNS2<eigen_space_ptrtype, ml_space_ptrtype, vec_space_ptrtype>::build(mesh, Xh, Mh, Vh);
    eigenModes = solverEigen->solve();
}

void
SolverNS2::setA()
{
    auto solvera = SolverA<vec_space_ptrtype>::build(mesh, Vh);
    a = solvera->solve();
}

void
SolverNS2::solveSP()
{
    auto solversp = SolverSpectralProblem<vec_space_ptrtype, eigentuple_type>::build(mesh, Vh);
    solversp->setA(a);
    solversp->setEigen(eigenModes);

    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
        std::cout << " ---------- init R coeff ----------\n";
    solversp->init();

    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
        std::cout << " ---------- solve spectral problem ----------\n";
    u = solversp->solve();
}

void
SolverNS2::post()
{
    v = vf::project( _range=elements(mesh), _space=Vh, _expr=idv(a) + idv(u));
    vex = Vh->element(expr<3,1>(soption("solverns2.v_ex")));
    auto err = normL2(elements(mesh), idv(v)-idv(vex));
    if(Environment::isMasterRank())
        std::cout << "error(" << eigenModes.size() << ") : " << err << std::endl;
}

void
SolverNS2::logInfo(std::ostream& out)
{
    if(Environment::isMasterRank())
    {
        out << "[info] np = " << Environment::numberOfProcessors() << std::endl
            << "[info] geo = " << soption("gmsh.filename") << std::endl
            << "[info] h = " << doption("gmsh.hsize") << std::endl
            << "[info] elt = " << mesh->numGlobalElements() << std::endl
            << "[info] Nh dof = " << Nh->nDof() << std::endl
            << "[info] Vh dof = " << Vh->nDof() << std::endl;
    }
}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "solverns2.verbose", po::value<int>()->default_value( 0 ), "level of verbosity" )

        ( "solverns2.needEigen", po::value<bool>()->default_value( true ), "need eigenmodes" )
        ( "solverns2.markerList", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "solverns2.computeEigen", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "solverns2.nbMode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "solverns2.print", po::value<bool>()->default_value( false ), "print matrices" )

        ( "solverns2.needA0", po::value<bool>()->default_value( false ), "need relief a0" )
        ( "solverns2.computeA0", po::value<bool>()->default_value( true ), "need to compute a0, else load it" )
        ( "solverns2.radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "solverns2.speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "solverns2.alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius,speed" )

        ( "solverns2.needA1", po::value<bool>()->default_value( false ), "need relief a1" )
        ( "solverns2.computeA1", po::value<bool>()->default_value( true ), "need to compute a1, else load it" )
        ( "solverns2.alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "solverns2.needA2", po::value<bool>()->default_value( false ), "need relief a2" )
        ( "solverns2.computeA2", po::value<bool>()->default_value( true ), "need to compute a2, else load it" )
        ( "solverns2.alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "solverns2.needSP", po::value<bool>()->default_value( true ), "need to run the spectral problem" )
        ( "solverns2.nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "solverns2.computeRijk", po::value<bool>()->default_value( false ), "compute or load Rijk" )
        //( "solverns2.computeRaik", po::value<bool>()->default_value( false ), "compute or load Raik" )
        ( "solverns2.computeRiak", po::value<bool>()->default_value( true ), "compute or load Riak" )
        ( "solverns2.computeRfk", po::value<bool>()->default_value( true ), "compute or load Rfk" )
        ( "solverns2.computeRpk", po::value<bool>()->default_value( true ), "compute or load Rpk" )
        ( "solverns2.f", po::value<std::string>()->default_value( "{0,0,1}" ), "f" )

        ( "solverns2.v_ex", po::value<std::string>()->default_value( "{0,0,8*(1-(x*x + y*y))}:x:y"), "v exacte" )
        ;
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "a0" ) ).add( backend_options( "grada0" ) );
    libOptions.add( backend_options( "a2" ) );
    return libOptions.add( feel_options() );
}

