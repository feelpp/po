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
    std::vector<size_type> indexesToKeep;
    sparse_matrix_ptrtype C;
    vec_element_type a;
    scalar_element_type psi0;

    ned_element_type u;
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
    logTime(t, "mesh", ioption("solverns2.verbose") > 1);

    auto e = exporter( mesh );

    initSpaces();
    logTime(t, "spaces", ioption("solverns2.verbose") > 1);
    logInfo(LOG(INFO));
    logInfo(std::cout);

    if( boption("solverns2.needEigen") || boption("solverns2.needSP"))
    {
        setEigen();
        for( int i = 0; i < eigenModes.size(); i = i + 10)
        {
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), std::get<1>(eigenModes[i]) );
            e->add( ( boost::format( "psi-%1%" ) % i ).str(), std::get<2>(eigenModes[i]) );
        }
        logTime(t, "eigenmodes", ioption("solverns2.verbose") > 0);
    }

    if( boption("solverns2.needA0") || boption("solverns2.needA1") || boption("solverns2.needA2") || boption("solverns2.needSP"))
    {
        setA();
        e->add( "a", a);
        //e->add( "psi0", psi0);
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

    e->save();

    logTime(total, "total", ioption("solverns2.verbose") > 0);
}

void
SolverNS2::load_mesh()
{
    Feel::fs::path mypath(soption( _name="gmsh.filename" ));
    std::string mesh_name = ( boost::format( "%1%.msh" )
                              %mypath.stem().string() ).str();


    if( boption("solverns2.computeA0")
        && boption("solverns2.computeA1")
        && boption("solverns2.computeA2") )
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
    auto solverEigen = SolverEigenNS2<eigen_space_ptrtype, scalar_space_ptrtype>::build(mesh, Xh, Sh);
    eigenModes = solverEigen->solve();
    C = solverEigen->C;
    indexesToKeep = solverEigen->indexesToKeep;
}

void
SolverNS2::setA()
{
    auto solvera = SolverA<vec_space_ptrtype, ml_space_ptrtype>::build(mesh, Vh, Mh);
    a = solvera->solve();
    //psi0 = solvera->psi0;
}

void
SolverNS2::solveSP()
{
    auto solversp = SolverSpectralProblem<vec_space_ptrtype, eigen_space_ptrtype>::build(mesh, Vh, Xh);
    solversp->setA(a);
    solversp->setEigen(C, indexesToKeep);

    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
        std::cout << " ---------- solve spectral problem ----------\n";
    u = solversp->solve();
}

void
SolverNS2::post()
{
    v = Vh->element();
    auto form2V = form2(_test=Vh, _trial=Vh);
    form2V = integrate(elements(mesh), inner(idt(v),id(v)));
    auto form1V = form1(_test=Vh);
    form1V = integrate( elements(mesh), inner(idv(a) + idv(u), id(v)));
    form2V.solve(_rhs=form1V, _solution=v);

    vex = Vh->element(expr<3,1>(soption("solverns2.v_ex")));
    auto err = normL2(elements(mesh), idv(v)-idv(vex));
    auto errC = normL2(elements(mesh), curlv(v) - curlv(vex));

    LOG(INFO) << "[" << eigenModes.size() << "] ||v-vex|| = " << err;
    LOG(INFO) << "[" << eigenModes.size() << "] ||curl(v)-curl(vex)|| = " << errC;
    if(Environment::isMasterRank())
        std::cout << "[" << eigenModes.size() << "] ||v-vex|| = " << err << std::endl
                  << "[" << eigenModes.size() << "] ||curl(v)-curl(vex)|| = " << errC << std::endl;
}$

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
            << "[info] Vh dof = " << Vh->nDof() << std::endl
            << "[info] Sh dof = " << Sh->nDof() << std::endl;
    }
}
