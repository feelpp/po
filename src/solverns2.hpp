#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/exporter.hpp>

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
    typedef Lagrange<2, Vectorial> vec_fct_type;
    typedef RaviartThomas<0> rt_fct_type;

    typedef bases<ned_fct_type, scalar1_fct_type> eigen_basis_type;
    typedef FunctionSpace<mesh_type, eigen_basis_type> eigen_space_type;
    typedef boost::shared_ptr<eigen_space_type> eigen_space_ptrtype;

    typedef typename eigen_space_type::template sub_functionspace<0>::type ned_space_type;
    typedef boost::shared_ptr<ned_space_type> ned_space_ptrtype;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef bases<rt_fct_type, scalar1_fct_type> ml_basis_type;
    typedef FunctionSpace<mesh_type, ml_basis_type> ml_space_type;
    typedef boost::shared_ptr<ml_space_type> ml_space_ptrtype;

    typedef typename ml_space_type::template sub_functionspace<0>::type rt_space_type;
    typedef boost::shared_ptr<rt_space_type> rt_space_ptrtype;
    typedef typename rt_space_type::element_type rt_element_type;

    typedef bases<scalar2_fct_type> scalar2_basis_type;
    typedef FunctionSpace<mesh_type, scalar2_basis_type> scalar_space_type;
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
    rt_space_ptrtype RTh;
    scalar_space_ptrtype Sh;
    vec_space_ptrtype Vh;

    eigenmodes_type eigenModes;
    rt_element_type a;

    vec_element_type u;
    vec_element_type v;
    vec_element_type vex;
    vec_element_type verr;

    double err;

    void load_mesh();
    void initSpaces();
    void setEigen();
    void setA();
    void solveSP();
    void post();
    void logInfo();
    void logMesh();
};

void
SolverNS2::solve()
{
    tic();
    tic();

    load_mesh();
    toc("mesh", ioption("solverns2.verbose") > 1);
    tic();

    auto e = exporter( mesh );

    initSpaces();
    a = RTh->element();
    u = Vh->element();
    v = Vh->element();
    toc("spaces", ioption("solverns2.verbose") > 1);
    logInfo();

    if( boption("solverns2.needEigen") || boption("solverns2.needSP"))
    {
        tic();
        setEigen();
        for( int i = 0; i < eigenModes.size(); i = i + 10)
        {
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), std::get<1>(eigenModes[i]) );
            e->add( ( boost::format( "psi-%1%" ) % i ).str(), std::get<2>(eigenModes[i]) );
        }
        toc("eigenmodes", ioption("solverns2.verbose") > 0);
    }

    if( boption("solverns2.needA0") || boption("solverns2.needA1") || boption("solverns2.needA2") || boption("solverns2.needSP"))
    {
        tic();
        setA();
        e->add( "a", a);
        toc("a", ioption("solverns2.verbose") > 0);
    }

    if( boption("solverns2.needSP"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute spectral problem ----------\n";

        tic();
        solveSP();
        e->add("u", u);
        toc("sp", ioption("solverns2.verbose") > 0);
        tic();
        post();
        e->add("v", v);
        e->add("vex", vex);
        e->add("verr", verr);
        toc("post", ioption("solverns2.verbose") > 0);
    }

    e->save();

    toc("total", ioption("solverns2.verbose") > 0);
    // Environment::saveTimers(ioption("solverns2.verbose") > 1);
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
        && boption("solverns2.computeRaik")
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

    logMesh();
}

void
SolverNS2::initSpaces()
{
    Xh = eigen_space_type::New( mesh );
    Nh = Xh->template functionSpace<0>();
    Mh = ml_space_type::New( mesh );
    RTh = Mh->template functionSpace<0>();
    Sh = scalar_space_type::New( mesh );
    Vh = vec_space_type::New( mesh );
}

void
SolverNS2::setEigen()
{
    auto solverEigen = SolverEigenNS2<eigen_space_ptrtype, scalar_space_ptrtype>::build(mesh, Xh, Sh);
    eigenModes = solverEigen->solve();
}

void
SolverNS2::setA()
{
    auto solvera = SolverA<vec_space_ptrtype, ml_space_ptrtype>::build(mesh, Vh, Mh);
    a = solvera->solve();
}

void
SolverNS2::solveSP()
{
    auto solversp = SolverSpectralProblem<vec_space_ptrtype, eigentuple_type>::build(mesh, Vh);
    //solversp->setA(a);
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
    v = Vh->element();
    auto form2V = form2(_test=Vh, _trial=Vh);
    form2V = integrate(elements(mesh), inner(idt(v),id(v)));
    auto form1V = form1(_test=Vh);
    form1V = integrate( elements(mesh), inner(idv(a) + idv(u), id(v)));
    form2V.solve(_rhs=form1V, _solution=v);

    vex = Vh->element(expr<3,1>(soption("solverns2.v_ex")));
    auto uex = Vh->element();
    uex = vf::project( _range=elements(mesh), _space=Vh, _expr=idv(vex) - idv(a));
    verr = vf::project( _range=elements(mesh), _space=Vh, _expr=idv(u) + idv(a) - idv(vex));
    err = normL2(elements(mesh), idv(u)+idv(a) - idv(vex));
    auto errU = normL2(elements(mesh), idv(u) - idv(uex));

    LOG(INFO) << "error(" << eigenModes.size() << ") : " << err;
    if(Environment::isMasterRank())
    {
        std::cout << "errorV(" << eigenModes.size() << ") : " << err << std::endl;
        std::cout << "errorU(" << eigenModes.size() << ") : " << errU << std::endl;
        auto filename = (boost::format("%1%-info.md") %Environment::about().appName()).str();
        std::fstream s;
        s.open (filename, std::fstream::out|std::fstream::app);
        s << "\n#Results\n"
          << "error(" << eigenModes.size() << ") : " << err << std::endl;
        s.close();
    }



    auto g_s = soption("solverns2.alpha0");
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", doption( "solverns2.radius" ) },
                { "speed", doption( "solverns2.speed" ) } } );

    auto divvex = normL2(elements(mesh), divv(vex));
    auto vexn = normL2(markedfaces(mesh,1), inner(idv(vex), N()) + g );
    auto cvexn = normL2(boundaryfaces(mesh), inner(curlv(vex), N()));
    auto divvh = normL2(elements(mesh), divv(a)+divv(u));
    auto vn = normL2(markedfaces(mesh,1), inner(idv(a), N()) + inner(idv(u), N()) + g );
    auto cvn = normL2(boundaryfaces(mesh), inner(curlv(a), N()) + inner(curlv(u), N()));
    auto diva = normL2(elements(mesh), divv(a));
    auto an = normL2(markedfaces(mesh,1), inner(idv(a), N()) + g );
    auto can = normL2(boundaryfaces(mesh), inner(curlv(a), N()));
    auto divu = normL2(elements(mesh), divv(u));
    auto un = normL2(boundaryfaces(mesh), inner(idv(uex), N()));
    auto cun = normL2(boundaryfaces(mesh), inner(curlv(uex), N()));

    if( Environment::isMasterRank() )
        std::cout << "||div vex|| = " << divvex << std::endl
                  << "||vex.n-a0||= " << vexn << std::endl
                  << "||cvex.n||  = " << cvexn << std::endl
                  << "||div v||   = " << divvh << std::endl
                  << "||v.n-a0||  = " << vn << std::endl
                  << "||cv.n||    = " << cvn << std::endl
                  << "||div a||   = " << diva << std::endl
                  << "||a.n-a0||  = " << an << std::endl
                  << "||ca.n||    = " << can << std::endl
                  << "||div u||   = " << divu << std::endl
                  << "||u.n||     = " << un << std::endl
                  << "||cu.n||    = " << cun << std::endl;
}

void
SolverNS2::logInfo()
{
    LOG(INFO) << "[info] np = " << Environment::numberOfProcessors() << std::endl
              << "[info] geo = " << soption("gmsh.filename") << std::endl
              << "[info] h = " << doption("gmsh.hsize") << std::endl
              << "[info] elt = " << mesh->numGlobalElements() << std::endl
              << "[info] Nh dof = " << Nh->nDof() << std::endl
              << "[info] Vh dof = " << Vh->nDof() << std::endl
              << "[info] Sh dof = " << Sh->nDof() << std::endl;

    if(Environment::isMasterRank())
    {
        auto filename = (boost::format("%1%-info.md") %Environment::about().appName()).str();
        std::fstream s;
        s.open (filename, std::fstream::out|std::fstream::app);
        s << "\n#Environment\n"
          << "x | Environment\n"
          << ":-: | :-:\n"
          << "geo | " << soption("gmsh.filename") << "\n"
          << "h | " << doption("gmsh.hsize") << "\n"
          << "np | " << Environment::numberOfProcessors() << "\n";
        if( boption("solverns2.computeEigen"))
            s << "modes | " << ioption("solvereigen.nev") << "\n";
        else
            s << "modes | " << ioption("solverns2.nbMode") << "\n";
        s << "alpha0 | " << soption("solverns2.alpha0") << "\n"
          << "alpha2 | " << soption("solverns2.alpha2") << "\n";
        s << "\n#Spaces\n"
          << "x | Nh | Vh | Sh\n"
          << ":-: | :-: | :-: | :-:\n"
          << "global dof | " << Nh->nDof() << " | " << Vh->nDof() << " | " << Sh->nDof() << "\n"
          << "local dof | " << Nh->nLocalDof() << " | " << Vh->nLocalDof() << " | " << Sh->nLocalDof() << std::endl;
        s.close();
    }

}

void SolverNS2::logMesh()
{
    LOG(INFO) << "[mesh]  - mesh entities" << std::endl
              << "[mesh]  number of elements : " << mesh->numGlobalElements() << std::endl
              << "[mesh]  number of faces : " << mesh->numGlobalFaces() << std::endl
              << "[mesh]  number of edges : " << mesh->numGlobalEdges() << std::endl
              << "[mesh]  number of points : " << mesh->numGlobalPoints() << std::endl
              << "[mesh]  number of vertices : " << mesh->numGlobalVertices() << std::endl
              << "[mesh]  - mesh sizes" << std::endl
              << "[mesh]  h max : " << mesh->hMax() << std::endl
              << "[mesh]  h min : " << mesh->hMin() << std::endl
              << "[mesh]  h avg : " << mesh->hAverage() << std::endl
              << "[mesh]  measure : " << mesh->measure() << std::endl;

    for( auto marker: mesh->markerNames() )
        LOG(INFO) << "[mesh]  - Marker " << marker.first << std::endl;

    if(Environment::isMasterRank())
    {
        auto filename = (boost::format("%1%-info.md") %Environment::about().appName()).str();
        std::fstream s;
        s.open (filename, std::fstream::out);
        s << "#Mesh\n"
          << "x | " << soption( _name="gmsh.filename" ) << "\n"
          << ":-: | :-:\n"
          << "number of elements | " << mesh->numGlobalElements() << "\n"
          << "number of faces | " << mesh->numGlobalFaces() << "\n"
          << "number of edges | " << mesh->numGlobalEdges() << "\n"
          << "number of points | " << mesh->numGlobalPoints() << "\n"
          << "number of vertices | " << mesh->numGlobalVertices() << "\n"
          << "h max | " << mesh->hMax() << "\n"
          << "h min | " << mesh->hMin() << "\n"
          << "h avg | " << mesh->hAverage() << "\n"
          << "measure | " << mesh->measure() << "\n";

        for( auto marker: mesh->markerNames() )
            s << "marker | " << marker.first << std::endl;
        s.close();
    }
}
