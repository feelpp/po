#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/exporter.hpp>

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
    typedef Lagrange<0, Vectorial> p0_fct_type;
    typedef Lagrange<1, Scalar> scalar1_fct_type;
    typedef Lagrange<2, Scalar> scalar2_fct_type;
    typedef Lagrange<2, Vectorial> vec_fct_type;
    typedef RaviartThomas<0> rt_fct_type;

    using p0_basis_type = bases<p0_fct_type>;
    using p0_space_type = FunctionSpace<mesh_type, p0_basis_type>;
    using p0_space_ptrtype = boost::shared_ptr<p0_space_type>;

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

    using solver_a_type = SolverA<vec_space_ptrtype, ml_space_ptrtype>;
    using solver_a_ptrtype = boost::shared_ptr<solver_a_type>;

    using solver_sp_type = SolverSpectralProblem<vec_space_ptrtype, ned_space_ptrtype, scalar_space_ptrtype, rt_space_ptrtype, eigentuple_type>;
    using solver_sp_ptrtype = boost::shared_ptr<solver_sp_type>;

    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = boost::shared_ptr<exporter_type>;

    void solve();

private:
    mesh_ptrtype mesh;

    eigen_space_ptrtype Xh;
    ned_space_ptrtype Nh;
    ml_space_ptrtype Mh;
    rt_space_ptrtype RTh;
    scalar_space_ptrtype Sh;
    vec_space_ptrtype Vh;
    p0_space_ptrtype P0;

    rt_element_type a;

    ned_element_type u;
    vec_element_type v;
    vec_element_type vex;
    vec_element_type verr;

    solver_a_ptrtype solverA;
    solver_sp_ptrtype solverSP;

    exporter_ptrtype e;

    int iter;
    double t;
    double dt;

    double err;

    void load_mesh();
    void initSpaces();
    void setA( double t );
    void solveSP( double t );
    void post( double t );
    void logInfo();
    void logMesh();
};

void
SolverNS2::solve()
{
    tic();

    load_mesh();

    initSpaces();

    solverA = SolverA<vec_space_ptrtype, ml_space_ptrtype>::build(mesh, Vh, Mh);
    solverSP = SolverSpectralProblem<vec_space_ptrtype, ned_space_ptrtype, scalar_space_ptrtype, rt_space_ptrtype, eigentuple_type>::build(mesh, Vh, Nh, Sh, RTh);
    solverSP->setEigen();
    solverSP->setRijk();

    auto Ih = I( _domainSpace=RTh, _imageSpace=P0);

    t = doption("solverns2.startTime");
    if( boption("solverns2.aSteady") )
    {
        setA( t );
        solverSP->setA( a );
        solverSP->init( t );
    }

    dt = doption( "solverns2.timeStep" );
    for( t = doption("solverns2.startTime"), iter = 0;
         t < doption("solverns2.finalTime");
         iter++, t += dt
         )
    {
        if( Environment::isMasterRank() && ioption("solverns2.verbose") > 0 )
            std::cout << "iteration " << iter << " at time " << t << "s" << std::endl;

        tic();
        if( !boption("solverns2.aSteady") )
            setA( t );
        e->step(t)->add("a", Ih(a));

        solveSP( t );
        e->step(t)->add("u", u);

        post( t );
        e->step(t)->add("v", v);
        e->step(t)->add("vex", vex);
        e->step(t)->add("verr", verr);

        e->save();
        toc("iteration", ioption("solverns2.verbose") > 0 );
    }

    toc("total", ioption("solverns2.verbose") > 0);
    // Environment::saveTimers(ioption("solverns2.verbose") > 1);
}

void
SolverNS2::load_mesh()
{
    tic();
    Feel::fs::path path(soption( _name="solverns2.path" ));
    auto meshname = path.append(soption( _name="gmsh.filename" ));
    mesh = loadGMSHMesh( _mesh=new mesh_type, _filename=meshname.string());

    e = exporter( mesh );

    logMesh();
    toc("mesh", ioption("solverns2.verbose") > 1);
}

void
SolverNS2::initSpaces()
{
    tic();
    Xh = eigen_space_type::New( mesh );
    Nh = Xh->template functionSpace<0>();
    Mh = ml_space_type::New( mesh );
    RTh = Mh->template functionSpace<0>();
    Sh = scalar_space_type::New( mesh );
    Vh = vec_space_type::New( mesh );
    P0 = p0_space_type::New( mesh );

    a = RTh->element();
    u = Nh->element();
    v = Vh->element();

    logInfo();
    toc("spaces", ioption("solverns2.verbose") > 1);
}

void
SolverNS2::setA( double t )
{
    tic();

    a = solverA->solve( t );

    toc("a", ioption("solverns2.verbose") > 1);
}

void
SolverNS2::solveSP( double t)
{
    tic();

    if( !boption("solverns2.aSteady") )
    {
        solverSP->setA(a);
        solverSP->init( t );
    }

    u = solverSP->solve( t );

    toc("sp", ioption("solverns2.verbose") > 1);
 }

void
SolverNS2::post( double t )
{
    tic();

    v = Vh->element();
    auto form2V = form2(_test=Vh, _trial=Vh);
    form2V = integrate(elements(mesh), inner(idt(v),id(v)));
    auto form1V = form1(_test=Vh);
    form1V = integrate( elements(mesh), inner(idv(a) + idv(u), id(v)));
    form2V.solve(_rhs=form1V, _solution=v);

    auto vex_expr = expr<3,1>(soption("solverns2.v_ex"));
    vex_expr.setParameterValues({{"t",t}});
    vex = Vh->element(vex_expr);

    auto uex = Vh->element();
    uex = vf::project( _range=elements(mesh), _space=Vh, _expr=idv(vex) - idv(a));
    verr = vf::project( _range=elements(mesh), _space=Vh, _expr=idv(u) + idv(a) - idv(vex));
    err = normL2(elements(mesh), idv(u)+idv(a) - idv(vex));
    auto errU = normL2(elements(mesh), idv(u) - idv(uex));

    LOG(INFO) << "error(" << solverSP->M << ") : " << err;
    if(Environment::isMasterRank())
    {
        std::cout << "errorV(" << solverSP->M << ") : " << err << std::endl;
        std::cout << "errorU(" << solverSP->M << ") : " << errU << std::endl;
        auto filename = (boost::format("%1%-info.md") %Environment::about().appName()).str();
        std::fstream s;
        s.open (filename, std::fstream::out|std::fstream::app);
        s << "\n#Results\n"
          << "error(" << solverSP->M << ") : " << err << std::endl;
        s.close();
    }


    // auto g_s = soption("solverns2.alpha0");
    // auto vars = Symbols{ "x", "y", "radius", "speed" };
    // auto g_e = parse( g_s, vars );
    // auto g = expr( g_e, vars );
    // auto g = expr(soption("solverns2.alpha0"));
    // g.setParameterValues(
    //     {
    //         { "radius", doption( "solverns2.radius" ) },
    //         { "speed", doption( "solverns2.speed" ) },
    //         { "t", t }
    //     } );

    // auto divvex = normL2(elements(mesh), divv(vex));
    // auto vexn = normL2(markedfaces(mesh,1), inner(idv(vex), N()) + g );
    // auto cvexn = normL2(boundaryfaces(mesh), inner(curlv(vex), N()));
    // auto divvh = normL2(elements(mesh), divv(a)+divv(u));
    // auto vn = normL2(markedfaces(mesh,1), inner(idv(a), N()) + inner(idv(u), N()) + g );
    // auto cvn = normL2(boundaryfaces(mesh), inner(curlv(a), N()) + inner(curlv(u), N()));
    // auto diva = normL2(elements(mesh), divv(a));
    // auto an = normL2(markedfaces(mesh,1), inner(idv(a), N()) + g );
    // auto can = normL2(boundaryfaces(mesh), inner(curlv(a), N()));
    // auto divu = normL2(elements(mesh), divv(u));
    // auto un = normL2(boundaryfaces(mesh), inner(idv(uex), N()));
    // auto cun = normL2(boundaryfaces(mesh), inner(curlv(uex), N()));

    // if( Environment::isMasterRank() )
    //     std::cout << "||div vex|| = " << divvex << std::endl
    //               << "||vex.n-a0||= " << vexn << std::endl
    //               << "||cvex.n||  = " << cvexn << std::endl
    //               << "||div v||   = " << divvh << std::endl
    //               << "||v.n-a0||  = " << vn << std::endl
    //               << "||cv.n||    = " << cvn << std::endl
    //               << "||div a||   = " << diva << std::endl
    //               << "||a.n-a0||  = " << an << std::endl
    //               << "||ca.n||    = " << can << std::endl
    //               << "||div u||   = " << divu << std::endl
    //               << "||u.n||     = " << un << std::endl
    //               << "||cu.n||    = " << cun << std::endl;

    toc("post", ioption("solverns2.verbose") > 1);
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
