#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelvf/vf.hpp>

#include "solvereigenns2.hpp"
#include "initcoeff.hpp"

using namespace Feel;

AboutData
makeAbout()
{
    Feel::AboutData about( "po_offline", "po_offline" );
    about.addAuthor( "Romain Hild", "", "hild.romain@gmail.com", "" );
    return about;
}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "offline.verbose", po::value<int>()->default_value( 0 ), "level of verbosity" )

        ( "eigen.compute", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "eigen.marker-list", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "eigen.nb-mode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "eigen.format", po::value<std::string>()->default_value( "hdf5" ), "format in which save eigenfunctions (hdf5, binary, text)" )
        ( "eigen.print", po::value<bool>()->default_value( false ), "print matrices" )
        ( "eigen.export", po::value<bool>()->default_value( false ), "export eigen modes" )
        ( "eigen.test", po::value<bool>()->default_value( false ), "test eigenmodes" )
        ;
    return myappOptions;
}

int
main( int argc, char **argv )
{
    using value_type = double;
    using mesh_type = Mesh<Simplex<3> >;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    using ned_fct_type = Nedelec<0, NedelecKind::NED1>;
    using scalar1_fct_type = Lagrange<1, Scalar>;
    using scalar2_fct_type = Lagrange<2, Scalar>;

    using eigen_basis_type = bases<ned_fct_type, scalar1_fct_type>;
    using eigen_space_type = FunctionSpace<mesh_type, eigen_basis_type>;
    using eigen_space_ptrtype = boost::shared_ptr<eigen_space_type>;

    using scalar2_basis_type = bases<scalar2_fct_type>;
    using scalar_space_type = FunctionSpace<mesh_type, scalar2_basis_type>;
    using scalar_space_ptrtype = boost::shared_ptr<scalar_space_type>;
    using scalar_element_type = typename scalar_space_type::element_type;


    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=feel_options(),
                     _about=makeAbout() );

    tic();
    auto filename = (boost::format("%1%-info.md") %Environment::about().appName()).str();
    std::fstream s;
    if(Environment::isMasterRank())
    {
        s.open (filename, std::fstream::out);
        s << "#Offline\n";
        for( int i = 0; i < argc; i++ )
            s << argv[i] << " ";
        s << std::endl;
    }

    tic();
    auto mesh = loadMesh( _mesh=new mesh_type, _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS );
    if(Environment::isMasterRank())
    {
        s << "\n##Mesh\n"
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
    }
    toc("mesh", ioption("offline.verbose") > 0 );

    tic();
    auto Xh = eigen_space_type::New( mesh );
    auto Nh = Xh->template functionSpace<0>();
    auto Sh = scalar_space_type::New( mesh );
    if(Environment::isMasterRank())
    {
        s << "\n##Space\n"
          << "space | dof | np | local dof\n"
          << ":-: | :-: | :-: | :-:\n"
          << "Nh | " << Nh->nDof() << " | " << Environment::numberOfProcessors() << " | " << Nh->nLocalDof() << std::endl;
    }
    toc("spaces", ioption("offline.verbose") > 0 );

    auto solverEigen = SolverEigenNS2<eigen_space_ptrtype, scalar_space_ptrtype>::build(mesh, Xh, Sh);
    auto eigenModes = solverEigen->solve();

    auto initCoeff = InitCoeff<decltype(eigenModes)>::build( eigenModes );
    initCoeff->initRijk();
    toc("total", ioption("offline.verbose") > 0 );

    auto gn = normL2(_range=boundaryfaces(mesh), _expr=trans(idv(std::get<1>(eigenModes[0])))*N());
    auto divg = normL2(_range=elements(mesh), _expr=divv(std::get<1>(eigenModes[0])));
    auto curlgn = normL2(_range=boundaryfaces(mesh), _expr=trans(curlv(std::get<1>(eigenModes[0])))*N());
    // auto g0xn = normL2(_range=boundaryfaces(mesh), _expr=cross(idv(eigenModes[0]),N()));
    // auto g0n = normL2(_range=boundaryfaces(mesh), _expr=trans(idv(eigenModes[0]))*N());
    // auto g0 = normL2(_range=boundaryfaces(mesh), _expr=idv(eigenModes[0]));
    // auto e = normL2(_range=elements(mesh), _expr=idv(std::get<1>(eigenModes[0]))-(idv(eigenModes[0])+trans(gradv(std::get<2>(eigenModes[0])))));
    auto e2 = normL2(_range=elements(mesh), _expr=curlv(std::get<1>(eigenModes[0]))-std::sqrt(std::get<0>(eigenModes[0]))*idv(std::get<1>(eigenModes[0])));
    auto curl2 = integrate(_range=elements(mesh), _expr=trans(curlv(std::get<1>(eigenModes[0])))*curlv(std::get<1>(eigenModes[0]))).evaluate()(0,0);
    auto norm = normL2(_range=elements(mesh), _expr=idv(std::get<1>(eigenModes[0])));
    auto psi = integrate(_range=boundaryfaces(mesh), _expr=idv(std::get<2>(eigenModes[0]))).evaluate()(0,0);
    auto psin = normL2(_range=boundaryfaces(mesh), _expr=gradv(std::get<2>(eigenModes[0]))*N());


    if(Environment::isMasterRank())
    {
        s << "\n##Eigenmodes\n"
          << "norm(g.n) | " << gn << "\n"
          << "norm(divg) | " << divg << "\n"
          << "norm(curlg.n) | " << curlgn << "\n"
          // << "norm(g0xn) | " << g0xn << "\n"
          // << "norm(g0.n) | " << g0n << "\n"
          // << "norm(g0) | " << g0 << "\n"
          // << "norm(err) | " << e << "\n"
          << "norm(curlg-lg) | " << e2 << "\n"
          << "(curl,curl) | " << curl2 << "\n"
          << "norm(gi) | " << norm << "\n"
          << "norm(Gpsi.n) | " << psin << "\n"
          << "int psi | " << psi << std::endl;

        s.close();
    }


    if( Environment::isMasterRank() )
        std::cout << "Path : " << Feel::fs::current_path() << std::endl;

    return 0;
}
