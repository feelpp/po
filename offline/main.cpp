#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>

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

        ( "eigen.marker-list", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "eigen.nb-mode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "eigen.print", po::value<bool>()->default_value( false ), "print matrices" )
        ( "eigen.export", po::value<bool>()->default_value( false ), "export eigen modes" )
        ( "eigen.test", po::value<bool>()->default_value( false ), "test eigenmodes" )

        // if loadMesh = false, all compute options must be true !!!!
        ( "offline.load-mesh", po::value<bool>()->default_value( true ), "load the mesh or create it" )
        ( "eigen.compute", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "coeff.compute", po::value<bool>()->default_value( true ), "compute or load Rijk" )
        // if loadMesh = false, all compute options must be true !!!!
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

    if( Environment::isMasterRank() )
        std::cout << "Path : " << Feel::fs::current_path() << std::endl;

    auto mesh = loadMesh( _mesh=new mesh_type, _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_PROPAGATE_MARKERS );
    auto Xh = eigen_space_type::New( mesh );
    auto Nh = Xh->template functionSpace<0>();
    auto Sh = scalar_space_type::New( mesh );

    auto solverEigen = SolverEigenNS2<eigen_space_ptrtype, scalar_space_ptrtype>::build(mesh, Xh, Sh);
    auto eigenModes = solverEigen->solve();

    auto initCoeff = InitCoeff<decltype(eigenModes)>::build( eigenModes );
    initCoeff->initRijk();

    return 0;
}
