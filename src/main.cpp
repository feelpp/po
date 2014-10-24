#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "psi0.hpp"
#include "eigenprob.hpp"
#include "spectralproblem.hpp"

using namespace Feel;
using namespace Eigen;


inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "needP0", po::value<bool>()->default_value( true ), "need to compute psi0" )
        ( "computeP0", po::value<bool>()->default_value( true ), "compute psi0 or load it" )
        ( "radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius" )

        ( "needEigen", po::value<bool>()->default_value( true ), "need the eigen modes" )
        ( "computeEigen", po::value<bool>()->default_value( true ), "compute the eigen modes or load them" )
        ( "useCurl", po::value<bool>()->default_value( true ), "use curl or grad form" )
        ( "usePresDiv", po::value<bool>()->default_value( false ), "use a pressure term in div form" )
        ( "usePresGrad", po::value<bool>()->default_value( true ), "use a pressure term in grad form" )
        ( "useDiric", po::value<bool>()->default_value( false ), "use Dirichlet condition" )
        ( "bccurln", po::value<bool>()->default_value( true ), "need boundary condition curl g.n (eigenlap)" )
        ( "bcn", po::value<bool>()->default_value( true ), "need boundary condition1 g.n (eigenlap)" )
        ( "divdiv", po::value<bool>()->default_value( true ), "need divdiv term" )
        ( "needDecomp", po::value<bool>()->default_value( true ), "need to decompose the eigen modes" )
        ( "computeDecomp", po::value<bool>()->default_value( true ), "compute the decomposition of the modes or load them" )
        ( "meanPsi", po::value<double>()->default_value( 1. ), "psi average" )
        ( "needDebug", po::value<bool>()->default_value( false ), "debug" )

        ( "needPS", po::value<bool>()->default_value( true ), "need to run the spectral problem" )
        ( "computeRijk", po::value<bool>()->default_value( false ), "compute or load Rijk" )
        ( "computeRiak", po::value<bool>()->default_value( false ), "compute or load Riak" )
        ( "computeRfk", po::value<bool>()->default_value( false ), "compute or load Rfk" )
        ( "f", po::value<std::string>()->default_value( "{0,0,1}" ), "f" )
        ( "computeRpk", po::value<bool>()->default_value( true ), "compute or load Rpk" )
        ( "nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" );
    return myappOptions;
}


po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ); // Poisson
    libOptions.add( backend_options( "gi0" ) ).add( backend_options( "psi" ) ).add( backend_options( "gradpsi" ) ).add( backend_options( "curl" ) ).add( backend_options( "curl2") ); // Eigen_Curl
    return libOptions.add( feel_options() );
}


AboutData
makeAbout()
{
    Feel::AboutData about( "po_app", "po_app" );
    about.addAuthor( "Romain Hild", "", "romain.hild@plasticomnium.com", "" );
    return about;
}


boost::shared_ptr<Mesh<Simplex<3> > >
load_mesh()
{
    Feel::fs::path mypath(soption( _name="gmsh.filename" ));
    std::string mesh_name = ( boost::format( "%1%.msh" )
                                 %mypath.stem().string() ).str();


    boost::shared_ptr<Mesh<Simplex<3> > > mesh;
    if(option(_name="computeEigen").as<bool>())
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _rebuild_partitions=(mypath.extension() == ".msh") );
    else
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _filename=mesh_name );

#if FEELPP_VERSION_GREATER_THAN(0,98,0)
    LOG(INFO) << " - mesh entities" << std::endl;
    LOG(INFO) << " number of elements : " << mesh->numGlobalElements() << std::endl;
    LOG(INFO) << " number of faces : " << mesh->numGlobalFaces() << std::endl;
    LOG(INFO) << " number of edges : " << mesh->numGlobalEdges() << std::endl;
    LOG(INFO) << " number of points : " << mesh->numGlobalPoints() << std::endl;
    LOG(INFO) << " number of vertices : " << mesh->numGlobalVertices() << std::endl;
    LOG(INFO) << " - mesh sizes" << std::endl;
    LOG(INFO) << " h max : " << mesh->hMax() << std::endl;
    LOG(INFO) << " h min : " << mesh->hMin() << std::endl;
    LOG(INFO) << " h avg : " << mesh->hAverage() << std::endl;
    LOG(INFO) << " measure : " << mesh->measure() << std::endl;

    if ( Environment::isMasterRank() )
    {
        std::cout << " - mesh entities" << std::endl;
        std::cout << " number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << " number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << " number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << " number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << " number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << " h max : " << mesh->hMax() << std::endl;
        std::cout << " h min : " << mesh->hMin() << std::endl;
        std::cout << " h avg : " << mesh->hAverage() << std::endl;
        std::cout << " measure : " << mesh->measure() << std::endl;
    }
#endif

    return mesh;
}

int
main( int argc, char **argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=makeLibOptions(),
                     _about=makeAbout() );


    boost::shared_ptr<Mesh<Simplex<3> > > mesh = load_mesh();


    Psi0 p0 = Psi0( mesh, soption( _name="alpha0" ) );
    EigenProb eig = EigenProb( mesh );
    SpectralProblem sp = SpectralProblem( mesh );

    auto e = exporter( _mesh=mesh );


    if( boption(_name="needP0") || boption(_name="needPS") ){
        p0.run();
        e->add( "a", p0.gradu );
    }

    if( boption(_name="needEigen") || boption(_name="needPS") ){
        eig.run();
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++){
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig.g[i] );
            if( boption( _name="needDecomp") ){
                e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig.psi[i] );
            }
            if( boption(_name="needDebug") ){
                e->add( ( boost::format( "gradpsi-%1%" ) % i ).str(), eig.gradu[i] );
                e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig.g0[i] );
                e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig.modebis[i] );
            }
        }
    }

    if( boption(_name="needPS") ){
        sp.init( eig.g, eig.psi, eig.lambda, p0.gradu );
        sp.run();
        e->add( "u", sp.u );
        e->add( "v", sp.v );
    }

    if( boption(_name="needPS") || boption(_name="needEigen") || boption(_name="needP0") )
        e->save();


    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End -----" << std::endl;


    return 0;
}
