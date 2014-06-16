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
        ( "radius", po::value<double>()->default_value( 0.05 ), "cylinder's radius" )
        ( "speed", po::value<double>()->default_value( 0.015 ), "average speed" )
        ( "nu", po::value<double>()->default_value( 1.56e-5 ), "viscosity" )
        ( "f", po::value<std::string>()->default_value( "{0,0,1}:x:y:z" ), "f" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius,speed" )
        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )
        ( "alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )
        ( "penal", po::value<double>()->default_value( 1e-6 ), "penalization" )
        ( "needP0", po::value<bool>()->default_value( true ), "need to compute psi0" )
        ( "needEigen", po::value<bool>()->default_value( true ), "need to compute the eigen modes or to load them" )
        ( "needPS", po::value<bool>()->default_value( true ), "need to run the spectral problem" )
        ( "needDecomp", po::value<bool>()->default_value( true ), "need to decompose the eigen modes" )
        ( "needDebug", po::value<bool>()->default_value( false ), "debug" )
        ( "bccurln", po::value<bool>()->default_value( true ), "need boundary condition curl g.n (eigenlap)" )
        ( "bcn", po::value<bool>()->default_value( true ), "need boundary condition1 g.n (eigenlap)" )
        ( "divdiv", po::value<bool>()->default_value( true ), "need divdiv term" );
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ); // Poisson
    libOptions.add( backend_options( "gi0" ) ).add( backend_options( "psi" ) ).add( backend_options( "gradpsi" ) ).add( backend_options( "curl" ) ); // Eigen_Curl
    return libOptions.add( feel_options() );
}

AboutData
makeAbout()
{
    Feel::AboutData about( "po_app", "po_app" );
    about.addAuthor( "Romain Hild", "", "romain.hild@plasticomnium.com", "" );
    return about;
}

int
main( int argc, char **argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=makeLibOptions(),
                     _about=makeAbout() );

    fs::path mypath(soption( _name="gmsh.filename" ));
    std::string meshPartName = ( boost::format( "%1%/%2%_part%3%.msh" )
                                 %Environment::localGeoRepository()
                                 %mypath.stem().string()
                                 %Environment::numberOfProcessors() ).str();

    boost::shared_ptr<Mesh<Simplex<3> > > mesh;
    bool isMsh = mypath.extension() == ".msh";
    if(option(_name="needEigen").as<bool>())
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _rebuild_partitions=isMsh,
                         _rebuild_partitions_filename=meshPartName );
    else
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _filename=meshPartName );

    Psi0 p0 = Psi0( mesh, soption( _name="alpha0" ) );
    // Eigen_Curl eig = Eigen_Curl( mesh );
    EigenProb eig = EigenProb( mesh );
    // Darcy d = Darcy( mesh, soption( _name="alpha0" ) );
    SpectralProblem sp = SpectralProblem( mesh );

    auto e = exporter( _mesh=mesh );

    if( boption(_name="needP0") || boption(_name="needPS") ){
        p0.run();
        e->add( "psi0", p0.gradu );
    }
    if( boption(_name="needEigen") || boption(_name="needPS") ){
        eig.run();
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++){
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig.g[i] );
            if( boption( _name="needDecomp") ){
                e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig.g0[i] );
                e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig.psi[i] );
                if( boption(_name="needDebug") ){
                    e->add( ( boost::format( "gradpsi-%1%" ) % i ).str(), eig.gradu[i] );
                    e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig.modebis[i] );
                }
            }
        }
    }
    if( boption(_name="needPS") ){
        sp.init( eig.g, eig.psi, eig.lambda, p0.gradu );
        sp.run();
        e->add( "v", sp.u );
    }
    e->save();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End -----" << std::endl;
}
