#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>

#include "solverns2.hpp"

using namespace Feel;

AboutData
makeAbout()
{
    Feel::AboutData about( "po_app", "po_app" );
    about.addAuthor( "Romain Hild", "", "romain.hild@plasticomnium.com", "" );
    return about;
}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "solverns2.verbose", po::value<int>()->default_value( 0 ), "level of verbosity" )

        ( "solverns2.markerList", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "solverns2.nbMode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "solverns2.print", po::value<bool>()->default_value( false ), "print matrices" )
        ( "solverns2.exportEigen", po::value<bool>()->default_value( false ), "export eigen modes" )

        ( "solverns2.aSteady", po::value<bool>()->default_value( false ), "is a steady or not" )

        ( "solverns2.radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "solverns2.speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "solverns2.alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius,speed" )

        ( "solverns2.needA1", po::value<bool>()->default_value( false ), "need relief a1" )
        ( "solverns2.alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "solverns2.needA2", po::value<bool>()->default_value( false ), "need relief a2" )
        ( "solverns2.alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "solverns2.nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "solverns2.f", po::value<std::string>()->default_value( "{0,0,1}" ), "f" )

        // if loadMesh = false, all compute options must be true !!!!
        ( "solverns2.loadMesh", po::value<bool>()->default_value( true ), "load the mesh or create it" )
        ( "solverns2.computeEigen", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "solverns2.computeA0", po::value<bool>()->default_value( true ), "need to compute a0, else load it" )
        ( "solverns2.computeA1", po::value<bool>()->default_value( true ), "need to compute a1, else load it" )
        ( "solverns2.computeA2", po::value<bool>()->default_value( true ), "need to compute a2, else load it" )
        ( "solverns2.computeRijk", po::value<bool>()->default_value( false ), "compute or load Rijk" )
        ( "solverns2.computeRaik", po::value<bool>()->default_value( false ), "compute or load Raik" )
        ( "solverns2.computeRiak", po::value<bool>()->default_value( true ), "compute or load Riak" )
        ( "solverns2.computeRfk", po::value<bool>()->default_value( true ), "compute or load Rfk" )
        ( "solverns2.computeRpk", po::value<bool>()->default_value( true ), "compute or load Rpk" )
        // if loadMesh = false, all compute options must be true !!!!

        ( "solverns2.stokes", po::value<bool>()->default_value( true ), "compute Stokes if true, else compute Navier-Stokes" )
        ( "solverns2.newton.maxIt", po::value<int>()->default_value( 20 ), "maximum iteration of Newton" )
        ( "solverns2.newton.tol", po::value<double>()->default_value( 1e-8 ), "tolerance for Newton" )

        ( "solverns2.startTime", po::value<double>()->default_value( 0.0 ), "start time" )
        ( "solverns2.timeStep", po::value<double>()->default_value( 0.1 ), "time step" )
        ( "solverns2.finalTime", po::value<double>()->default_value( 1.0 ), "final time" )

        ( "solverns2.v_ex", po::value<std::string>()->default_value( "{0,0,2*(1-4*(x*x + y*y))}:x:y"), "v exacte" )

        ( "solverns2.testEigs", po::value<bool>()->default_value( false ), "test eigenmodes" )
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

int
main( int argc, char **argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=makeLibOptions(),
                     _about=makeAbout() );

    SolverNS2 app;
    app.solve();

    return 0;
}
