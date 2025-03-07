#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>

#include "solverns2.hpp"

using namespace Feel;

AboutData
makeAbout()
{
    Feel::AboutData about( "po_online", "po_online" );
    about.addAuthor( "Romain Hild", "", "hild.romain@gmail.com", "" );
    return about;
}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "solverns2.verbose", po::value<int>()->default_value( 0 ), "level of verbosity" )
        ( "solverns2.path", po::value<std::string>()->default_value( "." ), "path to mesh, base and Rijk coeff" )
        ( "solverns2.nb-mode", po::value<int>()->default_value( 1 ), "number of modes to use" )
        ( "solverns2.format", po::value<std::string>()->default_value( "hdf5" ), "format in which load eigenfunctions (hdf5 binary, text)" )
        ( "solverns2.a-steady", po::value<bool>()->default_value( false ), "is a steady or not" )

        ( "solverns2.radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "solverns2.speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "solverns2.alpha0", po::value<std::string>()->default_value( "2.*speed*(1.-(x*x+y*y)/(radius*radius)):x:y:radius:speed:t" ), "alpha0, depends on x,y,radius,speed" )

        ( "solverns2.need-a1", po::value<bool>()->default_value( false ), "need relief a1" )
        ( "solverns2.alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "solverns2.need-a2", po::value<bool>()->default_value( false ), "need relief a2" )
        ( "solverns2.alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "solverns2.nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ( "solverns2.f", po::value<std::string>()->default_value( "{0,0,0}" ), "f" )

        ( "solverns2.compute-a0", po::value<bool>()->default_value( true ), "need to compute a0, else load it" )
        ( "solverns2.compute-a1", po::value<bool>()->default_value( true ), "need to compute a1, else load it" )
        ( "solverns2.compute-a2", po::value<bool>()->default_value( true ), "need to compute a2, else load it" )
        ( "solverns2.compute-raik", po::value<bool>()->default_value( false ), "compute or load Raik" )
        ( "solverns2.compute-riak", po::value<bool>()->default_value( true ), "compute or load Riak" )
        ( "solverns2.compute-rfk", po::value<bool>()->default_value( true ), "compute or load Rfk" )
        ( "solverns2.compute-rpk", po::value<bool>()->default_value( true ), "compute or load Rpk" )

        ( "solverns2.stokes", po::value<bool>()->default_value( true ), "compute Stokes if true, else compute Navier-Stokes" )
        ( "solverns2.newton-max-it", po::value<int>()->default_value( 20 ), "maximum iteration of Newton" )
        ( "solverns2.newton-tol", po::value<double>()->default_value( 1e-8 ), "tolerance for Newton" )

        ( "solverns2.start-time", po::value<double>()->default_value( 0.0 ), "start time" )
        ( "solverns2.time-step", po::value<double>()->default_value( 0.1 ), "time step" )
        ( "solverns2.final-time", po::value<double>()->default_value( 1.0 ), "final time" )

        ( "solverns2.v-exact", po::value<std::string>()->default_value( "{0,0,2*(1-4*(x*x + y*y))}:x:y:t"), "v exacte" )
        ;
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "a0" ) );//.add( backend_options( "grada0" ) );
    libOptions.add( backend_options( "a2" ) );
    libOptions.add( backend_options( "post" ) );
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
