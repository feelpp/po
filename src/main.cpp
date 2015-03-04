#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>

#include "solverns2.hpp"

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "needEigen", po::value<bool>()->default_value( true ), "need eigenmodes" )
        ( "markerList", po::value<std::vector<std::string> >()->multitoken(), "list of markers of the boundary" )
        ( "computeEigen", po::value<bool>()->default_value( true ), "need to compute eigenmodes, else load them" )
        ( "nbMode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ( "print", po::value<bool>()->default_value( false ), "print matrices" )

        ( "needA0", po::value<bool>()->default_value( true ), "need relief a0" )
        ( "computeA0", po::value<bool>()->default_value( true ), "need to compute a0, else load it" )
        ( "radius", po::value<double>()->default_value( 0.5 ), "cylinder's radius" )
        ( "speed", po::value<double>()->default_value( 1 ), "average speed" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends on x,y,radius,speed" )

        ( "needA1", po::value<bool>()->default_value( false ), "need relief a1" )
        ( "computeA1", po::value<bool>()->default_value( true ), "need to compute a1, else load it" )
        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "needA2", po::value<bool>()->default_value( false ), "need relief a2" )
        ( "computeA2", po::value<bool>()->default_value( true ), "need to compute a2, else load it" )
        ( "alpha2", po::value<std::string>()->default_value( "4.*speed/(radius*radius)" ), "alpha2, depends on speed and radius" )

        ( "needPS", po::value<bool>()->default_value( true ), "need to run the spectral problem" )
        ( "computeRijk", po::value<bool>()->default_value( false ), "compute or load Rijk" )
        ( "computeRiak", po::value<bool>()->default_value( false ), "compute or load Riak" )
        ( "computeRfk", po::value<bool>()->default_value( false ), "compute or load Rfk" )
        ( "f", po::value<std::string>()->default_value( "{0,0,1}" ), "f" )
        ( "computeRpk", po::value<bool>()->default_value( true ), "compute or load Rpk" )
        ( "nu", po::value<double>()->default_value( 1 ), "viscosity" )
        ;
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "a0" ) ).add( backend_options( "grada0" ) );
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

    SolverNS2 app;
    app.solve();

    return 0;
}
