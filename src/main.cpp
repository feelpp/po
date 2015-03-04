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
        ( "print", po::value<bool>()->default_value( false ), "print matrices" );
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
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
