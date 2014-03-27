#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "poisson.h"

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "rayon", po::value<double>()->default_value( 0.05 ), "rayon" )
        ( "vitesse", po::value<double>()->default_value( 0.015 ), "vitesse moyenne d'entree" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * vitesse * (1. - (x*x + y*y) / (rayon * rayon))" ), "alpha0, depend de x,y,rayon,vitesse" )
        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )
        ( "needEigen", po::value<bool>()->default_value( false ), "le besoin de recalculer les fonctions propres" )
        ( "pathEigen", po::value<std::string>()->default_value( "mode-" ), "path des modes" );
    return myappOptions;
}

AboutData
makeAbout()
{
    Feel::AboutData about( "po_app","po_app" );
    about.addAuthor( "Romain Hild", "", "romain.hild@plasticomnium.com", "" );
    return about;
}

int
main( int argc, char **argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=feel_options().add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ),
                     _about=makeAbout() );

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    Poisson p1 = Poisson(mesh, "0.", option(_name="alpha0").as<std::string>() );
    p1.run();

    auto e = exporter( _mesh=mesh );
    e->add( "u", p1.U.template element<0>() );
    e->add( "grad_u", p1.gradu );
    e->save();

}
