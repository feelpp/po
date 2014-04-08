#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "poisson.h"
#include "eigen_curl.h"
#include "darcy.h"
#include "spectralproblem.h"

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "radius", po::value<double>()->default_value( 0.05 ), "cylinder's radius" )
        ( "speed", po::value<double>()->default_value( 0.015 ), "average speed" )
        ( "alpha0", po::value<std::string>()->default_value( "2. * speed * (1. - (x*x + y*y) / (radius * radius))" ), "alpha0, depends of x,y,radius,speed" )
        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )
        ( "needEigen", po::value<bool>()->default_value( true ), "need to compute the eigen modes or to load them" )
        ( "nbApp", po::value<int>()->default_value( 0 ), "app to launch (0:Poisson, 1:Eigen problem, 2:Darcy)" );
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ); // Poisson
    libOptions.add( backend_options( "gi0" ) ).add( backend_options( "psi" ) ).add( backend_options( "gradpsi" ) ); // Eigen_Curl
    libOptions.add( backend_options( "psi0Div" ) ); // Darcy
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

    fs::path mypath(option( _name="gmsh.filename" ).as<std::string>());
    std::string meshPartName = ( boost::format( "%1%/%2%_part%3%.msh" )
                                 %Environment::localGeoRepository()
                                 %mypath.stem().string()
                                 %Environment::numberOfProcessors() ).str();

    boost::shared_ptr<Mesh<Simplex<3> > > mesh;
    if(option(_name="needEigen").as<bool>())
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _rebuild_partitions=true,
                         _rebuild_partitions_filename=meshPartName );
    else
        mesh = loadMesh( _mesh=new Mesh<Simplex<3> >,
                         _filename=meshPartName );

    //auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );
    auto e = exporter( _mesh=mesh );

    switch(option(_name="nbApp").as<int>()){
    case 0:{
      Poisson p1 = Poisson(mesh, option(_name="alpha0").as<std::string>() );
      p1.run();
      e->add( "grad_u", p1.gradu );
      e->save();
      break;
    }
    case 1:{
      Eigen_Curl eig = Eigen_Curl(option(_name="needEigen").as<bool>(), mesh);
      for(int i=0; i<option(_name="solvereigen.nev").as<int>(); i++){
        e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig.g[i] );
        e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig.g0[i] );
        e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig.gradu[i] );
        e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig.modebis[i] );
      }
      e->save();
      break;
    }
    case 2:{
      Darcy d = Darcy(mesh,  option(_name="alpha0").as<std::string>() );
      d.run();
      e->add( "psiDiv", d.U );
      e->save();
      break;
    }
    default:
        break;
    }
    if ( Environment::worldComm().isMasterRank() )
      std::cout << "-----End-----" << std::endl;
}
