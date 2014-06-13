#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "poisson.h"
#include "eigen_curl.h"
#include "eigenlap.h"
#include "darcy.h"
#include "spectralproblem.h"
#include "eigenlapZ.h"

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
        ( "gamma", po::value<double>()->default_value( 20. ), "penalisation for curl g.n (eigenlap)" )
        ( "bc", po::value<double>()->default_value( 0. ), "boundary condition for eigenlapZ" )
        ( "needEigen", po::value<bool>()->default_value( true ), "need to compute the eigen modes or to load them" )
        ( "bccurln", po::value<bool>()->default_value( true ), "need boundary condition curl g.n (eigenlap)" )
        ( "bcn", po::value<bool>()->default_value( true ), "need boundary condition1 g.n (eigenlap)" )
        ( "divdiv", po::value<bool>()->default_value( true ), "need divdiv term" )
        ( "needDecomp", po::value<bool>()->default_value( true ), "need to decompose the eigen modes" )
        ( "needRelev", po::value<bool>()->default_value( false ), "need a" )
        ( "needDebug", po::value<bool>()->default_value( false ), "debug" )
        ( "testBessel", po::value<bool>()->default_value( false ), "test the third composante using Bessel functions" )
        ( "k", po::value<double>()->default_value( 0.0 ), "k-th Bessel function" )
        ( "j", po::value<int>()->default_value( 1 ), "j-th Bessel's root" )
        ( "m", po::value<int>()->default_value( 1 ), "m-th triple index" )
        ( "nbApp", po::value<int>()->default_value( 0 ), "app to launch (0:Poisson, 1:eigencurl, 2:darcy 3:spectral problem, 4:eigenlap 5:eigenlapZ)" );
    return myappOptions;
}

po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ); // Poisson
    libOptions.add( backend_options( "gi0" ) ).add( backend_options( "psi" ) ).add( backend_options( "gradpsi" ) ).add( backend_options( "curl" ) ); // Eigen_Curl
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

    Poisson p0 = Poisson( mesh, soption( _name="alpha0" ) );
    Eigen_Curl eig = Eigen_Curl( mesh );
    EigenLap eig2 = EigenLap( mesh );
    EigenLapZ eig3 = EigenLapZ( mesh );
    Darcy d = Darcy( mesh, soption( _name="alpha0" ) );
    SpectralProblem sp = SpectralProblem( mesh );

    auto e = exporter( _mesh=mesh );

    switch(ioption(_name="nbApp")){
    case 0:{
        p0.run();
        e->add( "grad_u", p0.gradu );
        e->save();
        break;
    }
    case 1:{
        eig.run();
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++){
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig.g[i] );
            e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig.g0[i] );
            e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig.psi[i] );
            e->add( ( boost::format( "gradpsi-%1%" ) % i ).str(), eig.gradu[i] );
            e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig.modebis[i] );
        }
        e->save();
        break;
    }
    case 2:{
        d.run();
        e->add( "psiDiv", d.U );
        e->save();
        break;
    }
    case 3:{
        if( boption( _name="needRelev" ) ){
            p0.run();
            e->add( "a", p0.gradu );
        }
        eig3.run();
        sp.init( eig3.g, eig3.psi, eig3.lambda, p0.gradu );
        sp.run();
        e->add( "u", sp.u );
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++)
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig3.g[i] );
        e->save();
        break;
    }
    case 4:{
        eig2.run();
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++){
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig2.g[i] );
            if( boption( _name="needDecomp") ){
                e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig2.g0[i] );
                e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig2.psi[i] );
                e->add( ( boost::format( "gradpsi-%1%" ) % i ).str(), eig2.gradu[i] );
                e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig2.modebis[i] );
            }
        }
        e->save();
        break;
    }
    case 5:{
        eig3.run();
        for(int i=0; i<ioption(_name="solvereigen.nev"); i++){
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), eig3.g[i] );
            if( boption( _name="needDecomp") ){
                e->add( ( boost::format( "g0-%1%" ) % i ).str(), eig3.g0[i] );
                e->add( ( boost::format( "psi-%1%" ) % i ).str(), eig3.psi[i] );
                e->add( ( boost::format( "gradpsi-%1%" ) % i ).str(), eig3.gradu[i] );
                e->add( ( boost::format( "modebis-%1%" ) % i ).str(), eig3.modebis[i] );
            }
        }
        e->save();
        break;
    }
    default:
        break;
    }
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End -----" << std::endl;
}
