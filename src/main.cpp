/** \file main.cpp
    \brief Main function to resolves the problem
*/

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include "psi0.hpp"
#include "eigenprob.hpp"
#include "spectralproblem.hpp"
#include "testCurl.hpp"

using namespace Feel;
using namespace Eigen;

/** \brief Add the option to the application

    Here is the list of the options :
    - needP0 : need to compute psi0
    - computeP0 : compute psi0 or load it
    - radius : cylinder's radius
    - speed : average speed
    - alpha0 : alpha0, depends on x,y,radius
    - needEigen : need the eigen modes
    - computeEigen : compute the eigen modes or load them
    - useCurl : use curl or grad form
    - usePresDiv : use a pressure term in div form
    - usePresGrad : use a pressure term in grad form
    - useDiric : use Dirichlet condition
    - bccurln : need boundary condition curl g.n
    - bcn : need boundary condition1 g.n
    - divdiv : need divdiv term
    - needDecomp : need to decompose the eigen modes
    - computeDecomp : compute the decomposition of the modes or load them
    - meanPsi : psi average
    - needDebug : debug
    - needPS : need to run the spectral problem
    - computeRijk : compute or load Rijk
    - computeRiak : compute or load Riak
    - computeRfk : compute or load Rfk
    - f : f
    - computeRpk : compute or load Rpk
    - nu : viscosity
    - alpha2 : alpha2, depends on speed and radius
    - alpha1 : alpha1, (0.)
    - testCurl : test curl
    - t : pol order 2
    - t1 : curl(t)
    - t2 : curl2(t)
*/
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

        ( "alpha1", po::value<std::string>()->default_value( "0." ), "alpha1, (0.)" )

        ( "testCurl", po::value<bool>()->default_value( true ), "test curl")
        ( "t", po::value<std::string>()->default_value( "{x*x+2*y*y+3*z*z,2*x*x+3*y*y+z*z,3*x*x+y*y+2*z*z}:x:y:z" ), "pol order 2" )
        ( "t1", po::value<std::string>()->default_value( "{2*y-2*z,6*z-6*x,4*x-4*y}:x:y:z" ), "curl(t)" )
        ( "t2", po::value<std::string>()->default_value( "{-10,-6,-8}:x:y:z" ), "curl2(t)" );
    return myappOptions;
}

/// Add the options of Feel++ to the application
po::options_description
makeLibOptions()
{
    po::options_description libOptions( "Lib options" );
    libOptions.add( backend_options( "psi0" ) ).add( backend_options( "gradpsi0" ) ); // Poisson
    libOptions.add( backend_options( "gi0" ) ).add( backend_options( "psi" ) ).add( backend_options( "gradpsi" ) ).add( backend_options( "curl" ) ).add( backend_options( "curl2") ); // Eigen_Curl
    return libOptions.add( feel_options() );
}

/// Add the info to the application
AboutData
makeAbout()
{
    Feel::AboutData about( "po_app", "po_app" );
    about.addAuthor( "Romain Hild", "", "romain.hild@plasticomnium.com", "" );
    return about;
}

/** \brief Load the mesh

    If the option computeEigen is set to true and if the option gmsh.filename's extension is .msh, rebuid the mesh with the correct number of partition\n
    else, load the mesh with the extension .msh
    */
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

#ifndef PO_TRAVIS
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

/** \brief Main function of the application

    Put in place the environment for the application\n
    Load the mesh

    If the option needP0 is set to true, run Psi0\n
    If the option needEigen is set to true, run EigenProb\n
    If the option needPS is set to true, run Psi0, EigenProb and SpectralProblem
    If the option testCurl is set to true, run TestCurl

    Export the elements needed
*/
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
    TestCurl tc = TestCurl( mesh );

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

    if( boption(_name="testCurl") ){
        tc.run();
        e->add( "t", tc.test[0] );
        e->add( "curl(t)", tc.test[1] );
        e->add( "curl2(t)", tc.test[2] );
        e->add( "projCurl", tc.test[3] );
        e->add( "projCurl2", tc.test[4] );
        e->add( "sysCurl", tc.test[5] );
        e->add( "sysCurl2", tc.test[6] );
    }

    if( boption(_name="needPS") || boption(_name="needEigen") || boption(_name="needP0") || boption("testCurl") )
        e->save();


    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- End -----" << std::endl;


    return 0;
}

/** \mainpage Non Standard Navier-Stokes equations for large scale computational fluid dynamics in Feel++.

 Notations :
 \f{align*}
 gradient(v)&=(\partial_x v, \partial_y v, \partial_z v)=\grad v\\
 gradient(\mathbf{v})&=\begin{pmatrix}
\partial_x v_x & \partial_y v_x & \partial_z v_x\\
\partial_x v_y & \partial_y v_y & \partial_z v_y\\
\partial_x v_z & \partial_y v_z & \partial_z v_z
\end{pmatrix}=\grad\mathbf*{v}\\
divergence(\mathbf{v})&==\frac{\partial v_x}{\partial x}+\frac{\partial v_y}{\partial y}+\frac{\partial v_z}{\partial z}=\div \mathbf{v}\\
curl(\mathbf{v})&=\begin{pmatrix}
\partial_y v_z - \partial_z v_y\\
\partial_z v_x - \partial_x v_z\\
\partial_x v_y - \partial_y v_x
\end{pmatrix}=\rot \mathbf{v}\\
curl(curl(\mathbf{v}))&=\rott \mathbf{v}\\
H^1(\Omega) &= \{v \in L^2(\Omega)\;|\; \grad v\in L^2(\Omega)\}\\
H^1_0(\Omega) &= \{v \in H^1(\Omega)\; |\; v\restr = 0\}\\
H(\mathrm{div}) &= \{\mathbf{v} \in [L^2(\Omega)]^3\; |\; \div\mathbf{v} \in L^2(\Omega) \}\\
H(\mathrm{rot}) &= \{\mathbf{v} \in [L^2(\Omega)]^3\; |\; \rot\mathbf{v} \in L^2(\Omega) \}\\
L^2_\sigma(\Omega) &= \{\mathbf{v} \in [L^2(\Omega)]^3\; |\; \div \mathbf{v} = 0\text{ et }\mathbf{v}\cdot \mathbf{n}\restr = 0 \}\\
D^1(\Omega) &= \{\mathbf{v} \in [H^1(\Omega)]^3\cap L^2_\sigma(\Omega)\; |\; (\rot \mathbf{v}\cdot \mathbf{n})\restr = 0  \}
 \f}

We are looking for the couple \f$(\mathbf{v},p)\f$, respectively the speed and the pressure which are solutions of the incompressible Navier-Stokes equation in \f$Q_T=\Omega\times[0,T]\f$, where \f$\Omega\f$ is an open set in \f$\R^3\f$ and \f$\partial\Omega\f$ is its boundary.
\f{pb}
We are looking for \f$(\mathbf{v},p)\f$ such as :
\begin{equation*}
\left\{\begin{aligned}
&\frac{\partial\mathbf{v}}{\partial t} + (\rot \mathbf{v})\times \mathbf{v} + \grad q + \frac{1}{Re}\rott \mathbf{v}=\mathbf{f} = 0\\

\end{aligned}\right.
\end{equation*}
where \f$q=\frac{\mathbf{v}\cdot\mathbf{v}}{2}+p\f$
\f}
*/
