/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description gridoptions( "Eigen options" );
    gridoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "kappa", po::value<double>()->default_value( 1 ), "coefficient" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ;
    return gridoptions.add( Feel::feel_options() );
}

template<int Dim>
class Eigen
:
public Simget
{
    typedef Simget super;
public:
    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;

    /**
     * Constructor
     */
    Eigen()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        eigen( SolverEigen<value_type>::build( this->vm() ) )
        {
        }

    void run();
private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
    std::vector<int> flags;
    boost::shared_ptr<SolverEigen<value_type> > eigen;
}; // Eigen

template<int Dim> const uint16_type Eigen<Dim>::Order;

template<int Dim>
void
Eigen<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Eigen<" << Dim << ">\n";
    Environment::changeRepository( boost::format( "po/%1%/%2%D-P%3%/h_%4%/" )
                                  % this->about().appName()
                                  % Dim
                                  % Order
                                  % meshSize );

    auto mesh = loadMesh(_mesh = new mesh_type );

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ), gradt(u)*trans(grad(v)) );

    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), idt( u )*id( v ) );


    int maxit = option(_name="solvereigen-maxiter").template as<int>();
    int tol = option(_name="solvereigen-tol").template as<double>();

    int nev = option(_name="solvereigen-nev").template as<int>();

    int ncv = option(_name="solvereigen-ncv").template as<int>();;

    double eigen_real, eigen_imag;

    SolverEigen<double>::eigenmodes_type modes;

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "nev= " << nev <<std::endl;
        std::cout << "ncv= " << ncv <<std::endl;
    }

    modes=
    eigs( _matrixA=a.matrixPtr(),
         _matrixB=b.matrixPtr(),
         _nev=nev,
         _ncv=ncv,
         _transform=SINVERT,
         _spectrum=SMALLEST_MAGNITUDE,
         _verbose = true );

    auto femodes = std::vector<decltype( Xh->element() )>( modes.size(), Xh->element() );

    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";

        int i = 0;
        BOOST_FOREACH( auto mode, modes )
        {
            std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second.get<1>() << ")\n";
            femodes[i++] = *mode.second.get<2>();
        }
    }

    auto exporter =  export_type::New( this->vm(),
                                      ( boost::format( "%1%-%2%-%3%" )
                                       % this->about().appName()
                                       % shape
                                       % Dim ).str() ) ;

    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );

        int i = 0;
        BOOST_FOREACH( auto mode, femodes )
        {
            exporter->step( 0 )->add( ( boost::format( "mode-%1%" ) % i++ ).str(), mode );
        }

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }

}

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                    _desc=makeOptions(),
                    _about=about(_name="po_eigen",
                                 _author="",
                                 _email="") );

    Application app;

    app.add( new Eigen<2>() );
    app.run();
}





