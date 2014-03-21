#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

using namespace Feel;
using namespace Feel::vf;

class Psi0Div
    :
public Application
{
    typedef Application super;

public:

    typedef double value_type;
    typedef Simplex<3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    /// [space]
    typedef bases<RaviartThomas<1>, Lagrange<1, Scalar> > prod_basis_type;
    typedef FunctionSpace<mesh_type, prod_basis_type> prod_space_type;
    /// [space]

    void run();
};

/// [option]
inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "psi0Div options" );
    myappOptions.add_options()
        ( "rayon", po::value<double>()->default_value( 0.05 ), "rayon du cylindre" )
        ( "vitesse", po::value<double>()->default_value( 0.015 ), "vitesse moyenne d'entree" )
        ( "profil", po::value<std::string>()->default_value( "2. * vitesse * (1. - (x * x + y * y) / (rayon * rayon)):x:y:vitesse:rayon" ), "alpha0" );
    return myappOptions.add( feel_options() );
}
/// [option]


void
Psi0Div::run()
{
    /// [mesh]
    auto mesh = loadMesh(_mesh = new mesh_type );
    /// [mesh]

    /// [alpha0]
    double rayon = option(_name="rayon").template as<double>();
    double vitesse = option( _name="vitesse").template as<double>();
    auto alpha0 = 2. * vitesse * (1. - (Px() * Px() + Py() * Py()) / (rayon * rayon));
    /// [alpha0]

    /// [element]
    auto Xh = prod_space_type::New( mesh );
    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(v,q)" );

    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "v" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "q" );
    /// [element]

    /// [rhs]
    auto l = form1( _test=Xh );
    l = integrate( _range=markedfaces(mesh, 1), //entree
                   _expr=-alpha0*id(q) );
    l += integrate( _range=markedfaces(mesh, 2), //sortie
                    _expr=alpha0*id(q) );
    l += integrate( _range=markedfaces(mesh, 3), //tour
                    _expr=cst(0.)*id(q) );
    /// [rhs]

    /// [bilinear]
    auto a = form2( _test=Xh, _trial=Xh );
    a += integrate( _range=elements(mesh),
                    _expr = -trans( idt( u )) * id( v )
                    + trans(gradt( p )) * id( v )
                    + trans(idt( u )) * grad( q ) );

    a.solve( _rhs=l, _solution=U );
    /// [bilinear]
    auto e = exporter( _mesh=mesh );
    e->add( "grad_u", U );
    e->save();

}

int
main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                    _desc=makeOptions(),
                    _about=about(_name="po_psi0Div",
                                 _author="Romain Hild",
                                 _email="romain.hild@plasticomnium.com") );

    Psi0Div app;
    app.run();
}

