/*
  laplace(psi) = 0
  grad(psi).n = alpha

  int_O grad(u)*grad(v) + u*nu + v*lambda = int_pO alpha*v
  use of Lagrange multipliers int_0 u = 0
 */
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>

using namespace Feel;
using namespace Feel::vf;

// [option]
inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "psi0 options" );
    myappOptions.add_options()
        ( "rayon", po::value<double>()->default_value( 0.05 ), "rayon" )
        ( "vitesse", po::value<double>()->default_value( 0.015 ), "vitesse moyenne d'entree" )
        ( "profil", po::value<std::string>()->default_value( "2. * vitesse * (1. - (x*x + y*y) / (rayon * rayon))" ), "alpha0, depend de x,y,rayon,vitesse" );
    return myappOptions; // Add the default feel options to your list
}
// [option]

int main(int argc, char**argv )
{
    // initialize feel++
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=feel_options().add( backend_options( "psi0" ) ).add ( backend_options( "gradpsi0" ) ),
                     _about=about(_name="po_psi0",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com"));

    // [space]
    typedef FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<1, Scalar>, Lagrange<0, Scalar> > > space_type;
    // [space]

    // [alpha0]
    auto alpha0_s = option( _name="profil" ).as<std::string>();
    auto vars = Symbols{ "x", "y", "vitesse", "rayon" };
    auto alpha0_e = parse( alpha0_s, vars );
    auto alpha0 = expr( alpha0_e, vars );
    alpha0.setParameterValues( {
            { "rayon", option( _name="rayon" ).template as<double>() },
                { "vitesse", option( _name="vitesse" ).template as<double>() },
                    {"test", 3} } );
    // [alpha0]

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    // [element]
    auto Vh = space_type::New( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.template element<0>() ;
    auto lambda = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto nu = V.template element<1>() ;
    // [element]

    // [bilinear]
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=gradt(u)*trans(grad(v)) + id( v )*idt( lambda ) + idt( u )*id( nu ) );
    // [bilinear]

    // [rhs]
    auto l = form1( _test=Vh );
    l = integrate( _range=markedfaces(mesh, 1), //entree
                   _expr=-alpha0*id(v) );
    l += integrate( _range=markedfaces(mesh, 2), //sortie
                    _expr=alpha0*id(v) );
    l += integrate( _range=markedfaces(mesh, 3), //tour
                    _expr=cst(0.)*id(v) );
    // [rhs]

    a.solve( _name="psi0", _rhs=l, _solution=U );

    // [gradpsi0]
    auto Xh = Pchv<1>( mesh );
    auto gradu = Xh->element();
    auto b = form2( _trial=Xh, _test=Xh );
    b = integrate( _range=elements(mesh), _expr=trans(idt(gradu))*id(gradu));
    auto f = form1( _test=Xh );
    f = integrate( _range=elements(mesh), _expr=gradv( U.template element<0>() )*id(gradu));
    // gradu is the L2 projection of grad(psi0) over Xh
    b.solve( _name="gradpsi0", _rhs=f, _solution=gradu );

    // [gradpsi0]

    auto e = exporter( _mesh=mesh );
    e->add( "u", U.template element<0>() );
    e->add( "grad_u", gradu );
    e->save();
}
