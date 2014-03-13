#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>

using namespace Feel;
using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "My app options" );
    myappOptions.add_options()
        ( "rayon", po::value<double>()->default_value( 0.05 ), "rayon" )
        ( "vitesse", po::value<double>()->default_value( 0.015 ), "vitesse moyenne d'entree" )
        ;
    return myappOptions.add( feel_options() ); // Add the default feel options to your list
}

int main(int argc, char**argv )
{
    // initialize feel++
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="po_psi0",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com"));

    typedef FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<1, Scalar>, Lagrange<0, Scalar> > > space_type;

    double rayon = option(_name="rayon").template as<double>();
    double vitesse = option( _name="vitesse").template as<double>();

    auto alpha0 = 2. * vitesse * (1. - (Px() * Px() + Pz() * Pz()) / (rayon * rayon));

    // create mesh
    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    auto Vh = space_type::New( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.template element<0>() ;
    auto lambda = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto nu = V.template element<1>() ;

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=gradt(u)*trans(grad(v)) + id( v )*idt( lambda ) + idt( u )*id( nu ) );

    auto l = form1( _test=Vh );
    l = integrate( _range=markedfaces(mesh, 1), //entree
                   _expr=alpha0*id(v) );
    l += integrate( _range=markedfaces(mesh, 2), //sortie
                    _expr=-alpha0*id(v) );
    l += integrate( _range=markedfaces(mesh, 3), //tour
                    _expr=cst(0.) );

    a.solve( _rhs=l, _solution=U );

    auto Xh = Pchv<1>( mesh );
    auto gradu = Xh->element();
    gradu = vf::project( _space=Xh, _range=elements( mesh ),
                         _expr=trans(gradv( U.template element<0>() )) );

    auto e = exporter( _mesh=mesh );
    e->add( "u", U.template element<0>() );
    e->add( "grad_u", gradu );
    e->save();
}
