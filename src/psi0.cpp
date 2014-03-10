#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>

using namespace Feel;

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
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="po_psi0",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com"));
    // create mesh
    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    // function space
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    double rayon = option(_name="r").template as<double>();
    double vitesse = option( _name="v").template as<double>();

    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );

    // right hand side
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v) );

    // apply the boundary condition
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
            _expr=constant(0.) );

    //l = integrate(_range=boundaryfaces(mesh),
                  //_expr=(2. * vitesse * (1. - (Px()*Px() + Pz()*Pz())/(rayon*rayon)))*id(v) );


    //a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
                  //_expr= 2 * vitesseMoyIn * (1 - (Px()*Px() + Pz()*Pz())/(ray*ray)) );

    // solve the equation a(u,v) = l(v)
    a.solve(_rhs=l,_solution=u);

    auto gu = gradv(u);
    // export results
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    //e->add( "gu", gu );
    e->save();
}
