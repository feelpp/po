#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>

using namespace Feel;

#define FEELPP_DIM 3

template<typename F>
void
test( F Vh, std::string exporterName )
{
    auto u = Vh->element();
    // do stuff here with u
    auto e = exporter( _mesh=Vh->mesh(), _name=exporterName );
    e->add( "u", u );
    e->save();
}

template<typename F>
typename F::element_type::element_type
test2( F Vh )
{
    // auto alpha = expr(alphaS);
    auto u = Vh->element();
    auto a = form2( _test=Vh, _trial=Vh );
    // a = integrate( elements(Vh->mesh()), ...);
    auto l = form1( _test=Vh );
    // l = ...;
    // a .solve();
    // do stuff here with u
    return u;
}

int main( int argc, char** argv )
{
    typedef Simplex<FEELPP_DIM> convex_type;
    typedef Mesh<convex_type> mesh_type;

    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="templateExample", _author="R.H", _email="hild.romain@gmail.com" ) );

    auto mesh = loadMesh( _mesh=new mesh_type);
    auto Sh = Pch<1>( mesh );
    test<decltype(Sh)>( Sh, "pch1" );
    auto u = test2<decltype(Sh)>( Sh );

    auto Vh = Pchv<2>( mesh );
    auto v = test2<decltype(Vh)>( Vh );
    auto e = exporter( mesh, "pchv2" );
    e->add( "v", v );
    e->save();

    return 0;
}
