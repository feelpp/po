#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_test",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(new Mesh<Simplex<3> >);

#define LM
#ifndef LM
    typedef FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<1,Scalar> > > space_type;
#else
    typedef FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<1,Scalar>, Lagrange<0,Scalar> > > space_type;
#endif

    auto Vh = space_type::New( mesh );

#ifndef LM
    auto u = Vh->element();
    auto v = Vh->element();
#else
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto lambda = U.element<1>();
    auto v = V.element<0>();
    auto nu = V.element<1>();
#endif

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=inner(gradt(u),grad(v)) );
#ifdef LM
    a += integrate( _range=elements(mesh),
                    _expr=id(v)*idt(lambda) + idt(u)*id(nu) );
#endif

    auto f = expr(soption("functions.f"));
    auto g_s = soption("functions.g");
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", 0.5 },
                { "speed", 1 } } );

    auto d = integrate(markedfaces(mesh, 1), -g).evaluate()(0,0);
    d += integrate(markedfaces(mesh, 2), g).evaluate()(0,0);
    std::cout << " g : " << g_s << "\nint : " << d << std::endl;

    auto l = form1( _test=Vh );
    // l = integrate( _range=elements(mesh),
    //                _expr=f*id(v) );
    l = integrate( _range=markedfaces(mesh, 1),
                   _expr=-g*id(v) );
    l = integrate( _range=markedfaces(mesh, 2),
                   _expr=g*id(v) );
    a.solve(_rhs=l,
#ifndef LM
            _solution=u
#else
            _solution=U
#endif
            );
    auto e = exporter(mesh);
    e->add("u", u);
    e->save();

    return 0;
}
