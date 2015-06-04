#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

typedef Mesh<Simplex<3> > meshdim;
typedef bases<Lagrange<1,Scalar> > baseP1;
typedef bases<Lagrange<1,Scalar>, Lagrange<0,Scalar> > baseP1P0;
typedef FunctionSpace<meshdim, baseP1> FEspaceP1;
typedef FunctionSpace<meshdim, baseP1P0> FEspaceP1P0;



int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="AAA",
                                   _author="BS-RH",
                                   _email="benjamin.surowiec@plasticomnium.com"));

    auto Th = loadMesh(new meshdim);

    auto XVh2 = FEspaceP1P0::New( Th );

    auto U = XVh2->element();
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto V = XVh2->element();
    auto v = V.element<0>();
    auto q = V.element<1>();

    auto a = form2( _trial=XVh2, _test=XVh2 );

    a = integrate( _range=elements(Th),
                   _expr=inner(gradt(u),grad(v)) );

    a += integrate( _range=elements(Th),
                    _expr=idt(u)*id(q) );

    a += integrate( _range=elements(Th),
                    _expr=id(v)*idt(p) );

    auto f = expr(soption("functions.f"));
    auto h_s = soption("functions.h");
    auto vars = Symbols{ "x", "y", "radius", "meanvel" };
    auto h_e = parse( h_s, vars );
    auto alpha0 = expr( h_e, vars );
    alpha0.setParameterValues( {
            { "radius", 0.5 },
                { "meanvel", 1 } } );

    auto d = integrate(markedfaces(Th, 1), -alpha0).evaluate()(0,0);
    d += integrate(markedfaces(Th, 2), alpha0).evaluate()(0,0);
    std::cout << " alpha0 : " << h_s << "\nint : " << d << std::endl;

    auto l = form1( _test=XVh2);

    l = integrate( _range=markedfaces(Th, 1),
                   _expr=-alpha0*id(v) );

    l = integrate( _range=markedfaces(Th, 2),
                   _expr=alpha0*id(v) );

    a.solve( _rhs=l,
             _solution=U );

    auto e = exporter(Th);
    e->add("u", u);
    e->add("p", p);
    e->save();

    return 0;
}
