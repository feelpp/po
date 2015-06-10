#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
//#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pch.hpp>
//#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

#define DIM 2

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_test",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(new Mesh<Simplex<DIM> >);

    // typedef FunctionSpace<Mesh<Simplex<DIM> >, bases<Lagrange<2, Vectorial, Continuous>, Lagrange<2, Scalar> > > space_type;
    // typedef FunctionSpace<Mesh<Simplex<DIM> >, bases<Lagrange<2, Vectorial, Discontinuous>, Lagrange<2, Scalar> > > space_type;
    typedef FunctionSpace<Mesh<Simplex<DIM> >, bases<RaviartThomas<0>, Lagrange<2, Scalar, Continuous> > > space_type;
    // typedef FunctionSpace<Mesh<Simplex<DIM> >, bases<RaviartThomas<0> > > space_type;
    auto Sh = space_type::New( mesh );

    // auto f = expr<DIM,1>(soption("functions.f"));
    // auto u1 = Sh->element();
    // u1.on(_range=elements(mesh), _expr=f);
    // auto err = normL2(elements(mesh), idv(u1)-f);
    // err += normL2(elements(mesh), divv(u1)-div(f));
    // std::cout << "elements : " << err << std::endl;
    // auto u5 = Sh->element();
    // u5.on(_range=boundaryfaces(mesh), _expr=f);
    // err = normL2(boundaryfaces(mesh), trans(idv(u5))*N()-trans(f)*N());
    // std::cout << "boundaryfaces : " << err << std::endl;


    auto g = expr(soption("functions.g"));
    auto h = grad<DIM>(g);
    auto lhs = laplacian(g);
    auto intf = integrate(elements(mesh), lhs).evaluate()(0,0);
    auto intgn = integrate(boundaryfaces(mesh), h*N()).evaluate()(0,0);
    std::cout << g << std::endl
              << h << std::endl
              << lhs << std::endl
              << intf << std::endl
              << intgn << std::endl;

    auto U = Sh->element();
    auto V = Sh->element();
    auto u = U.template element<0>();
    auto v = V.template element<0>();
    auto p = U.template element<1>();
    auto q = V.template element<1>();

    auto a = form2( _trial=Sh, _test=Sh );
    a = integrate(elements(mesh), -trans(idt(u))*id(v));
    a += integrate(elements(mesh), -divt(u)*id(q));
    a += integrate(elements(mesh), -idt(p)*div(v));
    a += integrate(elements(mesh), -1./2*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q))));
    a += integrate(elements(mesh), -1./2*divt(u)*div(v));
    a += integrate(elements(mesh), -1./2*trans(curlt(u))*curl(v));

    auto l = form1( _test=Sh );
    l = integrate( elements(mesh), -lhs*id(q));
    l += integrate( elements(mesh), -1./2*lhs*div(v));

    // a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=trans(h));

    a.solve(_rhs=l, _solution=U);

    auto err = normL2( elements(mesh), idv(u) - trans(gradv(p)));
    auto errU = normL2( elements(mesh), idv(u) - trans(h));
    auto errP = normL2( elements(mesh), idv(p) - g);
    auto divu = normL2( elements(mesh), divv(u) - lhs);
    auto lapp = normL2(elements(mesh), laplacianv(p) - lhs);
    auto errUn = normL2( boundaryfaces(mesh), trans(idv(u))*N() - h*N());
    auto un = normL2( boundaryfaces(mesh), trans(idv(u))*N());
    auto hn = normL2( boundaryfaces(mesh), trans(h)*N());
    auto errGpn = normL2( boundaryfaces(mesh), gradv(p)*N() - h*N());
    auto intun = integrate( boundaryfaces(mesh), trans(idv(u))*N()).evaluate()(0,0);
    if( Environment::isMasterRank())
        std::cout << "err : " << err << std::endl
                  << "errU : " << errU << std::endl
                  << "errP : " << errP << std::endl
                  << "lapp : " << lapp << std::endl
                  << "divu : " << divu << std::endl
                  << "errUn : " << errUn << std::endl
                  << "un : " << un << std::endl
                  << "hn : " << hn << std::endl
                  << "intun : " << intun << std::endl
                  << "gpn : " << errGpn << std::endl;


    auto e = exporter(mesh);
    e->add("u", u);
    e->add("p", p);
    e->save();

    return 0;
}
