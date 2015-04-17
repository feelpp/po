#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_stokesmixed",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef Mesh<Simplex<3 > > mesh_type;

    typedef bases<RaviartThomas<0>, Nedelec<0,NedelecKind::NED1>, Lagrange<1,Scalar > > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;

    // create the mesh
    auto mesh = loadMesh(_mesh=new mesh_type );

    auto Xh = space_type::New( mesh );

    auto U = Xh->element();

    auto u = U.element<0>();
    auto w = U.element<1>();
    auto p = U.element<2>();
  
    auto a = form2( _test=Xh, _trial=Xh );
    a = integrate( elements(mesh),
		   trans(idt(w))*id(w) - trans(idt(u))*curl(w)
		   + trans(curlt(w))*id(u) - idt(p)*div(u)
		   + divt(u)*id(p) );

    auto f = expr<3,1>(soption("functions.f"));
    auto sigma = expr<3,1>(soption("functions.sigma"));
    auto pi = expr(soption("functions.pi"));

    auto l = form1( _test=Xh );
    l = integrate( elements(mesh),
    		   trans(f)*id(u) );
    // l += integrate( boundaryfaces(mesh),
    // 		    trans(sigma)*cross(id(w),N())
    // 		    - pi*inner(id(u),N()) );

    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u, _expr=cst(0.) );
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=w, 
	     _expr=vec(cst(0.),cst(0.),cst(0.)) );
	     

    a.solve( _rhs=l, _solution=U );

    auto e = exporter( mesh );
    e->add( "u", u );
    e->add( "w", w );
    e->add( "p", p );
    e->save();
}
