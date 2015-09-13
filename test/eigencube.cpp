#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_eigencube",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );

    auto Vh = Pchv<2>( mesh );
    auto Sh = Pch<1>( mesh );

    if ( Environment::isMasterRank() )
        std::cout << "[dof] number of dof: " << Vh->nDof() << std::endl;

    auto u = Vh->element("0");
    auto cu = Vh->element();
    auto curlAna = Vh->element("0");
    auto ccurlAna = Vh->element("0");
    auto divuAna = Sh->element("0");

    auto a1 = vec(cst(0.0113),cst(0.),cst(-0.0219));
    auto a2 = vec(cst(-0.0004),cst(0.),cst(0.0219));
    auto a3 = vec(cst(0.0004),cst(0.),cst(-0.0219));
    auto a4 = vec(cst(-0.0113),cst(0.),cst(0.0219));
    decltype(a1) a[4] = {a1,a2,a3,a4};

    auto b1 = vec(cst(0.0175),cst(-0.0089),cst(0.0111));
    auto b2 = vec(cst(-0.0175),cst(0.0320),cst(-0.0111));
    auto b3 = vec(cst(-0.0175),cst(0.0089),cst(-0.0034));
    auto b4 = vec(cst(-0.0175),cst(0.0320),cst(-0.0034));
    decltype(b1) b[4] = {b1,b2,b3,b4};

    auto k1 = vec(cst(2.),cst(3.),cst(1.));
    auto k2 = vec(cst(2.),cst(3.),cst(-1.));
    auto k3 = vec(cst(2.),cst(-3.),cst(1.));
    auto k4 = vec(cst(-2.),cst(3.),cst(1.));
    decltype(k1) k[4] = {k1,k2,k3,k4};


    for( int i = 0; i < 4; i++ ) {
        u += project(_range=elements(mesh), _space=Vh,
                      _expr=2*(a[i]*cos(inner(P(),k[i]))
			       - b[i]*sin(inner(P(),k[i]))) );

	divuAna += project(_range=elements(mesh), _space=Sh,
		      _expr=-2*(inner(k[i],a[i])*sin(inner(P(),k[i]))
				+ inner(k[i],b[i])*cos(inner(P(),k[i]))) );
	curlAna += project(_range=elements(mesh), _space=Vh,
		      _expr=2*(cross(a[i],k[i])*sin(inner(P(),k[i]))
			       + cross(b[i],k[i])*cos(inner(P(),k[i]))) );
	ccurlAna += project(_range=elements(mesh), _space=Vh,
		       _expr=2*(cross(cross(k[i],a[i]),k[i]) * cos(inner(P(),k[i]))
				+ cross(cross(b[i],k[i]),k[i]) * sin(inner(P(),k[i])) ) );
    }

    cu = vf::project(_range=elements(mesh), _space=Vh, _expr=curlv(u) );


    auto nor = normL2(_range=elements(mesh), _expr=inner(idv(u),N()) );

    auto div = normL2(_range=elements(mesh), _expr=divv(u) );
    // auto divAna = normL2(_range=elements(mesh), _expr=idv(divuAna) );
    // auto errD = normL2(_range=elements(mesh),
    // 		       _expr=divv(u)-idv(divuAna) );

    auto err = normL2(_range=elements(mesh),
		      _expr=curlv(u)+std::sqrt(14)*idv(u) );
    auto errC = normL2(_range=elements(mesh),
		       _expr=curlv(u)-idv(cu) );
    auto errC2 = normL2(_range=elements(mesh),
		       _expr=idv(cu)-idv(curlAna) );

    auto errCurl = normL2(_range=elements(mesh),
			  _expr=curlv(cu) - idv(ccurlAna) );
    auto err2 = normL2(_range=elements(mesh),
		       _expr=idv(ccurlAna) - 14*idv(u) );

    if ( Environment::isMasterRank() ) {
        std::cout << "normale : " << nor << std::endl;
        std::cout << "divergence : " << div << std::endl;
	// std::cout << "divergence analytique : " << divAna << std::endl;
        // std::cout << "err(divv-divAna) : " << errD << std::endl;
        std::cout << "err(curlv-eigen) : " << err << std::endl;
        std::cout << "err(cu-curlAna) : " << errC << std::endl;
        std::cout << "err(curlv-curlAna) : " << errC2 << std::endl;
        std::cout << "err(curlv(cu)-curl2Ana) : " << errCurl << std::endl;
        std::cout << "err(curlAna-14u) : " << err2 << std::endl;
    }


    auto e = exporter( mesh );
    e->add( "u", u);
    e->add( "curlu", cu);
    e->save();

}
