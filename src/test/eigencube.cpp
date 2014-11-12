#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feeldiscr/pchv.hpp>
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
    auto u = Vh->element();

    auto a1 = vec(cst(0.0113),cst(0.),cst(-0.0219));
    auto a2 = vec(cst(-0.0004),cst(0.),cst(0.0219));
    auto a3 = vec(cst(0.0004),cst(0.),cst(-0.0219));
    auto a4 = vec(cst(-0.0113),cst(0.),cst(0.0219));

    auto b1 = vec(cst(0.0175),cst(-0.0089),cst(0.0111));
    auto b2 = vec(cst(-0.0175),cst(0.0320),cst(-0.0111));
    auto b3 = vec(cst(-0.0175),cst(-0.0089),cst(0.0034));
    auto b4 = vec(cst(-0.0175),cst(0.0320),cst(-0.0034));

    auto k1 = vec(cst(2.),cst(3.),cst(1.));
    auto k2 = vec(cst(2.),cst(3.),cst(-1.));
    auto k3 = vec(cst(2.),cst(-3.),cst(1.));
    auto k4 = vec(cst(-2.),cst(3.),cst(1.));

    u = project(_range=elements(mesh), _space=Vh,
                _expr=a1*cos(inner(P(),k1))-b1*sin(inner(P(),k1)) );
    u += project(_range=elements(mesh), _space=Vh,
                 _expr=a2*cos(inner(P(),k2))-b2*sin(inner(P(),k2)) );
    u += project(_range=elements(mesh), _space=Vh,
                 _expr=a3*cos(inner(P(),k3))-b3*sin(inner(P(),k3)) );
    u += project(_range=elements(mesh), _space=Vh,
                 _expr=a4*cos(inner(P(),k4))-b4*sin(inner(P(),k4)) );

    auto err = normL2(_range=elements(mesh), _expr=curlv(u)+std::sqrt(14)*idv(u) );
    auto nor = normL2(_range=elements(mesh), _expr=inner(idv(u),N()) );

    if ( Environment::isMasterRank() )
	std::cout << "n : " << nor << "\ne : " << err << std::endl;

    auto e = exporter( mesh );
    e->add( "u", u);
    e->save();

}
