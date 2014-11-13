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

    if ( Environment::isMasterRank() )
        std::cout << "[dof] number of dof: " << Vh->nDof() << std::endl;

    auto ur = Vh->element("0");
    auto ui = Vh->element("0");

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
        ur += project(_range=elements(mesh), _space=Vh,
                      _expr=a[i]*cos(inner(P(),k[i])) - b[i]*sin(inner(P(),k[i])) );
        ui += project(_range=elements(mesh), _space=Vh,
                      _expr=a[i]*sin(inner(P(),k[i])) + b[i]*cos(inner(P(),k[i])) );
    }


    auto errR = normL2(_range=elements(mesh), _expr=curlv(ur)-std::sqrt(14)*idv(ur) );
    auto errI = normL2(_range=elements(mesh), _expr=curlv(ui)-std::sqrt(14)*idv(ui) );
    auto norR = normL2(_range=elements(mesh), _expr=inner(idv(ur),N()) );
    auto norI = normL2(_range=elements(mesh), _expr=inner(idv(ui),N()) );
    auto divR = normL2(_range=elements(mesh), _expr=divv(ur) );
    auto divI = normL2(_range=elements(mesh), _expr=divv(ui) );

    if ( Environment::isMasterRank() ) {
        std::cout << "nR : " << norR << "\nnI : " << norI << std::endl;
        std::cout << "eR : " << errR << "\neI : " << errI << std::endl;
        std::cout << "dR : " << divR << "\ndI : " << divI << std::endl;
    }


    auto e = exporter( mesh );
    e->add( "Re(u)", ur);
    e->add( "Im(u)", ui);
    e->save();

}
