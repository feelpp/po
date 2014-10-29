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
                     _about=about(_name=
                                  "po_spheromak",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );

    auto Vh = Pchv<1>( mesh );
    std::cout << "dim Vh : " << Vh->nDof() << std::endl;

    auto u = Vh->element();

    auto eta1 = sqrt(Px()*Px()+Py()*Py()+Pz()*Pz());
    auto eta2 = acos(Pz()/eta1);
    auto eta3 = atan(Py()/Px());

    auto l = doption("parameters.l");

    auto uext = vec( 2*l/eta1*(l/eta1*sin(eta1/l)-cos(eta1/l))*cos(eta2),
		     -1/eta1*(l/eta1*cos(eta1/l)-l*l/(eta1*eta1)*sin(eta1/l)+sin(eta1/l))*sin(eta2),
		     1/eta1*(l/eta1*sin(eta1/l)-cos(eta1/l))*sin(eta2) );
    
    auto uint = vec( -eta1/(12*l*l*l)*cos(eta2),
		     (2/(3*l)-2*eta1*eta1/(15*l*l*l))*sin(eta2),
		     (-eta1/(3*l*l)-eta1*eta1*eta1/(30*l*l*l*l))*sin(eta2) );

    u = project( _range=markedelements(mesh, "ext"), _space=Vh,
		 _expr=uext );
    u += project( _range=markedelements(mesh, "int"), _space=Vh,
		  _expr=uint );

    auto normal = normL2( _range=boundaryfaces(mesh), _expr=inner(idv(u),N()) );
    auto err = normL2( _range=elements(mesh), _expr=l*idv(u)-curlv(u) );

    std::cout << "normal : " << normal << "\nerreur : " << err << std::endl;

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();

    return 0;
}
