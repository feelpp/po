//#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>

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

    //auto mesh = unitSphere<1>();
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );

    auto Vh = Pch<2>( mesh );

    auto u = Vh->element();
    auto v = Vh->element();
    auto w = Vh->element();

    auto eta1 = sqrt(Px()*Px()+Py()*Py()+Pz()*Pz());
    auto eta2 = acos(Pz()/eta1);

    auto l = doption("parameters.l");

    u = project( _range=markedelements(mesh,"int"), _space=Vh,
		 _expr=(2/(3*l)-Px()*Px()/(15*l*l*l))*cos(eta2) );
    v = project( _range=markedelements(mesh,"int"), _space=Vh,
		 _expr=(-2/(3*l)+2*Px()*Px()/(15/l*l*l))*sin(eta2) );
    w = project( _range=markedelements(mesh,"int"), _space=Vh,
		 _expr=(Px()/(3*l*l)-Px()*Px()*Px()/(30*l*l*l*l))*sin(eta2) );

    u = project( _range=markedelements(mesh,"ext"), _space=Vh,
    		 _expr=2*l/(eta1*eta1)*(l/eta1*sin(eta1/l)-cos(eta1/l))*cos(eta2) );
    v = project( _range=markedelements(mesh,"ext"), _space=Vh,
    		 _expr=-1/eta1*(l/eta1*cos(eta1/l)-l*l/(eta1*eta1)*sin(eta1/l)+sin(eta1/l))*sin(eta2) );
    w = project( _range=markedelements(mesh,"ext"), _space=Vh,
    		 _expr=1/eta1*(l/eta1*sin(eta1/l)-cos(eta1/l))*sin(eta2) );

    auto normal = normL2( _range=boundaryfaces(mesh), _expr=idv(u) );
    if ( Environment::isMasterRank() )
    	std::cout << "u|Gamma = " << normal << std::endl;


    auto Xh = Pchv<2>( mesh );
    auto V = Xh->element();

    V = vf::project( _range=elements(mesh), _space=Xh,
    		 _expr=vec(idv(u)*sin(idv(v))*cos(idv(w)),
    			   idv(u)*sin(idv(v))*sin(idv(w)),
    			   idv(u)*cos(idv(v)) ) );

    auto n = normL2( _range=boundaryfaces(mesh), _expr=inner(idv(V),N()) );
    auto err = normL2( _range=elements(mesh), _expr=l*idv(V)-curlv(V) );
    if ( Environment::isMasterRank() )
	std::cout << "normal : " << n << "\nerr : " << err << std::endl;

    auto e = exporter( mesh );
    e->add( "u", u );
    e->add( "v", v );
    e->add( "w", w );
    e->add( "V", V );
    e->save();

}
