#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="po_testD1",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    
    auto f = expr<3,1>(soption("functions.f"));
    f.setParameterValues( {
	    { "lambda", doption( _name="parameters.l" ) } } );
    std::cout << f << std::endl;

    auto mesh = unitSphere();

    auto Vh = Pchv<1>( mesh );

    // auto eta1 = sqrt(Px()*Px()+Py()+Py()+Pz()+Pz());
    // auto eta2 = acos(Pz()/sqrt(Px()*Px()+Py()+Py()+Pz()+Pz()));
    // auto eta3 = atan(Py()/Px());

    // auto l = 4.493409458;

    // auto v = vec(2*l/(eta1*eta1)*(l/eta1*sin(eta1/l)-cos(eta1/l))*cos(eta2),
    // 		      -1/eta1*(l/eta1*cos(eta1/l)-l*l/(eta1*eta1)*sin(eta1/l)+sin(eta1/l))*sin(eta2),
    // 		      1/eta1*(l/eta1*sin(eta1/l)-cos(eta1/l))*sin(eta2));

    auto u = Vh->element(f);

    auto e = exporter( mesh );
    e->add("u",u);
    e->save();
    
    return 0;
}
