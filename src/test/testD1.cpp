#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
//#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
		     _desc=feel_options(),
                     _about=about(_name=
#if PO_CYL
				  "po_testD1_cyl",
#else
				  "po_testD1_car",
#endif
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );

    auto Vh = Pch<1>( mesh );
    std::cout << "dim Vh : " << Vh->nDof() << std::endl;

    auto u = Vh->element();
    auto v = Vh->element();

#if PO_CYL
    auto eta1 = Px();
    auto eta2 = Py();
    auto eta3 = Pz();
    auto Jacmat = mat<3,3> ( cos(eta2),-eta1*sin(eta2),cst(0.),
			     sin(eta2),eta1*cos(eta2),cst(0.),
			     cst(0.),cst(0.),cst(1.)            );
#else
    auto eta1 = sqrt(Px()*Px()+Py()*Py());
    //auto eta2 = atan(Py()/Px());
    auto eta2 = acos(Px()/eta1);
    auto eta3 = Pz();
    auto Jacmat = mat<3,3> ( cst(1.),cst(0.),cst(0.),
			     cst(0.),cst(1.),cst(0.),
			     cst(0.),cst(0.),cst(1.)  );
#endif

    auto Jac = det(Jacmat);
    auto Jacmatinv = inv(Jacmat);
    auto Jacmattrans = trans(Jacmat);
    auto Jacmattransinv = inv(Jacmattrans);


    auto lambda = doption( "parameters.l" );

    //auto uexact = pow(eta1,4./3.)*sin(4*eta2/3);
    auto uexact = pow(eta1,doption("parameters.p"))*sin(4*eta2/3);
    std::cout << "int(u) = " << integrate(elements(mesh), uexact*Jac).evaluate() << std::endl;


    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
		   _expr=trans(Jacmattransinv*trans(gradt(u)))*(Jacmattransinv*trans(grad(v)))*Jac );

    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u,
	     _expr=uexact );

    a.solve( _rhs=l, _solution=u );

    std::cout << "Error L2 norm : " << normL2(_range=elements(mesh), _expr=uexact-idv(u)) << std::endl;

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();


    // auto f = expr<3,1>(soption("functions.f"));
    // f.setParameterValues( {
    // 	    { "lambda", doption( _name="parameters.l" ) } } );
    // std::cout << f << std::endl;


    // auto eta1 = sqrt(Px()*Px()+Py()+Py()+Pz()+Pz());
    // auto eta2 = acos(Pz()/sqrt(Px()*Px()+Py()+Py()+Pz()+Pz()));
    // auto eta3 = atan(Py()/Px());

    // auto l = 4.493409458;

    // auto v = vec(2*l/(eta1*eta1)*(l/eta1*sin(eta1/l)-cos(eta1/l))*cos(eta2),
    // 		      -1/eta1*(l/eta1*cos(eta1/l)-l*l/(eta1*eta1)*sin(eta1/l)+sin(eta1/l))*sin(eta2),
    // 		      1/eta1*(l/eta1*sin(eta1/l)-cos(eta1/l))*sin(eta2));

    // auto u = Vh->element(f);

    // auto e = exporter( mesh );
    // e->add("u",u);
    // e->save();
    
    return 0;
}
