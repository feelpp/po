#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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
#if PO_SPH
                                  "po_test_sph_sph",
#else
                                  "po_test_sph_car",
#endif
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3> > );

    auto Vh = Pch<1>( mesh );
    std::cout << "dim Vh : " << Vh->nDof() << std::endl;

    auto u = Vh->element();
    auto v = Vh->element();

#if PO_SPH
    auto eta1 = Px();
    auto eta2 = Py();
    auto eta3 = Pz();
    auto Jacmat = mat<3,3> ( sin(eta2)*cos(eta3), eta1*cos(eta2)*cos(eta3), -eta1*sin(eta2)*sin(eta3),
                             sin(eta2)*sin(eta3), eta1*cos(eta2)*sin(eta3), eta1*sin(eta2)*cos(eta3),
                             cos(eta2),           -eta1*sin(eta2),          cst(0.) );
#else
    auto eta1 = sqrt(Px()*Px()+Py()*Py()+Pz()*Pz());
    auto eta2 = acos(Pz()/eta1);
    auto eta3 = atan(Py()/Px());
    auto Jacmat = mat<3,3> ( cst(1.), cst(0.), cst(0.),
                             cst(0.), cst(1.), cst(0.),
                             cst(0.), cst(0.), cst(1.)  );
#endif

    auto Jac = det(Jacmat);
    auto Jacmatinv = inv(Jacmat);
    auto Jacmattrans = trans(Jacmat);
    auto Jacmattransinv = inv(Jacmattrans);


    auto uexact = 1./6*pow(eta1,doption("parameters.p"))*sin(eta2);
    std::cout << "int(u) = " << integrate(elements(mesh), uexact*Jac).evaluate() << std::endl;


    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=gradt(u)*Jacmatinv*Jacmattransinv*trans(grad(v))*Jac );

    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u,
             _expr=uexact );

    a.solve( _rhs=l, _solution=u );

    std::cout << "Error L2 norm : " << normL2(_range=elements(mesh), _expr=uexact-idv(u)) << std::endl;

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();

    return 0;
}
