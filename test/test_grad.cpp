#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

class TestGrad
:
public Simget
{
    typedef Simget super;
    static const uint16_type Dim = 3;
    static const uint16_type Order = 3;
 public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    typedef bases<Lagrange<Order, Vectorial> > basis_vtype;
    typedef FunctionSpace<mesh_type, basis_vtype > space_vtype;
    typedef typename space_vtype::element_type element_vtype;
    typedef boost::shared_ptr<space_vtype> space_vptrtype;

    typedef bases<Lagrange<Order, Scalar> > basis_stype;
    typedef FunctionSpace<mesh_type, basis_stype > space_stype;
    typedef typename space_stype::element_type element_stype;
    typedef boost::shared_ptr<space_stype> space_sptrtype;


    TestGrad( mesh_ptrtype );
    void run();
    void prob();
    void test();
 private:
    mesh_ptrtype mesh;
    space_sptrtype Vh;
    space_vptrtype Vh2;
};

TestGrad::TestGrad( mesh_ptrtype mesh)
{
    this->mesh = mesh;
    this->Vh = space_stype::New( mesh );
    LOG(INFO) << "----- Vh -----\n";
    LOG(INFO) << "[dof] number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Vh -----\n";
        std::cout << "[dof] number of dof: " << Vh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
    }

    this->Vh2 = space_vtype::New( mesh );
    LOG(INFO) << "----- Vh2 -----\n";
    LOG(INFO) << "[dof] number of dof: " << Vh2->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Vh2->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Vh2 -----\n";
        std::cout << "[dof] number of dof: " << Vh2->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Vh2->nLocalDof() << "\n";
    }
}

void
TestGrad::run()
{
    if( boption("testG") )
        test();
    else
        prob();
}

void
TestGrad::prob()
{
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    auto a = form2( _test=Vh, _trial=Vh);
    auto b = form2( _test=Vh, _trial=Vh);

    a = integrate(elements(mesh), inner(gradt(u),grad(v)));

    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.) );

    b = integrate( elements(mesh), inner(idt(u),id(v)) );

    auto modes = veigs( _formA=a, _formB=b);

    auto cu = Vh2->element();
    auto acu = form2( _test=Vh2, _trial=Vh2 );
    auto bcu = form1( _test=Vh2 );

    acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );

    auto ccu = Vh->element();
    auto accu = form2( _test=Vh, _trial=Vh );
    auto bccu = form1( _test=Vh );

    accu  = integrate( elements(mesh), idt(ccu)*id(ccu) );

    auto hessu = Vh->element();
    auto ahess = form2( _test=Vh, _trial=Vh );
    auto bhess = form1( _test=Vh );

    ahess = integrate( elements(mesh), trans(idt(hessu))*id(hessu) );


    auto e =  exporter( _mesh=mesh );

    int i = 0;
    double moyH=0, moyL=0;

    for( auto const& mode : modes )
    {
        bcu = integrate( elements(mesh), gradv(mode.second)*id(cu) );
        acu.solve(_rhs=bcu, _solution=cu, _name="grad");

        bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
        accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

        bhess = integrate( elements(mesh), trace(hessv(mode.second))*id(hessu) );
        ahess.solve( _rhs=bhess, _solution=hessu, _name="hess" );

        e->add( ( boost::format( "mode-%1%" ) % i ).str(), mode.second );
        e->add( ( boost::format( "lapl-%1%" ) % i ).str(), ccu );
        e->add( ( boost::format( "hess-%1%" ) % i ).str(), hessu );

        auto erreur2 = normL2(elements(mesh), divv(cu) + mode.first*idv(mode.second) );
        auto erreur = normL2(elements(mesh), idv(ccu) + mode.first*idv(mode.second) );
        auto erreur3 = normL2(elements(mesh), idv(hessu) + mode.first*idv(mode.second) );
        moyL += erreur;
        moyH += erreur3;

        if ( Environment::worldComm().isMasterRank() )
            std::cout << "err(" << i << ") = " << erreur << "     err2(" << i << ") = " << erreur2 << "     err3(" << i << ") = " << erreur3 << std::endl;
        i++;
    }
    moyL = moyL/i;
    moyH = moyH/i;

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "moyL = " << moyL << "\tmoyH = " << moyH << std::endl;

    e->save();
}

void
TestGrad::test()
{
    auto ext1 = expr(soption("t1"));
    auto laplacian_t1 = laplacian(ext1);
    auto grad_t1 = grad<3>(ext1);
    auto exgt1 = expr<3,1>(soption("gt1"));
    auto exlt1 = expr(soption("lt1"));
    auto vt1 = Vh->element(ext1);
    auto gt1 = Vh2->element(exgt1);
    auto lt1 = Vh->element(exlt1);

    auto ext2 = expr(soption("t2"));
    auto laplacian_t2 = laplacian(ext2);
    auto grad_t2 = grad<3>(ext2);
    auto exgt2 = expr<3,1>(soption("gt2"));
    auto exlt2 = expr(soption("lt2"));
    auto vt2 = Vh->element(ext2);
    auto gt2 = Vh2->element(exgt2);
    auto lt2 = Vh->element(exlt2);

    auto ext3 = expr(soption("t3"));
    auto laplacian_t3 = laplacian(ext3);
    auto grad_t3 = grad<3>(ext3);
    auto exgt3 = expr<3,1>(soption("gt3"));
    auto exlt3 = expr(soption("lt3"));
    auto vt3 = Vh->element(ext3);
    auto gt3 = Vh2->element(exgt3);
    auto lt3 = Vh->element(exlt3);

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "t1 = " << ext1 << "\tgrad(t1) = " << grad_t1 << "\tlap(t1) = " << laplacian_t1 << std::endl;
        std::cout << "t2 = " << ext2 << "\tgrad(t2) = " << grad_t2 << "\tlap(t2) = " << laplacian_t2 << std::endl;
        std::cout << "t3 = " << ext3 << "\tgrad(t3) = " << grad_t3 << "\tlap(t3) = " << laplacian_t3 << std::endl;
    }

    // [sysGrad]
    auto cu = Vh2->element();
    auto acu = form2( _test=Vh2, _trial=Vh2 );
    auto bcu = form1( _test=Vh2 );

    acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );
    // [sysGrad]
    // [sysDiv]
    auto ccu = Vh->element();
    auto accu = form2( _test=Vh, _trial=Vh );
    auto bccu = form1( _test=Vh );

    accu  = integrate( elements(mesh), idt(ccu)*id(ccu) );
    // [sysDiv]
    // [sysHess]
    auto hessu = Vh->element();
    auto ahess = form2( _test=Vh, _trial=Vh );
    auto bhess = form1( _test=Vh );

    ahess = integrate( elements(mesh), trans(idt(hessu))*id(hessu) );
    // [sysHess]

    double err1, err2, err3, moyL=0, moyH=0;

    // [rhsGrad]
    bcu = integrate( elements(mesh), gradv(vt1)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");
    // [rhsGrad]
    // [rhsDiv]
    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );
    // [rhsDiv]
    // [rhsHess]
    bhess = integrate( elements(mesh), trace(hessv(vt1))*id(hessu) );
    ahess.solve( _rhs=bhess, _solution=hessu, _name="hess" );
    // [rhsHess]

    err1 = normL2(elements(mesh), idv(cu) - idv(gt1) );
    err2 = normL2(elements(mesh), idv(ccu) - idv(lt1) );
    err3 = normL2(elements(mesh), idv(hessu) - idv(lt1) );
    moyL+=err2;
    moyH+=err3;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "errGrad1 = " << err1 << "     errLap1 = " << err2 << "     errHess1 = " << err3 << std::endl;


    bcu = integrate( elements(mesh), gradv(vt2)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");

    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

    bhess = integrate( elements(mesh), trace(hessv(vt2))*id(hessu) );
    ahess.solve( _rhs=bhess, _solution=hessu, _name="hess" );

    err1 = normL2(elements(mesh), idv(cu) - idv(gt2) );
    err2 = normL2(elements(mesh), idv(ccu) - idv(lt2) );
    err3 = normL2(elements(mesh), idv(hessu) - idv(lt2) );
    moyL+=err2;
    moyH+=err3;

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "errGrad2 = " << err1 << "     errLap2 = " << err2 << "     errHess2 = " << err3 << std::endl;


    bcu = integrate( elements(mesh), gradv(vt3)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");

    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

    bhess = integrate( elements(mesh), trace(hessv(vt3))*id(hessu) );
    ahess.solve( _rhs=bhess, _solution=hessu, _name="hess" );

    err1 = normL2(elements(mesh), idv(cu) - idv(gt3) );
    err2 = normL2(elements(mesh), idv(ccu) - idv(lt3) );
    err3 = normL2(elements(mesh), idv(hessu) - idv(lt3) );
    moyL+=err2;
    moyH+=err3;

    moyL=moyL/3;
    moyH=moyH/3;

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "errGrad3 = " << err1 << "     errLap3 = " << err2 << "     errHess3 = " << err3 << std::endl;
        std::cout << "moyL = " << moyL << "\tmoyH = " << moyH <<std::endl;
    }



}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "testG", po::value<bool>()->default_value( true ), "use test" )
        ( "t1", po::value<std::string>()->default_value( "sin(z)+cos(y):y:z" ), "test function 1" )
        ( "t2", po::value<std::string>()->default_value( "sin(z)+cos(x):x:z" ), "test function 2" )
        ( "t3", po::value<std::string>()->default_value( "sin(x*y):x:y" ), "test function 3" )
        ( "gt1", po::value<std::string>()->default_value( "{0,-sin(y),cos(z)}:y:z" ), "grad 1" )
        ( "gt2", po::value<std::string>()->default_value( "{-sin(x),0,cos(z)}:x:z" ), "grad 2" )
        ( "gt3", po::value<std::string>()->default_value( "{y*cos(x*y),x*cos(x*y),0}:x:y" ), "grad 3" )
        ( "lt1", po::value<std::string>()->default_value( "-cos(y)-sin(z):y:z" ), "laplacian 1" )
        ( "lt2", po::value<std::string>()->default_value( "-cos(x)-sin(z):x:z" ), "laplacian 2" )
        ( "lt3", po::value<std::string>()->default_value( "-(x*x+y*y)*sin(x*y):x:y" ), "laplacian 3" );
    return myappOptions;
}


int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=feel_options().add( backend_options( "grad" ) ).add ( backend_options( "div" ) ).add( backend_options( "hess" ) ),
                     _about=about(_name="po_test_grad",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com") );

    Application app;

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    app.add( new TestGrad(mesh) );
    app.run();

    return 0;
}
