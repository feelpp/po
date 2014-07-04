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
    static const uint16_type Order = 2;
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
    this->Vh2 = space_vtype::New( mesh );
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

    LOG(INFO) << "----- Vh -----\n";
    LOG(INFO) << "[dof] number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Vh -----\n";
        std::cout << "[dof] number of dof: " << Vh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
    }


    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    auto a = form2( _test=Vh, _trial=Vh);
    auto b = form2( _test=Vh, _trial=Vh);

    a = integrate(elements(mesh), inner(gradt(u),grad(v)));

    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.) );

    b = integrate( elements(mesh), inner(idt(u),id(v)) );
    //b += integrate( boundaryfaces(mesh), 2*id(v) );

    auto modes = veigs( _formA=a, _formB=b);

    auto cu = Vh2->element();
    auto acu = form2( _test=Vh2, _trial=Vh2 );
    auto bcu = form1( _test=Vh2 );

    acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );

    auto ccu = Vh->element();
    auto accu = form2( _test=Vh, _trial=Vh );
    auto bccu = form1( _test=Vh );

    accu  = integrate( elements(mesh), idt(ccu)*id(ccu) );

    auto e =  exporter( _mesh=mesh );

    int i = 0;

    for( auto const& mode : modes )
    {
        bcu = integrate( elements(mesh), gradv(mode.second)*id(cu) );
        //acu += on(boundaryfaces(mesh), _rhs=bcu, _element=cu, _expr=zero<3>());
        acu.solve(_rhs=bcu, _solution=cu, _name="grad");

        bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
        //accu += on(boundaryfaces(mesh), _rhs=bccu, _element=ccu, _expr=cst(0.) );
        accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

        e->add( ( boost::format( "mode-%1%" ) % i ).str(), mode.second );
        e->add( ( boost::format( "lapl-%1%" ) % i ).str(), ccu );

        auto erreur = normL2(elements(mesh), idv(ccu)+mode.first*idv(mode.second) );
        auto erreur2 = normL2(elements(mesh), divv(cu)+mode.first*idv(mode.second) );
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "err(" << i << ") = " << erreur << "     err2(" << i << ") = " << erreur2 << std::endl;
        i++;
    }

    e->save();
}

void
TestGrad::test()
{
    auto t1 = expr(soption("t1"));
    auto laplacian_t1=laplacian(t1);
    auto exgt1 = expr<3,1>(soption("gt1"));
    auto exlt1 = expr(soption("lt1"));
    auto vt1 = Vh->element(t1);
    auto gt1 = Vh2->element(exgt1);
    auto lt1 = Vh->element(exlt1);

    auto ext2 = expr(soption("t2"));
    auto laplacian_t2=laplacian(ext2);
    auto exgt2 = expr<3,1>(soption("gt2"));
    auto exlt2 = expr(soption("lt2"));
    auto vt2 = Vh->element(ext2);
    auto gt2 = Vh2->element(exgt2);
    auto lt2 = Vh->element(exlt2);

    auto ext3 = expr(soption("t3"));
    auto laplacian_t3=laplacian(ext3);
    auto exgt3 = expr<3,1>(soption("gt3"));
    auto exlt3 = expr(soption("lt3"));
    auto vt3 = Vh->element(ext3);
    auto gt3 = Vh2->element(exgt3);
    auto lt3 = Vh->element(exlt3);

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "t1 = " << t1 << "\tlap(t1) = " << laplacian_t1 << std::endl;
        std::cout << "t2 = " << ext2 << "\tlap(t2) = " << laplacian_t2 << std::endl;
        std::cout << "t3 = " << ext3 << "\tlap(t3) = " << laplacian_t3 << std::endl;
    }


    auto cu = Vh2->element();
    auto acu = form2( _test=Vh2, _trial=Vh2 );
    auto bcu = form1( _test=Vh2 );

    acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );

    auto ccu = Vh->element();
    auto accu = form2( _test=Vh, _trial=Vh );
    auto bccu = form1( _test=Vh );

    accu  = integrate( elements(mesh), idt(ccu)*id(ccu) );

    double err1, err2;


    bcu = integrate( elements(mesh), gradv(vt1)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");

    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );


    err1 = normL2(elements(mesh), idv(cu)-idv(gt1) );
    err2 = normL2(elements(mesh), idv(ccu)-idv(lt1) );

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "errGrad1 = " << err1 << "     errLap1 = " << err2 << std::endl;


    bcu = integrate( elements(mesh), gradv(vt2)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");

    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

    err1 = normL2(elements(mesh), idv(cu)-idv(gt2) );
    err2 = normL2(elements(mesh), idv(ccu)-idv(lt2) );

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "errGrad2 = " << err1 << "     errLap2 = " << err2 << std::endl;


    bcu = integrate( elements(mesh), gradv(vt3)*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="grad");

    bccu = integrate( elements(mesh), divv(cu)*id(ccu) );
    accu.solve( _rhs=bccu, _solution=ccu, _name="div" );

    err1 = normL2(elements(mesh), idv(cu)-idv(gt3) );
    err2 = normL2(elements(mesh), idv(ccu)-idv(lt3) );

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "errGrad3 = " << err1 << "     errLap3 = " << err2 << std::endl;



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
                     _desc_lib=feel_options().add( backend_options( "grad" ) ).add ( backend_options( "div" ) ),
                     _about=about(_name="po_testgrad",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com") );

    Application app;

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    app.add( new TestGrad(mesh) );
    app.run();

    return 0;
}
