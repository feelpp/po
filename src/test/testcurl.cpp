#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

class TestCurl:public Simget
{
    static const uint16_type Dim = 3;
    static const uint16_type Order = 3;

public:
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //typedef bases<Nedelec<0,NedelecKind::NED1> > basis_nedtype;
    typedef bases<Lagrange<Order, Vectorial> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    std::vector<element_type> test;

    TestCurl( mesh_ptrtype );
    void run();
 private:

    mesh_ptrtype mesh;
    space_ptrtype Vh;
};


TestCurl::TestCurl( mesh_ptrtype mesh )
{
    this->mesh = mesh;
    this->Vh = space_type::New( mesh );
    this->test = std::vector<element_type>(15, Vh->element() );    
}

void
TestCurl::run()
{
    LOG(INFO) << "----- Vh -----\n";
    LOG(INFO) << "[dof] number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- Vh -----\n";
        std::cout << "[dof] number of dof: " << Vh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
    }

    // [exacte]
    auto t = expr<3,1>( soption("t") );
    auto t1 = expr<3,1>( soption("t1") );
    auto t2 = expr<3,1>( soption("t2") );

    // solution exacte
    auto vt = Vh->element(t);     // t
    // [exacte]

    auto vt1 = Vh->element(t1);   // curl(t)
    auto vt2 = Vh->element(t2);   // curl2(t)

    auto curl_t = curl(t);             //ne fonctionne pas
    auto curl2_t = curl(curl_t);

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "t : " << t << "\nt1 : " << t1 << "\nt2 : " << t2/* << "\ncurl(t) : " << curl_t << "\ncurl2(t) : " << curl2_t*/ << std::endl;
    }


    // [Nedelec]
    // auto vt = Vh->element();

    // auto ap = form2( _test=Vh, _trial=Vh);
    // auto bp = form1( _test=Vh );
    // ap = integrate( elements(mesh), trans(idt(vt))*id(vt) );
    // bp = integrate( elements(mesh), trans(t)*id(vt) );
    // ap.solve(_rhs=bp, _solution=vt, _name="gi0" );
    // [Nedelec]

    // utilisation de project et curlv
    auto curl_vt = vf::project(_space=Vh, _range=elements(mesh),
                               _expr=curlv(vt) );
    auto curl2_vt = vf::project(_space=Vh, _range=elements(mesh),
                                _expr=curlv(curl_vt) );

    // [masse]
    auto cu = Vh->element();  //curl
    auto ccu = Vh->element(); //curl2
    auto a = form2( _test=Vh, _trial=Vh );
    auto b = form1( _test=Vh );
    a = integrate( elements(mesh), trans(idt(cu))*id(cu) );
    // [masse]

    // [curl]
    b = integrate(elements(mesh),
                  trans(curlv(vt))*id(cu) );
    a.solve(_rhs=b, _solution=cu, _name="curl" );
    // [curl]

    // [curl2]
    b = integrate(elements(mesh),
                  trans(curlv(cu))*id(cu) );
    a.solve(_rhs=b, _solution=ccu, _name="curl2" );
    // [curl2]

    auto norm2id = normL2( elements(mesh), idv(vt) - t );

    auto norm2curl_1 = normL2( elements(mesh), curlv(vt) - t1 );
    auto norm2curl2_1 = normL2( elements(mesh), curlv(cu) - t2 );

    auto norm2curl_2 = normL2( elements(mesh), idv(curl_vt) - t1 );
    auto norm2curl2_2 = normL2( elements(mesh), idv(curl2_vt) - t2 );

    // [erreur]
    auto norm2curl_3 = normL2( elements(mesh), idv(cu) - t1 );
    auto norm2curl2_3 = normL2( elements(mesh), idv(ccu) - t2 );
    // [erreur]

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "\t\t\t\tcurlv \t project \t systeme" << std::endl;
        std::cout << "norm(vt - t) =\t\t" << norm2id << std::endl;
        std::cout << "norm(curlv(t)-t1)=\t" << norm2curl_1 << "\t" << norm2curl_2 << "\t" << norm2curl_3 << std::endl;
        std::cout << "norm(curl2v(t)-t2)=\t" << norm2curl2_1 << "\t" << norm2curl2_2 << "\t" << norm2curl2_3 << std::endl;
    }

    test[0] = vt;
    test[1] = vt1;
    test[2] = vt2;
    test[3] = curl_vt;
    test[4] = curl2_vt;
    test[5] = cu;
    test[6] = ccu;

}

inline
po::options_description
makeOptions()
{
    po::options_description myappOptions( "PlasticOmnium options" );
    myappOptions.add_options()
        ( "t", po::value<std::string>()->default_value( "{x*x+2*y*y+3*z*z,2*x*x+3*y*y+z*z,3*x*x+y*y+2*z*z}:x:y:z" ), "pol order 2" )
        ( "t1", po::value<std::string>()->default_value( "{2*y-2*z,6*z-6*x,4*x-4*y}:x:y:z" ), "curl(t)" )
        ( "t2", po::value<std::string>()->default_value( "{-10,-6,-8}:x:y:z" ), "curl2(t)" );
    return myappOptions;
}


int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _desc_lib=feel_options().add( backend_options("curl") ).add ( backend_options("curl2") ),
                     _about=about(_name="po_testcurl",
                                  _author="Romain Hild",
                                  _email="romain.hild@plasticomnium.com") );

    Application app;

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>> );

    app.add( new TestCurl(mesh) );
    app.run();

    return 0;
}
