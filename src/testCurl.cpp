#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/mpi/timer.hpp>

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/ginac.hpp>

#include "testCurl.hpp"

TestCurl::TestCurl( mesh_ptrtype mesh )
{
    this->mesh = mesh;
    this->Vh = space_vtype::New( mesh );
    this->test = std::vector<element_vtype>(15, Vh->element() );
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

    auto t = expr<3,1>( soption("t") );
    auto t1 = expr<3,1>( soption("t1") );
    auto t2 = expr<3,1>( soption("t2") );

    //    auto curl_t = curl(t);             //ne fonctionne pas
    //auto curl2_t = curl(curl_t);

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "t : " << t << "\nt1 : " << t1 << "\nt2 : " << t2/* << "\ncurl(t) : " << curl_t << "\ncurl2(t) : " << curl2_t*/ << std::endl;
    }

    // solution exacte
    auto vt = Vh->element(t);     // t
    auto vt1 = Vh->element(t1);   // curl(t)
    auto vt2 = Vh->element(t2);   // curl2(t)


    // utilisation de project et curlv
    auto curl_vt = vf::project(_space=Vh, _range=elements(mesh),
                               _expr=curlv(vt) );
    auto curl2_vt = vf::project(_space=Vh, _range=elements(mesh),
                                _expr=curlv(curl_vt) );


    // resolution d'un systeme
    auto cu = Vh->element();
    auto ccu = Vh->element();
    auto acu = form2( _test=Vh, _trial=Vh );
    auto bcu = form1( _test=Vh );
    acu = integrate( elements(mesh), trans(idt(cu))*id(cu) );

    // premier systeme pour curl(t)
    bcu = integrate(elements(mesh),
                    trans(curlv(vt))*id(cu) );
    acu.solve(_rhs=bcu, _solution=cu, _name="curl" );

    // second systeme pour curl2(t)
    bcu = integrate(elements(mesh),
                    trans(curlv(cu))*id(cu) );
    acu.solve(_rhs=bcu, _solution=ccu, _name="curl2" );

    // norme avec project
    auto norm2curl_1 = normL2( elements(mesh), curlv(vt) - idv(vt1) );
    auto norm2curl_2 = normL2( elements(mesh), idv(curl_vt) - idv(vt1) );
    auto norm2curl_3 = normL2( elements(mesh), idv(cu) - idv(vt1) );

    // norme avec systeme
    auto norm2curl2_1 = normL2( elements(mesh), curlv(curl_vt) - idv(vt2) );
    auto norm2curl2_2 = normL2( elements(mesh), idv(curl2_vt) - idv(vt2) );
    auto norm2curl2_3 = normL2( elements(mesh), idv(ccu) - idv(vt2) );

    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "\t\t\tcurlv\t\t project \t\t systeme" << std::endl;
        std::cout << "norm(curlv(t)-idv(vt1))=\t" << norm2curl_1 << "\t" << norm2curl_2 << "\t" << norm2curl_3 << std::endl;
        std::cout << "norm(curl2v(t)-idv(vt1))=\t" << norm2curl2_1 << "\t" << norm2curl2_2 << "\t" << norm2curl2_3 << std::endl;
    }

    test[0] = vt;
    test[1] = vt1;
    test[2] = vt2;
    test[3] = curl_vt;
    test[4] = curl2_vt;
    test[5] = cu;
    test[6] = ccu;

}
