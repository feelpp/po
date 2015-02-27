#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelalg/solvereigen.hpp>


using namespace Feel;
using namespace Feel::vf;

int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_test",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = unitSphere();

    auto Vh = Ned1h<0>(mesh);
    auto u = Vh->element();
    auto v = Vh->element();

    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    // a = integrate( _range=elements( mesh ), _expr=inner(gradt(u),grad(v)));
    auto matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Vh, _trial=Vh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    matB->close();

    auto modes = eigs( _matrixA=matA,
                       _matrixB=matB,
                       _solver=KRYLOVSCHUR,
                       _problem=GHEP,
                       _transform=SINVERT,
                       _spectrum=SMALLEST_MAGNITUDE
                       );

    int i = 0;
    for( auto const& mode: modes )
    {
        if ( Environment::isMasterRank() )
            std::cout << "eigenvalue " << i++ << " = (" << boost::get<0>(mode.second) << "," <<  boost::get<1>(mode.second) << ")" << std::endl;

        auto tmp1 = backend()->newVector( matA->mapRowPtr());
        auto tmp2 = backend()->newVector( matA->mapRowPtr());
        matA->multVector(boost::get<2>(mode.second), tmp1);
        matB->multVector(boost::get<2>(mode.second), tmp2);
        tmp1->add(-boost::get<0>(mode.second), tmp2);
        auto n = tmp1->linftyNorm();
        if( Environment::isMasterRank())
            std::cout << "norm = " << n << std::endl;
    }

    return 0;
}
