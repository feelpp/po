#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsphere.hpp>
//#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/ginac.hpp>
//#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>

using namespace Feel;

int main(int argc, char** argv)
{
	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="po_testD1",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = unitSphere();

    auto lambda = doption("lambda");
    auto V = expr<3,1>(soption("f"));

    return 0;
}
