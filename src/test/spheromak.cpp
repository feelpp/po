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

    auto e = exporter( mesh );
    e->add( "u", u );
    e->save();

}
