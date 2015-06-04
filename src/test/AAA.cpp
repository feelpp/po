#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

typedef Mesh<Simplex<3> > meshdim;
typedef bases<Lagrange<1,Scalar> > baseP1;
typedef bases<Lagrange<1,Scalar>, Lagrange<0,Scalar> > baseP1P0;
typedef FunctionSpace<meshdim, baseP1> FEspaceP1;
typedef FunctionSpace<meshdim, baseP1P0> FEspaceP1P0;



int main(int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="AAA",
                                   _author="BS-RH",
                                   _email="benjamin.surowiec@plasticomnium.com"));

    auto mesh = loadMesh(new meshdim);







}
