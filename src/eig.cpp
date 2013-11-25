// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_eig",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #

    //# marker2 #
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    //# endmarker2 #

    //# marker3 #

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _expr );
    
    auto b = form2(_trial=Vh, _test=Vh);
    b = integrate(_range=elements(mesh), _expr=idt(u)*id(v));
    //# endmarker3 #

    //# marker4 #
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    return 0;
    //# endmarker4 #
}
