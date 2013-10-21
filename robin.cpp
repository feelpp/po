/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-14

  Copyright (C) 2013 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file robin.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-14
 */
// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
	using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="po_laplacian_robin",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org"));
    auto mesh = unitSquare();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+= integrate(_range=boundaryfaces(mesh),
                  _expr=0.01*idt(u)*id(v) );
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));
    a.solve(_rhs=l,_solution=u);

    auto e = exporter( _mesh=mesh, _name="robin" );
    e->add( "u", u );
    e->save();
}
