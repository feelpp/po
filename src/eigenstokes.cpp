/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-04-17

  Copyright (C) 2014 Feel++ Consortium

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
#include <feel/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feeldiscr/pchv.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	po::options_description stokesoptions( "Stokes options" );
	stokesoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc=stokesoptions,
                     _about=about(_name="ge_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

    //auto Vh = THch<1>( mesh );
    //auto Vh = FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<2,Scalar>, Lagrange<2,Scalar>, Lagrange<2,Scalar>, Lagrange<1,Scalar> > >::New( mesh );
    auto Vh = FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<2,Scalar>, Lagrange<1,Scalar> > >::New(mesh);
    //auto Vh = FunctionSpace<Mesh<Simplex<3> >, bases<Lagrange<2,Scalar> > >::New(mesh);
    auto U = Vh->element();
    auto V = Vh->element();

    // auto u1 = U.template element<0>();
    // auto u2 = U.template element<1>();
    // auto u3 = U.template element<2>();
    // auto v1 = V.template element<0>();
    // auto v2 = V.template element<1>();
    // auto v3 = V.template element<2>();

    // auto p = U.template element<3>();
    // auto q = V.template element<3>();


    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    // auto deft = gradt( u );
    // auto def = grad( v );

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh);
    // a = integrate( _range=elements( mesh ), _expr=inner( deft,def ) );
    // a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );
    // a += integrate( _range=elements( mesh ), _expr=1e-6*idt(p)*id(q) );
    // a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=zero<3>() ) ;

    // a = integrate( elements( mesh ),
    //                _expr=dxt(u1)*dx(v1) + dyt(u1)*dy(v1) + dzt(u1)*dz(v1)
    //                + dxt(u2)*dx(v2) + dyt(u2)*dy(v2) + dzt(u2)*dz(v2)
    //                + dxt(u3)*dx(v3) + dyt(u3)*dy(v3) + dzt(u3)*dz(v3)
    //                + (dxt(u1)+dyt(u2)+dzt(u3))*id(q)
    //                + (dx(v1)+dy(v2)+dz(v3))*idt(p)
    //                + 1e-20*idt(p)*id(q) );

    // a += on(_range=markedfaces(mesh, 1), _rhs=l, _element=u3, _expr=cst(0.));
    // a += on(_range=markedfaces(mesh, 2), _rhs=l, _element=u3, _expr=cst(0.));
    // a += on(_range=markedfaces(mesh, 3), _rhs=l, _element=u1, _expr=cst(0.));
    // a += on(_range=markedfaces(mesh, 3), _rhs=l, _element=u2, _expr=cst(0.));

    // auto b = form2( _trial=Vh, _test=Vh);
    // b = integrate( _range=elements( mesh ),
    //                _expr=idt(u1)*id(v1)
    //                +idt(u2)*id(v2)
    //                +idt(u3)*id(v3)
    //                +0.*idt(p)*id(p) );

    a = integrate(_range=elements(mesh),
                  _expr= dxt(u)*dx(v) + dyt(u)*dy(v) + dzt(u)*dz(v) - 0.*idt(u)*id(v)
                  + 1e-20*idt(p)*id(q) );
    a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.) );
    a += on(_range=markedfaces(mesh,2), _rhs=l, _element=p, _expr=cst(0.) );

    auto b = form2( _trial=Vh, _test=Vh);
    b = integrate( _range=elements( mesh ),
                   _expr=idt(u)*id(v)
                   + 0.*idt(p)*id(q) );

    auto modes= veigs( _formA=a, _formB=b );

    auto Xh = Pchv<2>(mesh);
    auto g = Xh->element();
    auto modeTmp = Vh->element();
    auto e = exporter( _mesh=mesh );
    int i = 0;
    for( auto const& mode: modes )
    {
        // LOG_IF(WARNING, (i==0)&&( math::abs( mode.first - 52.34 ) > 1e-1 ) )
        //     << "Invalid stokes first eigen value " << mode.first << " should be " << 52.34;
        // modeTmp = mode.second;
        // g = vf::project(_space=Vh, _range=elements(mesh),
        //                 _expr=vec(idv(modeTmp.template element<0>()),
        //                           idv(modeTmp.template element<1>()),
        //                           idv(modeTmp.template element<2>()) ) );
        p = mode.second.element<1>();
        g = vf::project(_space=Xh, _range=elements(mesh),
                        _expr=vec(cst(0.),
                                  cst(0.),
                                  idv(mode.second.element<0>()) ) );
        // u = mode.second.element<0>();
        // e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), u );
        e->add( ( boost::format( "mode-p-%1%" ) % i ).str(), p );
        e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), g );
        auto normL2Div = normL2( _range=elements(mesh), _expr=divv(g) );
        auto normL2Bord = normL2( _range=elements(mesh), _expr=trans(curlv(g))*N() );
        if ( Environment::isMasterRank() )
        {
            std::cout << "Lambda_" << i << " = " <<  mode.first << " Divergence = " << normL2Div << " Curl(g).n = " << normL2Bord << "\n";
        }
        i++;
    }
    e->save();

    return 0;
}
