#include <feel/feelvf/vf.hpp>

#include "darcy.h"

Darcy::Darcy(mesh_ptrtype mesh, std::string g_s):super()
{
    this->mesh = mesh;
    this->g_s = g_s;
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "-----Darcy-----" << std::endl;
        std::cout << g_s << std::endl;
    }

}

void
Darcy::run()
{
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( this->g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", option( _name="radius" ).template as<double>() },
                { "speed", option( _name="speed" ).template as<double>() } } );

    auto Vh = space_type::New( mesh );
    U = Vh->element();
    auto V = Vh->element();
    auto u = U.template element<0>() ;
    auto v = V.template element<0>() ;
    auto p = U.template element<1>() ;
    auto q = U.template element<1>() ;
    // auto lambda = U.template element<2>() ;
    // auto nu = V.template element<2>() ;

    auto l = form1( _test=Vh );
    l = integrate( _range=markedfaces(mesh, 1), // inflow
                   _expr=-g*id(q) );
    l += integrate( _range=markedfaces(mesh, 2), // outflow
                    _expr=g*id(q) );
    l += integrate( _range=markedfaces(mesh, 3), // wall
                    _expr=cst(0.)*id(q) );
    /// [rhs]

    /// [bilinear]
    auto a = form2( _test=Vh, _trial=Vh );
    a = integrate( _range=elements(mesh),
                   _expr = -inner(idt( u ),id( v )) );
    a = integrate( _range=elements(mesh),
                   _expr = inner(trans(gradt( p )),id( v )) );
    a = integrate( _range=elements(mesh),
                   _expr = inner(trans(grad( q )),idt( u )) );
    // a = integrate( _range=elements(mesh),
    //                _expr = idt( p ) * id( nu ) );
    // a = integrate( _range=elements(mesh),
    //                _expr = idt( lambda ) * id( q ) );

    a.solve( _name="psi0Div", _rhs=l, _solution=U );
    /// [bilinear]

}
