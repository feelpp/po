#include <feel/feelvf/vf.hpp>

#include "poisson.h"

Poisson::Poisson(mesh_ptrtype mesh, std::string g_s):super()
{
    this->mesh = mesh;
    this->g_s = g_s;
    Xh = space_type::New( mesh );
    gradu = Xh->element();
}

void
Poisson::run()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Poisson -----" << std::endl;
        std::cout << g_s << std::endl;
    }

    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( this->g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", option( _name="radius" ).template as<double>() },
                { "speed", option( _name="speed" ).template as<double>() } } );

    auto Vh = mlSpace_type::New( mesh );
    U = Vh->element();
    auto V = Vh->element();
    auto u = U.template element<0>() ;
    auto lambda = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto nu = V.template element<1>() ;

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=inner(gradt(u),grad(v)) + id( v )*idt( lambda ) + idt( u )*id( nu ) );

    auto l = form1( _test=Vh );
    l = integrate( _range=markedfaces(mesh, 1), // inflow
                   _expr=-g*id(v) );
    l += integrate( _range=markedfaces(mesh, 2), // outflow
                    _expr=g*id(v) );
    l += integrate( _range=markedfaces(mesh, 3), // wall
                    _expr=cst(0.)*id(v) );

    a.solve( _name="psi0", _rhs=l, _solution=U );

    auto b = form2( _trial=Xh, _test=Xh );
    b = integrate( _range=elements(mesh),
                   _expr=inner(id(gradu),idt(gradu)) );
    auto k = form1( _test=Xh );
    k = integrate( _range=elements(mesh),
                   _expr=inner(trans(gradv( U.template element<0>() )),id(gradu)) );
    // gradu is the L2 projection of grad(psi0) over Xh
    b.solve( _name="gradpsi0", _rhs=k, _solution=gradu );
}
