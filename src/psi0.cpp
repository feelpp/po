/** \file psi0.cpp
    \brief Source file for the class Psi0
*/

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/ginac.hpp>

#include "psi0.hpp"


Psi0::Psi0(mesh_ptrtype mesh, std::string g_s):super()
{
    this->mesh = mesh;
    this->g_s = g_s;

    Xh = space_type::New( mesh );
    gradu = Xh->element();
}


void
Psi0::run()
{
    if ( Environment::worldComm().isMasterRank() ){
        std::cout << "----- Psi0 -----" << std::endl;
        std::cout << g_s << std::endl;
    }

    if( boption("computeP0") )
        compute_psi0();
    else
        load_psi0();
}

void
Psi0::compute_psi0() {
    // [option]
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( this->g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", doption( _name="radius" ) },
                { "speed", doption( _name="speed" ) } } );
    // [option]


    auto Vh = mlSpace_type::New( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.template element<0>() ;
    auto lambda = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto nu = V.template element<1>() ;

    LOG(INFO) << "----- PSI0 -----\n";
    LOG(INFO) << "[dof] number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(U): " << Vh->template functionSpace<0>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Vh->template functionSpace<0>()->nLocalDof() << "\n";
    LOG(INFO) << "[dof] number of dof(P): " << Vh->template functionSpace<1>()->nDof() << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Vh->template functionSpace<1>()->nLocalDof() << "\n";

    if ( Environment::isMasterRank() ){
        std::cout << "----- PSI0 -----\n";
        std::cout << "[dof] number of dof: " << Vh->nDof() << "\n";
        std::cout << "[dof] number of dof/proc: " << Vh->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(U): " << Vh->template functionSpace<0>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(U): " << Vh->template functionSpace<0>()->nLocalDof() << "\n";
        std::cout << "[dof] number of dof(P): " << Vh->template functionSpace<1>()->nDof() << "\n";
        std::cout << "[dof] number of dof/proc(P): " << Vh->template functionSpace<1>()->nLocalDof() << "\n";
    }


    // [bilinearA]
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=inner(gradt(u),grad(v))
                   // [bilinearA]
                   // [bilinearB]
                   + id( v )*idt( lambda ) + idt( u )*id( nu ) );
    // [bilinearB]

    // [rhs]
    auto l = form1( _test=Vh );
    l = integrate( _range=markedfaces(mesh, 1), // inflow
                   _expr=-g*id(v) );
    l += integrate( _range=markedfaces(mesh, 2), // outflow
                    _expr=g*id(v) );
    l += integrate( _range=markedfaces(mesh, 3), // wall
                    _expr=cst(0.)*id(v) );
    // [rhs]

    a.solve( _name="psi0", _rhs=l, _solution=U );


    // [gradpsi0]
    auto b = form2( _trial=Xh, _test=Xh );
    b = integrate( _range=elements(mesh),
                   _expr=inner(id(gradu),idt(gradu)) );
    auto k = form1( _test=Xh );
    k = integrate( _range=elements(mesh),
                   _expr=inner(trans(gradv( U.template element<0>() )),id(gradu)) );
    // [gradpsi0]

    // gradu is the L2 projection of grad(psi0) over Xh
    b.solve( _name="gradpsi0", _rhs=k, _solution=gradu );

    std::string path = "psi0";
    gradu.save(_path=path);

}

void
Psi0::load_psi0() {
    std::string path = "psi0";
    gradu.load(_path=path);
}
