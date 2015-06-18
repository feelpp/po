// MANDATORY INCLUDE AND USED SPACE

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

// TYPE DEFINITON

#define FEELPP_DIM 3

typedef Mesh<Simplex<FEELPP_DIM>> MeshDim;
typedef Lagrange<1,Scalar> BaseP1;
typedef RaviartThomas<0> BaseRT0;
typedef bases<BaseRT0,BaseP1> BaseRT0P1;
typedef FunctionSpace<MeshDim,BaseRT0P1> FESpaceRT0P1;

// EXTENDED OPTION

inline po::options_description
MakeOptions( )
{
    po::options_description EXPRoptions( "CylinderStokes options" );
        EXPRoptions.add_options( )
        ( "MeanVel", po::value<double>( )->default_value( 1. ),
          "Mean Velovity of the Fluid Flow" )
        ( "Radius", po::value<double>( )->default_value( 0.5 ),
          "Radius of the Domain Cylinder" )
        ;
    return EXPRoptions;
}

// MAIN PROGRAM

int main(int argc, char** argv)
{
    // ENVIRONMENT
    Environment env( _argc=argc,
                     _argv=argv,
                     _desc=MakeOptions( ),
                     _about=about( _name="CylinderStokes",
                                   _author="B.Surowiec",
                                   _email="benjamin.surowiec@plasticomnium.com" ) );

    // MESH LOADING
    auto Th = loadMesh( new MeshDim );

    // FINITE ELEMENT SPACE
    auto Xh = FESpaceRT0P1::New( Th );

    // UNKOWN VARIABLE
    auto Q2 = Xh->element( );
    auto q2 = Q2.element<0>( );
    auto psi2 = Q2.element<1>( );

    // TEST FUNCTION
    auto V = Xh->element( );
    auto v = V.element<0>( );
    auto phiv = V.element<1>( );

    // BOUNDARY CONDITION : RIGHT HAND SIDE : F, ALPHA0 AND ALPHA2
    auto f = expr( soption( "functions.f" ) );
    auto alpha2 = expr( soption( "functions.g" ) );
    //auto alpha0 = expr( soption( "functions.h" ) );

    // LINEAR FORM
    auto l2 = from1( _test=Xh );

    l2 = integrate( _range=elements( Th ),
                    _expr=trans( f )*id( v ) );

    // BILINEAR FORM -- SEEKING SOLUTION IN HDIV0
    auto adarcy = form2( _trial=Xh,
                         _test=Xh );

    adarcy = integrate ( _range=elements( Th ),
                         _expr=-idt( q2 )*id( v ) );

    adarcy += integrate ( _range=elements( Th ),
                          _expr=-idt( psi2 )*div( v ) );

    adarcy += integrate ( _range=elements( Th ),
                          _expr=divt( q2 )*id( phiv ) );

    adarcy += on ( _range=markedfaces( Th, 1 ), rhs=l2,
                   _element=q2, _expr=-alpha2 );

    adarcy += on ( _range=markedfaces( Th, 2 ), rhs=l2,
                   _element=q2, _expr=alpha2 );

    adarcy += on ( _range=markedfaces( Th, 3 ), rhs=l2,
                   _element=q2, _expr=cst(0.) );

    // SOLVING PROBLEM
    adarcy.solve( _rhs=l2,
                  _solution=Q2 );

    // RESULT STORAGE
    auto export = exporter( Th );
    export->add( "q2", q2 );
    export->add( "psi2", psi2 );
    export->save( );

    // END OF MAIN PROGRAM
    return 0;
}
