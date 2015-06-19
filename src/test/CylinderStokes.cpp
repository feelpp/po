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
typedef Lagrange<2,Scalar> BaseP2;
typedef Lagrange<1,Scalar> BaseP1;
typedef bases<BaseP2,BaseP1> BaseP2P1;
typedef FunctionSpace<MeshDim,BaseP2P1> FESpaceP2P1;

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
    auto Xh = FESpaceP2P1::New( Th );

    // UNKOWN VARIABLE
    auto PSI2 = Xh->element( );
    auto psi2 = PSI2.element<0>( );
    auto prpsi2 = PSI2.element<1>( );

    // TEST FUNCTION
    auto V = Xh->element( );
    auto v = V.element<0>( );
    auto prv = V.element<1>( );

    // RIGHT HAND SIDE : F AND ALPHA2
    auto f = expr( soption( "functions.f" ) );
    auto alpha2 = expr( soption( "functions.g" ) );

    // LINEAR FORM
    auto l2 = form1( _test=Xh );

    l2 = integrate( _range=elements( Th ),
                    _expr=trans( f )*id( v ) );

    l2 += integrate( _range=markedfaces( Th, 1 ),
                     _expr=-alpha2*id( v ) );

    l2 += integrate( _range=markedfaces( Th, 2 ),
                     _expr=alpha2*id( v ) );

    // BILINEAR FORM
    auto a2 = form2( _trial=Xh,
                     _test=Xh );

    a2 = integrate ( _range=elements( Th ),
                     _expr=inner(gradt( psi2 ),gradd( v )) );

    a2 += integrate ( _range=elements( Th ),
                      _expr=idt( psi2 )*id( prv ) );

    a2 += integrate ( _range=elements( Th ),
                      _expr=id( v )*idt( prpsi2 ) );

    // SOLVING PROBLEM
    a2.solve( _rhs=l2,
              _solution=PSI2 );

    // RESULT STORAGE
    auto e = exporter( Th );
    e->add( "psi2", psi2 );
    e->save( );

    // END OF MAIN PROGRAM
    return 0;
}
