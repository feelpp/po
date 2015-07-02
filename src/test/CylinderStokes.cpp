// MANDATORY INCLUDE AND USED SPACE

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;


// TYPE DEFINITON

#define FEELPP_DIM 3

typedef Mesh<Simplex<FEELPP_DIM>> MeshDim;
typedef Lagrange<0,Scalar> LagrP0;
typedef Lagrange<1,Scalar> LagrP1;
typedef Lagrange<2,Scalar> LagrP2;
typedef Lagrange<1,Vectorial> LagrVecP1;
typedef Lagrange<2,Vectorial> LagrVecP2;
typedef bases<LagrP2,LagrP1> BaseP2P1;
typedef bases<LagrVecP1> BaseVP1;
typedef bases<LagrVecP1,LagrP0> BaseVP1P0;
typedef bases<LagrVecP2,LagrP1> BaseVP2P1;
typedef FunctionSpace<MeshDim,BaseP2P1> FESpaceP2P1;
typedef FunctionSpace<MeshDim,BaseVP1> FESpaceVP1;
typedef FunctionSpace<MeshDim,BaseVP1P0> FESpaceVP1P0;
typedef FunctionSpace<MeshDim,BaseVP2P1> FESpaceVP2P1;


// EXTENDED OPTION

inline po::options_description
MakeOptions( )
{
    po::options_description EXPRoptions( "CylinderStokes EXPR options" );
        EXPRoptions.add_options( )
        ( "MeanVel", po::value<double>( )->default_value( 1. ),
          "Mean Velovity of the Fluid Flow" )
        ( "Radius", po::value<double>( )->default_value( 0.5 ),
          "Radius of the Domain Cylinder" );
    return EXPRoptions;
}
inline po::options_description
MakeLibOptions( )
{
    po::options_description LIBoptions( "CylinderStokes LIB options" );
        LIBoptions.add( backend_options( "PBalpha2" ) )
                  .add( backend_options( "PBQ2" ) );
        LIBoptions.add( backend_options( "PBalpha0" ) )
                  .add( backend_options( "PBQ0" ) );
        LIBoptions.add( backend_options( "PBksi2" ) );
    return LIBoptions.add( feel_options( ) );
}


// MAIN PROGRAM

int main(int argc, char** argv)
{
    // ENVIRONMENT
    Environment env( _argc=argc,
                     _argv=argv,
                     _desc=MakeOptions( ),
                     _desc_lib=MakeLibOptions( ),
                     _about=about( _name="CylinderStokes",
                                   _author="B.Surowiec",
                                   _email="benjamin.surowiec@plasticomnium.com" ) );

    // MESH LOADING
    auto Th = loadMesh( new MeshDim );

    // FINITE ELEMENT SPACE
    auto Xh = FESpaceP2P1::New( Th );
    auto Vh = FESpaceVP1::New( Th );
    //auto Zh = FESpaceVP1P0::New( Th );
    auto Zh = FESpaceVP2P1::New( Th );

    // UNKOWN VARIABLE
    auto PSI2 = Xh->element( );
    auto psi2 = PSI2.template element<0>( );
    auto prpsi2 = PSI2.template element<1>( );
    auto PSI0 = Xh->element( );
    auto psi0 = PSI0.template element<0>( );
    auto prpsi0 = PSI0.template element<1>( );
    auto Q2 = Vh->element( );
    auto Q0 = Vh->element( );
    auto KSI2 = Zh->element( );
    auto ksi2 = KSI2.template element<0>( );
    auto prksi2 = KSI2.template element<1>( );

    // TEST FUNCTION
    auto V = Xh->element( );
    auto v = V.template element<0>( );
    auto prv = V.template element<1>( );
    auto W = Vh->element( );
    auto Y = Zh->element( );
    auto y = Y.template element<0>( );
    auto pry = Y.template element<1>( );

    // RIGHT HAND SIDE : F, ALPHA0, ALPHA2
    auto f = expr( soption( "functions.f" ) );
    auto alpha2 = expr( soption( "functions.g" ) );
    auto alpha0 = expr( soption( "functions.h" ) );

    // LINEAR FORM
    auto la2 = form1( _test=Xh );
    la2 = integrate( _range=elements( Th ),
                     _expr=trans( f )*id( v ) );
    la2 += integrate( _range=markedfaces( Th, 1 ),
                      _expr=-alpha2*id( v ) );
    la2 += integrate( _range=markedfaces( Th, 2 ),
                      _expr=alpha2*id( v ) );
    auto la0 = form1( _test=Xh );
    la0 = integrate( _range=elements( Th ),
                     _expr=trans( f )*id( v ) );
    la0 += integrate( _range=markedfaces( Th, 1 ),
                      _expr=-alpha0*id( v ) );
    la0 += integrate( _range=markedfaces( Th, 2 ),
                      _expr=alpha0*id( v ) );

    // BILINEAR FORM
    auto aa2 = form2( _trial=Xh, _test=Xh );
    aa2 = integrate( _range=elements( Th ),
                     _expr=inner( gradt( psi2 ),grad( v ) ) );
    aa2 += integrate( _range=elements( Th ),
                      _expr=idt( psi2 )*id( prv ) );
    aa2 += integrate( _range=elements( Th ),
                      _expr=id( v )*idt( prpsi2 ) );
    auto aa0 = form2( _trial=Xh, _test=Xh );
    aa0 = integrate( _range=elements( Th ),
                     _expr=inner( gradt( psi0 ), grad( v ) ) );
    aa0 += integrate( _range=elements( Th ),
                      _expr=idt( psi0 )*id( prv ) );
    aa0 += integrate( _range=elements( Th ),
                      _expr=id( v )*idt( prpsi0 ) );

    // SOLVING PROBLEM
    aa2.solve( _name="PBalpha2",
               _rhs=la2,
               _solution=PSI2 );
    aa0.solve( _name="PBalpha0",
               _rhs=la0,
               _solution=PSI0 );

    // LINEAR FORM
    auto lq2 = form1( _test=Vh );
    lq2 = integrate( _range=elements( Th ),
                     _expr=inner( trans( gradv( psi2 ) ),id( W ) ) );
    auto lq0 = form1( _test=Vh );
    lq0 = integrate( _range=elements( Th ),
                     _expr=inner( trans( gradv( psi0 ) ),id( W ) ) );

    // BILINEAR FORM
    auto aq2 = form2( _trial=Vh, _test=Vh );
    aq2 = integrate( _range=elements( Th ),
                     _expr=inner( idt( Q2 ), id( W ) ) );
    auto aq0 = form2( _trial=Vh, _test=Vh );
    aq0 = integrate( _range=elements( Th ),
                     _expr=inner( idt( Q0 ), id( W ) ) );

    // SOLVING PROBLEM
    aq2.solve( _name="PBQ2",
               _rhs=lq2,
               _solution=Q2 );
    aq0.solve( _name="PBQ0",
               _rhs=lq0,
               _solution=Q0 );

    ///////////////////////////////////////////////////////////////////////
    //                   JUSQUE LA... TOUT VA BIEN :-D                   //
    ///////////////////////////////////////////////////////////////////////

    // LINEAR FORM
    auto lk2 = form1( _test=Zh );
    lk2 = integrate( _range=elements( Th ),
                     _expr=inner( idv( Q2 ), id( y ) ) );
    //lk2 += integrate( _range=elements( Th ),
    //                  _expr=1.0e-16*idt( prksi2 )*id( pry ) );

    ///////////////////////////////////////////////////////////////////////
    //                   JUSQUE LA... MOUAIS PAS SUR                     //
    ///////////////////////////////////////////////////////////////////////

    // BILINEAR FORM
    auto ak2 = form2( _trial=Zh, _test=Zh );
    ak2 = integrate( _range=elements( Th ),
                     _expr=inner( gradt( ksi2 ), grad( y ) ) );
    ak2 += integrate( _range=elements( Th ),
                      _expr=divt( ksi2 )*id( pry ) );
    ak2 += integrate( _range=elements( Th ),
                      _expr=idt( prksi2 )*div( y ) );
    ak2 += on( _range=boundaryfaces( Th ),
               _rhs=lk2,
               _element=ksi2,
               _expr=vec( cst( 0. ), cst( 0. ), cst( 0. ) ) );

    // SOLVING PROBLEM
    ak2.solve( _name="PBksi2",
               _rhs=lk2,
               _solution=KSI2 );

    // RESULT STORAGE
    auto e = exporter( Th );
    e->add( "psi2", psi2 );
    e->add( "psi0", psi0 );
    e->add( "Q2", Q2 );
    e->add( "Q0", Q0 );
    e->add( "ksi2", ksi2 );
    e->save( );

    // END OF MAIN PROGRAM
    if(Environment::isMasterRank())
        std::cout << "THE END" << std::endl;
    return 0;
}
