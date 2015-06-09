// MANDATORY INCLUDE AND USED SPACE

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

// TYPE DEFINITON

typedef Mesh<Simplex<FEELPP_DIM> > MeshDim;
typedef bases<Lagrange<1,Scalar> > BaseP1;
typedef bases<Lagrange<0,Scalar> > BaseP0;
typedef bases<BaseP1,BaseP0> BaseP1P0;
typedef FunctionSpace<MeshDim,bases<BaseP1>> FESpaceP1;
typedef FunctionSpace<MeshDim,BaseP1P0> FESpaceP1P0;

// EXTENDED OPTION

inline po::po_options_description
MakeOptions( )
{
  po:po_options_description EXPRoptions( "BBB optioons" );
  EXPRoptions.add_options( )
      ( "MeanVel", po::value<double>( ) -> default_value( 1. ), "Mean Velovity of the Flow" )
      ( "Radius", po::value<double>( ) -> default_value( 0.5 ), "Radius of the Cylinder" )
      ;
  return EXPRoptions;
}

// MAIN PROGRAM

int main(int argc, char** argv)
{

    // ENVIRONMENT
    Environment env( _argc=argc,
                     _argv=argv,
                     //_desc=feel_options( ),
                     // MakeOptions( ) = feel_options( ) + specific options
                     _desc=MakeOptions( ),
                     _about=about( _name="BBB",
                                   _author="B.Surowiec",
                                   _email="benjamin.surowiec@plasticomnium.com" ) );

    // MESH LOADING
    auto Th = loadmesh( new MeshDim );

    // FINITE ELEMENT SPACE
    auto XVh2 = FESpaceP1P0::New( Th );

    // UNKOWN VARIABLE
    auto U = XVh2 -> Element( );
    auto u = U.element<0>( );
    auto p = U.element<1>( );

    // TEST FUNCTION
    auto V = XVh2 -> Element( );
    auto v = U.element<0>( );
    auto q = U.element<1>( );

    // BILINEAR FORM
    auto a = form2( _trial=XVh2,
                    _test=XVh2 );

    a = integrate ( _range=element( Th ),
                    // inner() perform scalar product between transposed u and v
                    // _expr=inner( gradt( u ),grad( v ) ) );
                    // 2D integral from IBP = 0 because u.n = 0
                    _expr=trans( gradt( u ) )*grad( v ) );

    a += integrate( _range=element( Th ),
                    // flow rate conservation term
                    _expr=divt( u )*id( q ) );

    a += integrate( _range=element( Th ),
                    // flow rate conservation symetrical term
                    _expr=div( v )*idt( p ) );

    // BOUNDARY CONDITION : RIGHT HAND SIDE F
    auto f = expr( soption( "function.f" ) );

    // ALPHA_0 FUNCTION
    // old method before MakeOptions( )
    //    auto h_s = soption( "function.h" );
    //    auto vars = Symbols{ "x", "y", "Radius", "MeanVel" };
    //    auto h_e = parse( h_s, vars ) ;
    //    auto alpha0 = expr( h_e, vars );
    //    alpha0.setParameterValues( { { "Radius", 0.5 },
    //                                 { "MeanVel", 1 } } );
    auto alpha0 = expr( soption( "function.h" ) );

    // ALPHA_0 VERIFICATION
    auto d = integrate( markedfaces( Th, 1 ), -alpha0 ).evaluate( )( 0, 0 );
    d += integrate( markedfaces( Th, 2 ), alpha0 ).evaluate( )( 0, 0 );
    std::cout << " alpha_0 : " << alpha0 << "\nint : " << d << std::endl;

    // LINEAR FORM
    auto l = form1( _test=XVh2 );

    l += integrate( _range=markedfaces( Th, 1 ),
                    _expr=-alpha0*id( v ) );

    l += integrate( _range=markedfaces( Th, 2 ),
                    _expr=alpha0*id( v ) );

    // SOLVING PROBLEM
    a.solve( _rhs=l,
             _solution=U );

    // RESULT STORAGE
    auto e = exporter( Th );
    e -> add( "u", u );
    e -> add( "p", p );
    e -> save( );

    // END OF MAIN PROGRAM
    return 0;
}







