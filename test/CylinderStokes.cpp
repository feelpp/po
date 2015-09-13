// MANDATORY INCLUDE AND USED SPACE
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/dh.hpp>
using namespace Feel;
using namespace Feel::vf;

// TYPE DEFINITON
#define FEELPP_DIM 3
typedef Mesh<Simplex<FEELPP_DIM>> MeshDim;
typedef Lagrange<0,Scalar> LagrP0;
typedef Lagrange<1,Scalar> LagrP1;
typedef Lagrange<2,Scalar> LagrP2;
typedef Lagrange<3,Scalar> LagrP3;
typedef Lagrange<2,Vectorial> LagrVecP2;
typedef RaviartThomas<0> RaTh0;
typedef bases<LagrP1,LagrP0> BaseP1P0;
typedef bases<LagrP3,LagrP2> BaseP3P2;
typedef bases<LagrVecP2> BaseVP2;
typedef bases<LagrVecP2,LagrP0> BaseVP2P0;
typedef bases<RaTh0,LagrP1> BaseRT0P1;
typedef FunctionSpace<MeshDim,BaseP1P0> FESpaceP1P0;
typedef FunctionSpace<MeshDim,BaseP3P2> FESpaceP3P2;
typedef FunctionSpace<MeshDim,BaseVP2> FESpaceVP2;
typedef FunctionSpace<MeshDim,BaseVP2P0> FESpaceVP2P0;
typedef FunctionSpace<MeshDim,BaseRT0P1> FESpaceRT0P1;

// EXTENDED OPTION
inline po::options_description
MakeOptions( )
{
    po::options_description PERSOoptions( "CylinderStokes PERSO options" );
    PERSOoptions.add_options( )
        ( "Radius", po::value<double>( )->default_value( 0.05 ), "Radius of the Cylinder Domain" )
        ( "MeanVel", po::value<double>( )->default_value( 1. ), "Mean Velovity of the Fluid Flow" )
        ( "alift.avar1", po::value<bool>()->default_value( true ), "Lap+Grad+Eq PB : 3steps" )
        ( "alift.avar2", po::value<bool>()->default_value( false ), "Lap+Grad PB : 2steps" )
        ( "alift.avar3", po::value<bool>()->default_value( false ), "Mixed PB" )
        ( "alift.lph2", po::value<bool>()->default_value( true ), "Alpha 2 considered" )
        ( "alift.lph1", po::value<bool>()->default_value( false ), "Alpha 1 considered" )
        ( "alift.lph0", po::value<bool>()->default_value( true ), "Alpha 0 considered" )
        ( "alift.assy", po::value<bool>()->default_value( true ), "A Lift Recombination" )
        ( "alift.cvdisp", po::value<bool>()->default_value( true ), "Convergence Display" )
        ( "alift.store", po::value<bool>()->default_value( true ), "Storage" );
    return PERSOoptions;
}
inline po::options_description
MakeLibOptions( )
{
    po::options_description LIBoptions( "CylinderStokes LIB options" );
    LIBoptions.add( backend_options( "PBpsi2" ) );
    LIBoptions.add( backend_options( "PBQ2" ) );
    LIBoptions.add( backend_options( "PBksi2" ) );
    LIBoptions.add( backend_options( "PBpsi0" ) );
    LIBoptions.add( backend_options( "PBQ0" ) );
    LIBoptions.add( backend_options( "PBAA" ) );
    return LIBoptions.add( feel_options( ) );
}

// MAIN PROGRAM
int main(int argc, char** argv)
{
    // ENVIRONMENT
    Environment env( _argc=argc, _argv=argv,
                     _desc=MakeOptions( ), _desc_lib=MakeLibOptions( ),
                     _about=about( _name="CylinderStokes", _author="B.Surowiec", _email="benjamin.surowiec@plasticomnium.com" ) );

    // BOULEAN MANAGEMENT
    auto SCRWRT = Environment::isMasterRank( );
    auto AV1 = boption( "alift.avar1" );
    auto AV2 = boption( "alift.avar2" );
    auto AV3 = boption( "alift.avar3" );
    auto LPH2 = boption( "alift.lph2" );
    auto LPH1 = boption( "alift.lph1" );
    auto LPH0 = boption( "alift.lph0" );
    auto ASSY = boption( "alift.assy" );
    auto CVDISP = boption( "alift.cvdisp" );
    auto STORE = boption( "alift.store" );

    // OPTIONS GUARDRAIL
    if( (AV1&&AV2&&AV3)||(AV1&&AV2&&!AV3)||(AV1&&!AV2&&AV3)||(!AV1&&AV2&&AV3)||(!AV1&&!AV2&&!AV3) ) { AV1=1; AV2 = 0; AV3 = 0; }

    // MESH LOADING AND FINITE ELEMENT SPACE
    auto Th = loadMesh( new MeshDim );
    auto XhP1P0 = FESpaceP1P0::New( Th );
    auto XhP3P2 = FESpaceP3P2::New( Th );
    auto XhVP2 = FESpaceVP2::New( Th );
    auto XhVP2P0 = FESpaceVP2P0::New( Th );
    auto XhRT0P1 = FESpaceRT0P1::New( Th );

    // UNKOWN VARIABLE
    auto QP2w = XhRT0P1->element( );   auto Q2nw = QP2w.template element<0>( );   auto PSI2nw = QP2w.template element<1>( );
    auto QP0w = XhRT0P1->element( );   auto Q0nw = QP0w.template element<0>( );   auto PSI0nw = QP0w.template element<1>( );
    auto PSI2 = XhP1P0->element( );    auto psi2 = PSI2.template element<0>( );   auto prpsi2 = PSI2.template element<1>( );
    auto KSI2 = XhVP2P0->element( );   auto ksi2 = KSI2.template element<0>( );   auto prksi2 = KSI2.template element<1>( );
    auto PSI0 = XhP3P2->element( );    auto psi0 = PSI0.template element<0>( );   auto prpsi0 = PSI0.template element<1>( );
    auto Q2 = XhVP2->element( );
    auto Q0 = XhVP2->element( );
    auto AA = XhVP2->element( );

    // TEST FUNCTION
    auto ZZZZ = XhRT0P1->element( );   auto zzzz = ZZZZ.template element<0>( );   auto przzzz = ZZZZ.template element<1>( );
    auto VVVV = XhP1P0->element( );    auto vvvv = VVVV.template element<0>( );   auto prvvvv = VVVV.template element<1>( );
    auto YYYY = XhVP2P0->element( );   auto yyyy = YYYY.template element<0>( );   auto pryyyy = YYYY.template element<1>( );
    auto XXXX = XhP3P2->element( );    auto xxxx = XXXX.template element<0>( );   auto prxxxx = XXXX.template element<1>( );
    auto WW = XhVP2->element( );

    // KNOWN DATA : F, ALPHA0, ALPHA1, ALPHA2
    //auto rhsf = expr<3,1>( soption( "functions.f" ) );
    auto alpha2 = expr( soption( "functions.g" ) );
    //auto alpha1 = expr( soption( "functions.i" ) );
    auto alpha0 = expr( soption( "functions.h" ) );

    // GENERAL PARAMETERS & CONVERGENCE DATA
    auto precis = 4;
    auto delta = 0.5;
    auto valdivaa = normL2( _range=elements( Th ), _expr=cst( 0. ) );
    auto valdivasum = normL2( _range=elements( Th ), _expr=cst( 0. ) );
    auto valdivq0ksi2 = normL2( _range=elements( Th ), _expr=cst( 0. ) );
    auto valdivq2 = normL2( _range=elements( Th ), _expr=cst( 0. ) );
    auto valdivksi2 = normL2( _range=elements( Th ), _expr=cst( 0. ) );
    auto valdivq0 = normL2( _range=elements( Th ), _expr=cst( 0. ) );

    // START OF CALCULATIONS
    if( SCRWRT )
    {
        std::cout << std::endl << "    INITIALIZATION DONE !!" << std::endl;
        std::cout << std::endl << "    CURRENT RUNNING OPTIONS ARE :" << std::endl;
        std::cout << "    -----------------------------------------------------" << std::endl;
        std::cout << "    VARIANTS : " << "Variant 1" << "\t" << "Variant 2" << "\t" << "Variant 3" << std::endl;
        std::cout << "    VALUES   : " << AV1 << "\t\t" << AV2 << "\t\t" << AV3 << std::endl;
        std::cout << "    -----------------------------------------------------" << std::endl;
        std::cout << "    LIFTS    : " << "Alpha 2  " << "\t" << "Alpha 1  " << "\t" << "Alpha 0  " << std::endl;
        std::cout << "    VALUES   : " << LPH2 << "\t\t" << LPH1 << "\t\t" << LPH0 << std::endl;
        std::cout << "    -----------------------------------------------------" << std::endl;
        std::cout << std::endl << "    JOB IS CURRENTLY RUNNING..." << std::endl << std::endl;
    }

    // LIFT VARIANT 1
    if( AV1 && (LPH2||LPH1||LPH0) )
    {
        // Lift Alpha2
        if( LPH2 )
        {
            auto la2 = form1( _test=XhP1P0 );
            la2 = integrate( _range=markedfaces( Th, 1 ), _expr=-alpha2*id( vvvv ) );
            la2 += integrate( _range=markedfaces( Th, 2 ), _expr=alpha2*id( vvvv ) );
            auto aa2 = form2( _trial=XhP1P0, _test=XhP1P0 );
            aa2 = integrate( _range=elements( Th ), _expr=inner( gradt( psi2 ),grad( vvvv ) ) );
            aa2 += integrate( _range=elements( Th ), _expr=idt( psi2 )*id( prvvvv ) );
            aa2 += integrate( _range=elements( Th ), _expr=id( vvvv )*idt( prpsi2 ) );
            aa2.solve( _name="PBpsi2", _rhs=la2, _solution=PSI2 );
            if( SCRWRT ) std::cout << "    PSI2 SOLVED !!" << std::endl << "    [...]" << std::endl;

            auto lq2 = form1( _test=XhVP2 );
            lq2 = integrate( _range=elements( Th ), _expr=inner( trans( gradv( psi2 ) ),id( WW ) ) );
            auto aq2 = form2( _trial=XhVP2, _test=XhVP2 );
            aq2 = integrate( _range=elements( Th ), _expr=inner( idt( Q2 ), id( WW ) ) );
            aq2 += on( _range=markedfaces( Th, 1 ), _rhs=lq2, _element=Q2, _expr=vec( cst( 0. ), cst( 0. ), alpha2 ) );
            aq2 += on( _range=markedfaces( Th, 2 ), _rhs=lq2, _element=Q2, _expr=vec( cst( 0. ), cst( 0. ), alpha2 ) );
            aq2.solve( _name="PBQ2", _rhs=lq2, _solution=Q2 );
            if( SCRWRT ) std::cout << "    Q2   SOLVED !!" << std::endl << "    [...]" << std::endl;

            auto lk2 = form1( _test=XhVP2P0 );
            lk2 = integrate( _range=elements( Th ), _expr=-inner( idv( Q2 ), id( yyyy ) ) );
            auto ak2 = form2( _trial=XhVP2P0, _test=XhVP2P0 );
            ak2 = integrate( _range=elements( Th ), _expr=inner( gradt( ksi2 ), grad( yyyy ) ) );
            ak2 += integrate( _range=elements( Th ), _expr=-divt( ksi2 )*id( pryyyy ) );
            ak2 += integrate( _range=elements( Th ), _expr=idt( prksi2 )*div( yyyy ) );
            ak2 += on( _range=boundaryfaces( Th ), _rhs=lk2, _element=ksi2, _expr=vec( cst( 0. ), cst( 0. ), cst( 0. ) ) );
            ak2.solve( _name="PBksi2", _rhs=lk2, _solution=KSI2 );
            if( SCRWRT ) std::cout << "    KSI2 SOLVED !!" << std::endl << "    [...]" << std::endl;
        }

        // Lift Alpha1
        if( LPH1 ) { if( SCRWRT ) std::cout << "    Alpha 1 Not Considered ATM !!" << std::endl << "    [...]" << std::endl; }

        // Lift Alpha0
        if( LPH0 )
        {
            auto la0 = form1( _test=XhP3P2 );
            la0 = integrate( _range=markedfaces( Th, 1 ), _expr=-alpha0*id( xxxx ) );
            la0 += integrate( _range=markedfaces( Th, 2 ), _expr=alpha0*id( xxxx ) );
            auto aa0 = form2( _trial=XhP3P2, _test=XhP3P2 );
            aa0 = integrate( _range=elements( Th ), _expr=inner( gradt( psi0 ),grad( xxxx ) ) );
            aa0 += integrate( _range=elements( Th ), _expr=idt( psi0 )*id( prxxxx ) );
            aa0 += integrate( _range=elements( Th ), _expr=id( xxxx )*idt( prpsi0 ) );
            aa0.solve( _name="PBpsi0", _rhs=la0, _solution=PSI0 );
            if( SCRWRT ) std::cout << "    PSI0 SOLVED !!" << std::endl << "    [...]" << std::endl;

            auto lq0 = form1( _test=XhVP2 );
            lq0 = integrate( _range=elements( Th ), _expr=inner( trans( gradv( psi0 ) ),id( WW ) ) );
            auto aq0 = form2( _trial=XhVP2, _test=XhVP2 );
            aq0 = integrate( _range=elements( Th ), _expr=inner( idt( Q0 ), id( WW ) ) );
            aq0 += on( _range=markedfaces( Th, 1 ), _rhs=lq0, _element=Q0, _expr=vec( cst( 0. ), cst( 0. ), alpha0 ) );
            aq0 += on( _range=markedfaces( Th, 2 ), _rhs=lq0, _element=Q0, _expr=vec( cst( 0. ), cst( 0. ), alpha0 ) );
            aq0.solve( _name="PBQ0", _rhs=lq0, _solution=Q0 );
            if( SCRWRT ) std::cout << "    Q0   SOLVED !!" << std::endl << "    [...]" << std::endl;
        }

        // Assembly Of AA
        if( ASSY && (LPH2||LPH1||LPH0) )
        {
            auto laa = form1( _test=XhVP2 );
            laa = integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH2 ) laa += integrate( _range=elements( Th ), _expr=inner( idv( ksi2 ), id( WW ) ) );
            if( LPH1 ) laa += integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH0 ) laa += integrate( _range=elements( Th ), _expr=inner( idv( Q0 ), id( WW ) ) );
            auto aaa = form2( _trial=XhVP2, _test=XhVP2 );
            aaa = integrate( _range=elements( Th ), _expr=inner( idt( AA ), id( WW ) ) );
            aaa.solve( _name="PBAA", _rhs=laa, _solution=AA );
            if( SCRWRT ) std::cout << "    AA   RECOMBINED !!" << std::endl << "    [...]" << std::endl;
        }

        // Convergence Result Display
        if( CVDISP && (LPH2||LPH1||LPH0) )
        {
            if( LPH2 ) valdivq2 = normL2( _range=elements( Th ), _expr=divv( Q2 ) );
            if( LPH2 ) valdivksi2 = normL2( _range=elements( Th ), _expr=divv( ksi2 ) );
            if( LPH0 ) valdivq0 = normL2( _range=elements( Th ), _expr=divv( Q0 ) );
            if( LPH2 && LPH0 ) valdivq0ksi2 = valdivq0 + valdivksi2;
            if( LPH2 && LPH0 ) valdivasum = normL2( _range=elements( Th ), _expr=divv( ksi2 ) + divv( Q0 ) );
            if( ASSY ) valdivaa = normL2( _range=elements( Th ), _expr=divv( AA ) );
            if( SCRWRT )
            {
                std::cout.precision(precis);
                std::cout << std::endl << "    CONVERGENCE RESULTS :" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 ) std::cout << "    "  << std::scientific << valdivq2 << " = ||div(Q2)||_L2" << std::endl;
                if( LPH2 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 ) std::cout << "    " << std::scientific << valdivksi2 << " = ||div(ksi2)||_L2" << std::endl;
                if( LPH2 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH0 ) std::cout << "    " << std::scientific << valdivq0 << " = ||div(Q0)||_L2" << std::endl;
                if( LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    " << std::scientific << valdivq0ksi2 << " = ||div(Q0)||_L2 + ||div(ksi2)||_L2" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    " << std::scientific << valdivasum << " = ||div(Q0) + div(ksi2)||_L2" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( ASSY ) std::cout << "    " << std::scientific << valdivaa << " = ||div(AA)||_L2" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
            }
         }

        // Result Storage
        if( STORE && (LPH2||LPH1||LPH0) )
        {
            auto e = exporter( Th );
            if( LPH2 ) e->add( "psi2", psi2 );  e->add( "Q2", Q2 );  e->add( "ksi2", ksi2 );
            if( LPH0 ) e->add( "psi0", psi0 );  e->add( "Q0", Q0 );
            if( ASSY ) e->add( "AA", AA );
            e->save( );
            if( SCRWRT ) std::cout << std::endl << "    RESULTS HAVE BEEN STORED !! (PARAVIEW FORMAT)" << std::endl;
         }
    }

    // LIFT VARIANT 2
    if( AV2 && (LPH2||LPH1||LPH0) )
    {
        // Lift Alpha2
        if( LPH2 )
        {
            auto la2 = form1( _test=XhP1P0 );
            la2 = integrate( _range=markedfaces( Th, 1 ), _expr=-alpha2*id( vvvv ) );
            la2 += integrate( _range=markedfaces( Th, 2 ), _expr=alpha2*id( vvvv ) );
            auto aa2 = form2( _trial=XhP1P0, _test=XhP1P0 );
            aa2 = integrate( _range=elements( Th ), _expr=inner( gradt( psi2 ),grad( vvvv ) ) );
            aa2 += integrate( _range=elements( Th ), _expr=idt( psi2 )*id( prvvvv ) );
            aa2 += integrate( _range=elements( Th ), _expr=id( vvvv )*idt( prpsi2 ) );
            aa2.solve( _name="PBpsi2", _rhs=la2, _solution=PSI2 );
            if( SCRWRT ) std::cout << "    PSI2 SOLVED !!" << std::endl << "    [...]" << std::endl;

            auto lk2 = form1( _test=XhVP2P0 );
            lk2 = integrate( _range=elements( Th ), _expr=-inner( trans( gradv( psi2 ) ), id( yyyy ) ) );
            auto ak2 = form2( _trial=XhVP2P0, _test=XhVP2P0 );
            ak2 = integrate( _range=elements( Th ), _expr=inner( gradt( ksi2 ), grad( yyyy ) ) );
            ak2 += integrate( _range=elements( Th ), _expr=-divt( ksi2 )*id( pryyyy ) );
            ak2 += integrate( _range=elements( Th ), _expr=idt( prksi2 )*div( yyyy ) );
            ak2 += on( _range=boundaryfaces( Th ), _rhs=lk2, _element=ksi2, _expr=vec( cst( 0. ), cst( 0. ), cst( 0. ) ) );
            ak2.solve( _name="PBksi2", _rhs=lk2, _solution=KSI2 );
            if( SCRWRT ) std::cout << "    KSI2 SOLVED !!" << std::endl << "    [...]" << std::endl;
        }

        // Lift Alpha1
        if( LPH1 ) { if( SCRWRT ) std::cout << "    Alpha 1 Not Considered ATM !!" << std::endl << "    [...]" << std::endl; }

        // Lift Alpha0
        if( LPH0 )
        {
            auto la0 = form1( _test=XhP3P2 );
            la0 = integrate( _range=markedfaces( Th, 1 ), _expr=-alpha0*id( xxxx ) );
            la0 += integrate( _range=markedfaces( Th, 2 ), _expr=alpha0*id( xxxx ) );
            auto aa0 = form2( _trial=XhP3P2, _test=XhP3P2 );
            aa0 = integrate( _range=elements( Th ), _expr=inner( gradt( psi0 ),grad( xxxx ) ) );
            aa0 += integrate( _range=elements( Th ), _expr=idt( psi0 )*id( prxxxx ) );
            aa0 += integrate( _range=elements( Th ), _expr=id( xxxx )*idt( prpsi0 ) );
            aa0.solve( _name="PBpsi0", _rhs=la0, _solution=PSI0 );
            if( SCRWRT ) std::cout << "    PSI0 SOLVED !!" << std::endl << "    [...]" << std::endl;
        }

        // Assembly Of AA
        if( ASSY && (LPH2||LPH1||LPH0) )
        {
            auto laa = form1( _test=XhVP2 );
            laa = integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH2 ) laa += integrate( _range=elements( Th ), _expr=inner( idv( ksi2 ), id( WW ) ) );
            if( LPH1 ) laa += integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH0 ) laa += integrate( _range=elements( Th ), _expr=inner( trans( gradv( psi0 ) ), id( WW ) ) );
            auto aaa = form2( _trial=XhVP2, _test=XhVP2 );
            aaa = integrate( _range=elements( Th ), _expr=inner( idt( AA ), id( WW ) ) );
            aaa.solve( _name="PBAA", _rhs=laa, _solution=AA );
            if( SCRWRT ) std::cout << "    AA   RECOMBINED !!" << std::endl << "    [...]" << std::endl;
        }

        // Convergence Result Display
        if( CVDISP && (LPH2||LPH1||LPH0) )
        {
            if( LPH2 ) valdivksi2 = normL2( _range=elements( Th ), _expr=divv( ksi2 ) );
            if( ASSY ) valdivaa = normL2( _range=elements( Th ), _expr=divv( AA ) );
            if( SCRWRT )
            {
                std::cout.precision(precis);
                std::cout << std::endl << "    CONVERGENCE RESULTS :" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 ) std::cout << "    " << std::scientific << valdivksi2 << " = ||div(ksi2)||_L2" << std::endl;
                if( LPH2 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( ASSY ) std::cout << "    " << std::scientific << valdivaa << " = ||div(AA)||_L2" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
            }
         }

        // Result Storage
        if( STORE && (LPH2||LPH1||LPH0) )
        {
            auto e = exporter( Th );
            if( LPH2 ) e->add( "psi2", psi2 );  e->add( "ksi2", ksi2 );
            if( LPH0 ) e->add( "psi0", psi0 );
            if( ASSY ) e->add( "AA", AA );
            e->save( );
            if( SCRWRT ) std::cout << std::endl << "    RESULTS HAVE BEEN STORED !! (PARAVIEW FORMAT)" << std::endl;
         }
     }

    // LIFT VARIANT 3
    if( AV3 && (LPH2||LPH1||LPH0) )
    {
        // Lift Alpha2
        if( LPH2 )
        {
            auto la2q2 = form1( _test=XhRT0P1 );
            la2q2 = integrate( _range=elements( Th ), _expr=cst( 0. ) );
            auto aa2q2 = form2( _trial=XhRT0P1, _test=XhRT0P1 );
            aa2q2 = integrate( _range=elements( Th ), _expr=inner( idt( Q2nw ), id( zzzz ) ) );
            aa2q2 += integrate( _range=elements( Th ), _expr=idt( PSI2nw )*div( zzzz ) );
            aa2q2 += integrate( _range=elements( Th ), _expr=divt( Q2nw )*id( przzzz ) );
            aa2q2 += integrate( _range=elements( Th ), _expr=-delta*inner( idt( Q2nw )-trans( gradt( PSI2nw ) ), id( zzzz )+trans( grad( przzzz ) ) ) );
            aa2q2 += integrate( _range=elements( Th ), _expr=-delta*divt( Q2nw )*div( zzzz ) );
            aa2q2 += on( _range=markedfaces( Th, 1 ), _rhs=la2q2, _element=Q2nw, _expr=vec( cst( 0. ), cst( 0. ), alpha2 ) );
            aa2q2 += on( _range=markedfaces( Th, 2 ), _rhs=la2q2, _element=Q2nw, _expr=vec( cst( 0. ), cst( 0. ), alpha2 ) );
            aa2q2.solve( _name="PBQ2", _rhs=la2q2, _solution=QP2w );
            if( SCRWRT ) std::cout << "    Q2   SOLVED !!" << std::endl << "    [...]" << std::endl;

            auto lk2 = form1( _test=XhVP2P0 );
            lk2 = integrate( _range=elements( Th ), _expr=-inner( trans( gradv( psi2 ) ), id( yyyy ) ) );
            auto ak2 = form2( _trial=XhVP2P0, _test=XhVP2P0 );
            ak2 = integrate( _range=elements( Th ), _expr=inner( gradt( ksi2 ), grad( yyyy ) ) );
            ak2 += integrate( _range=elements( Th ), _expr=-divt( ksi2 )*id( pryyyy ) );
            ak2 += integrate( _range=elements( Th ), _expr=idt( prksi2 )*div( yyyy ) );
            ak2 += on( _range=boundaryfaces( Th ), _rhs=lk2, _element=ksi2, _expr=vec( cst( 0. ), cst( 0. ), cst( 0. ) ) );
            ak2.solve( _name="PBksi2", _rhs=lk2, _solution=KSI2 );
            if( SCRWRT ) std::cout << "    KSI2 SOLVED !!" << std::endl << "    [...]" << std::endl;
         }

        // Lift Alpha1
        if( LPH1 ) { if( SCRWRT ) std::cout << "    Alpha 1 Not Considered ATM !!" << std::endl << "    [...]" << std::endl; }

        // Lift Alpha0
        if( LPH0 )
        {
            auto la0q0 = form1( _test=XhRT0P1 );
            la0q0 = integrate( _range=elements( Th ), _expr=cst( 0. ) );
            auto aa0q0 = form2( _trial=XhRT0P1, _test=XhRT0P1 );
            aa0q0 = integrate( _range=elements( Th ), _expr=inner( idt( Q0nw ), id( zzzz ) ) );
            aa0q0 += integrate( _range=elements( Th ), _expr=idt( PSI0nw )*div( zzzz ) );
            aa0q0 += integrate( _range=elements( Th ), _expr=divt( Q0nw )*id( przzzz ) );
            aa0q0 += integrate( _range=elements( Th ), _expr=-delta*inner( idt( Q0nw )-trans( gradt( PSI0nw ) ), id( zzzz )+trans( grad( przzzz ) ) ) );
            aa0q0 += integrate( _range=elements( Th ), _expr=-delta*divt( Q0nw )*div( zzzz ) );
            aa0q0 += on( _range=markedfaces( Th, 1 ), _rhs=la0q0, _element=Q0nw, _expr=vec( cst( 0. ), cst( 0. ), alpha0 ) );
            aa0q0 += on( _range=markedfaces( Th, 2 ), _rhs=la0q0, _element=Q0nw, _expr=vec( cst( 0. ), cst( 0. ), alpha0 ) );
            aa0q0.solve( _name="PBQ0", _rhs=la0q0, _solution=QP0w );
            if( SCRWRT ) std::cout << "    Q0   SOLVED !!" << std::endl << "    [...]" << std::endl;
         }

        // Assembly Of AA
        if( ASSY && (LPH2||LPH1||LPH0) )
        {
            auto laa = form1( _test=XhVP2 );
            laa = integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH2 ) laa += integrate( _range=elements( Th ), _expr=inner( idv( ksi2 ), id( WW ) ) );
            if( LPH1 ) laa += integrate( _range=elements( Th ), _expr=cst( 0. ) );
            if( LPH0 ) laa += integrate( _range=elements( Th ), _expr=inner( idv( Q0 ), id( WW ) ) );
            auto aaa = form2( _trial=XhVP2, _test=XhVP2 );
            aaa = integrate( _range=elements( Th ), _expr=inner( idt( AA ), id( WW ) ) );
            aaa.solve( _name="PBAA", _rhs=laa, _solution=AA );
            if( SCRWRT ) std::cout << "    AA   RECOMBINED !!" << std::endl << "    [...]" << std::endl;
         }

        // Convergence Result Display
        if( CVDISP && (LPH2||LPH1||LPH0) )
        {
            if( LPH2 ) valdivq2 = normL2( _range=elements( Th ), _expr=divv( Q2 ) );
            if( LPH2 ) valdivksi2 = normL2( _range=elements( Th ), _expr=divv( ksi2 ) );
            if( LPH0 ) valdivq0 = normL2( _range=elements( Th ), _expr=divv( Q0 ) );
            if( LPH2 && LPH0 ) valdivq0ksi2 = valdivq0 + valdivksi2;
            if( LPH2 && LPH0 ) valdivasum = normL2( _range=elements( Th ), _expr=divv( ksi2 ) + divv( Q0 ) );
            if( ASSY ) valdivaa = normL2( _range=elements( Th ), _expr=divv( AA ) );
            if( SCRWRT )
            {
                std::cout.precision(precis);
                std::cout << std::endl << "    CONVERGENCE RESULTS :" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 ) std::cout << "    "  << std::scientific << valdivq2 << " = ||div(Q2)||_L2" << std::endl;
                if( LPH2 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 ) std::cout << "    " << std::scientific << valdivksi2 << " = ||div(ksi2)||_L2" << std::endl;
                if( LPH2 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH0 ) std::cout << "    " << std::scientific << valdivq0 << " = ||div(Q0)||_L2" << std::endl;
                if( LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    " << std::scientific << valdivq0ksi2 << " = ||div(Q0)||_L2 + ||div(ksi2)||_L2" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    " << std::scientific << valdivasum << " = ||div(Q0) + div(ksi2)||_L2" << std::endl;
                if( LPH2 && LPH0 ) std::cout << "    ----------------------------------------------" << std::endl;
                if( ASSY ) std::cout << "    " << std::scientific << valdivaa << " = ||div(AA)||_L2" << std::endl;
                std::cout << "    ----------------------------------------------" << std::endl;
            }
         }

        // Result Storage
        if( STORE && (LPH2||LPH1||LPH0) )
        {
            auto e = exporter( Th );
            if( LPH2 ) e->add( "psi2", psi2 );  e->add( "Q2", Q2 );  e->add( "ksi2", ksi2 );
            if( LPH0 ) e->add( "psi0", psi0 );  e->add( "Q0", Q0 );
            if( ASSY ) e->add( "AA", AA );
            e->save( );
            if( SCRWRT ) std::cout << std::endl << "    RESULTS HAVE BEEN STORED !! (PARAVIEW FORMAT)" << std::endl;
         }
    }

    // END OF MAIN PROGRAM
    if( SCRWRT )
        std::cout << std::endl << "    END OF JOB : CAN NOW BE KILLED !!" << std::endl << std::endl;
    return 0;
}

