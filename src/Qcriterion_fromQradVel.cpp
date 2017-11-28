#include <feel/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feelvf/vf.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description exportdboptions( "export database options" );
	exportdboptions.add_options()
        ( "ifile", po::value<std::string>(), "input database file" )
        ( "field", po::value<std::string>(), "field loaded in database " )
        ( "time-initial", po::value<double>(), "initial time used for pod" )
        ( "time-final", po::value<double>(), "final time used for pod" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=exportdboptions,
                   _about=about(_name="export_database",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string ifile;
    if ( Environment::vm().count("ifile") )
        ifile = Environment::vm()["ifile"].as<std::string>();
    std::string fieldPod;
    fieldPod = "grad_pressure";

//    if ( Environment::vm().count("field") )
//        fieldPod = Environment::vm()["field"].as<std::string>();

    if ( ifile.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile is missing\n";
        return 0;
    }
    if ( fieldPod.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --field is missing\n";
        return 0;
    }

    typedef Mesh<Simplex<3>> mesh_type;
    FeelppDatabase<mesh_type> myDb;
    myDb.setFilename( ifile );
    myDb.loadInfo();

    if ( !myDb.hasField( fieldPod ) )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because fieldPod does not found in databse\n";
        return 0;

   }

    CHECK( myDb.fieldInfo( fieldPod ).basisName() == "lagrange" ) << "only lagrange instantiated : " << myDb.fieldInfo( fieldPod ).basisName();
    CHECK( myDb.fieldInfo( fieldPod ).basisOrder() == 1 ) << "only order 1 instantiated : " << myDb.fieldInfo( fieldPod ).basisOrder();

    auto mesh = myDb.loadMesh( MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES );

    auto Xh = Pch<1>( mesh );
    auto Vh = Pchv<1>( mesh );

    auto uj = Vh->element();
//    auto Qi = Xh->element();

    auto Idiv = Div( _domainSpace = Vh, _imageSpace=Xh ) ;

    auto const& fieldsInfo = myDb.fieldsInfo();
    auto const& timeSetDb = myDb.timeSet();

    double ti = timeSetDb.front();
    double tf = timeSetDb.back();
    if ( Environment::vm().count("time-initial") )
        ti = doption(_name="time-initial");
    if ( Environment::vm().count("time-final") )
        tf = doption(_name="time-final");

    std::vector<double> timeSetIndex;
    for (int k=0; k<timeSetDb.size();++k)
    {
        if ( ( timeSetDb[k] >= (ti - 1e-9) ) &&
             ( timeSetDb[k] <= (tf + 1e-9 ) ) )
            timeSetIndex.push_back( k );
    }
    int nTimeStep = timeSetIndex.size();



    std::string geoExportType = "static";
    auto e = exporter( _mesh=Xh->mesh(),_geo=geoExportType);
    e->addRegions();

    for (int j=0;j<nTimeStep;j++)
    {
        if ( myDb.worldComm().isMasterRank() )
           std::cout << "Processing step : " << j << "/" << nTimeStep << "\n";

        myDb.load( timeSetIndex[j],fieldPod,uj );
        auto Qj = Idiv(uj) ;

        myDb.save( timeSetDb[j],"Qcriterion", Qj );

   
         // Writing result (ensight as default value)
         e->step( j )-> add( "Qcriterion", Qj );
         e->save();
    }

    return 0;
}
