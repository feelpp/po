#include <feel/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feelvf/vf.hpp>

template<typename DatabaseType,typename SpaceType>
void
runPOD( DatabaseType & myDb, boost::shared_ptr<SpaceType> const& space, std::string const& fieldPod )
{
    using namespace Feel;
    bool useScalarProductL2 = soption(_name="scalar-product") == "L2";
    auto ui = space->element();
    auto uj = space->element();
    auto const& fieldsInfo = myDb.fieldsInfo();
    auto const& timeSetDb = myDb.timeSet();

    double ti = timeSetDb.front();
    double tf = timeSetDb.back();
    if ( Environment::vm().count("time-initial") )
        ti = doption(_name="time-initial");
    if ( Environment::vm().count("time-final") )
        tf = doption(_name="time-final");
    int nOutModes;
    if ( Environment::vm().count("export-modes") )
        nOutModes = Environment::vm()["export-modes"].as<int>();

    std::vector<double> timeSetIndex;
    for (int k=0; k<timeSetDb.size();++k)
    {
        if ( ( timeSetDb[k] >= (ti - 1e-9) ) &&
             ( timeSetDb[k] <= (tf + 1e-9 ) ) )
            timeSetIndex.push_back( k );
    }
    int nTimeStep = timeSetIndex.size();
    if ( nTimeStep == 0 )
    {
        std::cout << "exit pod because timeset is empty\n";
        return;
    }


    // Time average fields computation
    if ( myDb.worldComm().isMasterRank() )
       std::cout<<"Time average computation"<<endl;

    double invNstep = 1/((double) nTimeStep) ;
    auto Um = space->element();

    for ( int i=0;i<nTimeStep;++i )
    {
        if (i==0)
        {
             myDb.load( timeSetIndex[i],fieldPod,ui );
             Um = invNstep * ui ;
        }
        else
        {
             myDb.load( timeSetIndex[i],fieldPod,ui ) ;
	     Um = Um + invNstep * ui ;
        }
    }


    if ( myDb.worldComm().isMasterRank() )
        std::cout << "start build matrix pod of size : " << nTimeStep << "," << nTimeStep << "\n";

    boost::shared_ptr<MatrixSparse<double>> scalarProductOperator;
    bool useScalarProductOperator = boption(_name="scalar-product.use-operator");
    if ( useScalarProductL2 && useScalarProductOperator )
    {
        scalarProductOperator = backend()->newMatrix(_test=space,_trial=space);
        form2(_test=space,_trial=space,_matrix=scalarProductOperator ) =
            integrate(_range=elements(space->mesh()),_expr=inner(idt(ui),id(ui)) );
    }


    //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pod( timeSet.size(), timeSet.size() );
    MatrixEigenDense<double> matrixPodFeel( nTimeStep,nTimeStep );
    auto & pod = matrixPodFeel.mat();

    for ( int i=0;i<nTimeStep;++i )
    {
        myDb.load( timeSetIndex[i],fieldPod,ui );
	ui = ui - Um ;

        for ( int j=0;j<i;++j )
        {
            myDb.load( timeSetIndex[j],fieldPod,uj );
	    uj = uj - Um ;

            if ( useScalarProductL2 )
            {
                if ( useScalarProductOperator )
                    pod( i,j ) = scalarProductOperator->energy( ui,uj );
                else
                    pod( i,j ) = integrate( _range=elements(space->mesh()),
                                            _expr=inner(idv(ui),idv(uj) ) ).evaluate()(0,0);
            }
            else
            {
                pod( i,j ) = inner_product( ui, uj );
            }
            pod( j,i ) = pod( i,j );
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "("<<i<<","<<j<<")"<<std::flush;
        }
        if ( useScalarProductL2 )
        {
            if ( useScalarProductOperator )
                pod( i,i ) = scalarProductOperator->energy( ui,ui );
            else
                pod( i,i ) = integrate( _range=elements(space->mesh()),
                                        _expr=inner(idv(ui),idv(ui), mpl::int_<InnerProperties::IS_SAME/*|InnerProperties::SQRT*/>() ) ).evaluate()(0,0);
        }
        else
        {
            pod( i,i ) = inner_product( ui,ui );
        }
        if ( myDb.worldComm().isMasterRank() )
            std::cout << "("<<i<<","<<i<<")\n"<<std::flush;
    }

    if ( myDb.worldComm().isMasterRank() )
    {
        //std::cout << "\npod=\n" << pod << "\n";
        std::ofstream file("pod.txt");
        if (file.is_open())
        {
            file << pod;
            file.close();
        }
        matrixPodFeel.printMatlab( "pod.m" );
    }

    //Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
    //eigenSolver.compute( pod );

   Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
       eigenSolver.compute( pod );


    auto eigenValues = eigenSolver.eigenvalues();

    //double myeigen = real(eigenValues[k]);
    auto eigenVectors = eigenSolver.eigenvectors();
    if ( myDb.worldComm().isMasterRank() )
        std::cout << "eigenValues=\n" << eigenValues << "\n";






    //Projection
    if ( myDb.worldComm().isMasterRank() )
       std::cout<<"Modes projection and normalization"<<endl;

    tic();
    std::string geoExportType = "static";
    auto e = exporter( _mesh=space->mesh(),_geo=geoExportType);

    e->addRegions();

    int nModes ;
    if ((nOutModes>nTimeStep) || (nOutModes == 0 ))
    {
         nModes = nTimeStep ;
    } else
    {
         nModes = nOutModes ;
    } 

    for ( int i=0;i<nModes;++i )
    {
         
	 auto PHI = space->element();
         

	 // Mode i computation
	 if ( myDb.worldComm().isMasterRank() )
	     std::cout << "Mode " << i << "; eigenValues=" << eigenValues(nTimeStep -1 - i) << "\n";

	 for ( int j=0;j<nTimeStep;++j )
         {
             double a = eigenVectors(j,nTimeStep -1 - i);
             myDb.load( timeSetIndex[j],fieldPod,uj );
             
	     
             if ( j==0 )
             {
		PHI = a*uj;
	     } 
	     else
	     {
	        PHI = PHI + a*uj;
             }
	      
         }

	 // Norm 
	 double nPHI = normL2( _range=elements(space->mesh()), _expr=idv(PHI) );
	 double ni = 1/nPHI ;
	 auto PHI_norm = space->element() ;
	 PHI_norm = ni*PHI ;


         // Writing result (ensight as default value)
	 double time = (double) i ;
	 //if ( myDb.worldComm().isMasterRank() )
         //   std::cout << "time = " << time << endl;
         e->step( time )-> add( "PHI", PHI_norm );
	 e->save();
    }
 
    //e->save();
    toc("Exporter");

}








int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description exportdboptions( "export database options" );
	exportdboptions.add_options()
        ( "ifile", po::value<std::string>(), "input database file" )
        ( "field", po::value<std::string>(), "field loaded in database " )
        ( "scalar-product", po::value<std::string>()->default_value( "L2" ), "scalar product used in pod : L2,euclidian " )
        ( "scalar-product.use-operator", po::value<bool>()->default_value( true ), "build  scalar product operator  " )
        ( "time-initial", po::value<double>(), "initial time used for pod" )
        ( "time-final", po::value<double>(), "final time used for pod" )
	( "export-modes", po::value<int>()->default_value( 0 ), "Number of saved modes in ensight exports" )
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
    if ( Environment::vm().count("field") )
        fieldPod = Environment::vm()["field"].as<std::string>();

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


    if ( myDb.isScalarField( fieldPod ) )
    {
        auto Vh = Pch<1>( mesh );
        runPOD( myDb, Vh, fieldPod );
    }
    else if ( myDb.isVectorialField( fieldPod ) )
    {
        auto Vh = Pchv<1>( mesh );
        runPOD( myDb, Vh, fieldPod );
    }

    return 0;
}
