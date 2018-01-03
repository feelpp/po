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

    double dt = (tf - ti) / nTimeStep ;

    // Time average fields computation
    /*if ( myDb.worldComm().isMasterRank() )
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
*/

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

    auto mem  = Environment::logMemoryUsage("memory usage after upate for use");
     std::cout << "[Correlation::updateForUse] resident memory before matrix product:     " << mem.memory_usage/1.e9  << "GBytes\n";


    for ( int i=0;i<nTimeStep;++i )
    {
	tic();
        myDb.load( timeSetIndex[i],fieldPod,ui );
	//ui = ui - Um ;
	toc("Loaded Ui");
        for ( int j=0;j<i;++j )
        {
	    tic();
            myDb.load( timeSetIndex[j],fieldPod,uj );
	    //uj = uj - Um ;
	    toc("Loaded Uj");

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

    auto mem2  = Environment::logMemoryUsage("memory usage after upate for use");
    std::cout << "[Correlation::updateForUse] resident memory after matrix product:     " << mem2.memory_usage/1.e9  << "GBytes\n";


    if ( myDb.worldComm().isMasterRank() )
    {
	    tic();
        //std::cout << "\npod=\n" << pod << "\n";
        std::ofstream file("pod.txt");
        if (file.is_open())
        {
            file << pod;
            file.close();
        }
        matrixPodFeel.printMatlab( "pod.m" );
	toc("Exporting POD matrix");
    }

    //Linearization
    std::cout<<"Linearization\n";
    auto R = pod.block(0,0,nTimeStep,nTimeStep-1);
    auto Rn = pod.block(0,nTimeStep-1,nTimeStep,1);

    //Eigen::ColPivHouseholderQR<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dec(R);
    //auto c = dec.solve(Rn);

    Eigen::ColPivHouseholderQR<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dec(R);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> c = dec.solve(Rn);
    std::cout<<"colPivQR"<<c.transpose()<<"\n";

    
    //Companion matrix
    std::cout<<"Companion matrix construction\n";
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> C(nTimeStep-1,nTimeStep-1);
    C.setZero();
    C.diagonal(-1).setOnes() ;
    C.block(0,nTimeStep-2,nTimeStep-1,1) = c ;

    
    // Solve Eigenvalues and eigenvectors
    std::cout<<"EigenSolver\n";
    Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
    eigenSolver.compute( C );

    auto eigenValues_ini = eigenSolver.eigenvalues();

    //double myeigen = real(eigenValues[k]);
    auto eigenVectors_ini = eigenSolver.eigenvectors();
    
    //ascending sort
    int nEigenValues = eigenValues_ini.size() ;
    std::vector<std::pair<double,int> > eigenValuesSortedWithInitialId;
    for (int k=0;k<nEigenValues;++k)
        eigenValuesSortedWithInitialId.push_back( std::make_pair(norm(eigenValues_ini[k]),k));
    std::sort(eigenValuesSortedWithInitialId.begin(),eigenValuesSortedWithInitialId.end(),
              [](std::pair<double,int> const& a, std::pair<double,int>const& b) { return b.first > a.first;});

    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> eigenValues(nEigenValues);
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> eigenVectors(nEigenValues,nEigenValues);

    for (int k1=0;k1<nEigenValues;++k1)
    {
         int k1Initial = eigenValuesSortedWithInitialId[k1].second;
         eigenValues(k1) = eigenValues_ini[k1Initial];
         for (int k2=0;k2<nEigenValues;++k2)
             eigenVectors(k2,k1) = eigenVectors_ini(k2,k1Initial);
    }
    std::cout<<"eigenValues sorted=\n"<<eigenValues<<"\n";

    auto mem3  = Environment::logMemoryUsage("memory usage after upate for use");
    std::cout << "[Correlation::updateForUse] resident memory after eigen:     " << mem3.memory_usage/1.e9  << "GBytes\n";


    //Projection
    if ( myDb.worldComm().isMasterRank() )
       std::cout<<"Modes projection and normalization"<<endl;

    tic();
    std::string geoExportType = "static";
    auto e = exporter( _mesh=space->mesh(),_geo=geoExportType);
    e->addRegions();

    toc("ini Exporter");

    tic();
    int nModes ;
    if ((nOutModes>nEigenValues) || (nOutModes == 0 ))
    {
         nModes = nEigenValues ;
    } else
    {
         nModes = nOutModes ;
    } 

    auto pi = M_PI ;
    for ( int i=0;i<nModes;++i )
    {

         auto PHI_real = space->element();
         auto PHI_imag = space->element();

	 double f ;
	 double signe ;
	 if (eigenValues(nEigenValues -1 -i ).imag()==0)
	 {
             signe = 1 ;
	     f = ( log(eigenValues(nEigenValues -1 -i )).imag() ) / ( 2*pi*dt ) ;
	 } else {
	     signe = eigenValues(nEigenValues -1 -i ).imag()/abs(eigenValues(nEigenValues -1 -i ).imag()) ;
	     f = ( log(eigenValues(nEigenValues -1 -i )).imag() ) / ( 2*pi*dt ) ;
         }
	 double G = ( log(eigenValues(nEigenValues -1 -i )).real() ) / dt ;

         // Mode i computation
         if ( myDb.worldComm().isMasterRank() )
             std::cout << "Mode " << i << "; eigenValues=" << eigenValues(nEigenValues -1 - i) <<"sign"<<signe<<" frequency=" << f << "Hz;  Growth rate=" << G <<"\n";

         for ( int j=0;j<nTimeStep;++j )
         {
             auto a = eigenVectors(j,nTimeStep -1 - i);
             myDb.load( timeSetIndex[j],fieldPod,uj );


             if ( j==0 )
             {
                PHI_real = (a.real()*uj);
                PHI_imag = (a.imag()*uj);
             }
             else
             {
                PHI_real = PHI_real + a.real()*uj;
                PHI_imag = PHI_imag + a.imag()*uj;
             }

         }

	 std::string real_path = "/data/scratch/atlas_abonnet/phi_real";
	 real_path += i;
	 std::string imag_path = "/data/scratch/atlas_abonnet/phi_imag";
	 imag_path += i;
	 
	 tic();
	 PHI_real.save( _path=real_path, _type="hdf5" );
	 PHI_imag.save( _path=imag_path, _type="hdf5" );
	 toc("HDF5 save");

	 tic();
	 tic();
         // Writing result (ensight as default value)
         e->step( i )-> add( "PHI_real", PHI_real );
	 e->step( i )-> add( "PHI_imag", PHI_imag );
	 toc("Ensight add");
	 tic();
	 e->save();
	 toc("Ensight save");
	 toc("Ensight export");
    }


    auto mem4  = Environment::logMemoryUsage("memory usage after upate for use");
    std::cout << "[Correlation::updateForUse] resident memory after 3D modes projection:     " << mem4.memory_usage/1.e9  << "GBytes\n";


    
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
	( "export-modes", po::value<int>()->default_value( 0 ), "Number of saved modes in ensight exports") 
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


    auto mem  = Environment::logMemoryUsage("memory usage after upate for use");
    std::cout << "[Mesh::updateForUse] resident memory before load mesh: " << mem.memory_usage/1.e9  << "GBytes\n";

    auto mesh = myDb.loadMesh( MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES );
    auto mem2  = Environment::logMemoryUsage("memory usage after upate fo    r use");
    std::cout << "[Mesh::updateForUse] resident memory before after mesh:     " << mem2.memory_usage/1.e9  << "GBytes\n";


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
