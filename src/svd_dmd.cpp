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

    auto pi = M_PI ;

    auto ui = space->element();
    auto uj = space->element();
    auto uj1 = space->element();

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

/*

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

    /* *****************************
    //Snapshot correlation matrix  *
    *******************************/

    MatrixEigenDense<double> matrixPodFeel( nTimeStep-1,nTimeStep-1 );
    auto & pod = matrixPodFeel.mat();

    MatrixEigenDense<double> matrixPodLFeel( nTimeStep-1,nTimeStep-1 );
    auto & podL = matrixPodLFeel.mat();

    for ( int i=0;i<nTimeStep-1;++i )
    {
        myDb.load( timeSetIndex[i],fieldPod,ui );
	//ui = ui - Um ;

        for ( int j=0;j<i;++j )
        {

            // Correlation U^T . U
            myDb.load( timeSetIndex[j],fieldPod,uj );


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

            // Correlation decalee U0^T . U1
            myDb.load( timeSetIndex[j+1],fieldPod,uj1 );

            if ( useScalarProductL2 )
            {
                if ( useScalarProductOperator )
                    podL( i,j ) = scalarProductOperator->energy( ui,uj1 );
                else
                    podL( i,j ) = integrate( _range=elements(space->mesh()),
                                            _expr=inner(idv(ui),idv(uj1) ) ).evaluate()(0,0);
            }
            else
            {
                podL( i,j ) = inner_product( ui, uj1 );
            }
            podL( j,i ) = podL( i,j );
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "("<<i<<","<<j<<")"<<std::flush;


        }

        myDb.load( timeSetIndex[i+1],fieldPod,uj1 );
        if ( useScalarProductL2 )
        {
            if ( useScalarProductOperator )
	    {
                pod( i,i ) = scalarProductOperator->energy( ui,ui );
                // Correlation decalee U0^T . U1
                podL( i,i ) = scalarProductOperator->energy( ui,uj1 );
	    } else
	    {
		pod( i,i ) = integrate( _range=elements(space->mesh()),
                                        _expr=inner(idv(ui),idv(ui), mpl::int_<InnerProperties::IS_SAME/*|InnerProperties::SQRT*/>() ) ).evaluate()(0,0);
                podL( i,i ) = integrate( _range=elements(space->mesh()),
                                        _expr=inner(idv(ui),idv(uj1), mpl::int_<InnerProperties::IS_SAME/*|InnerProperties::SQRT*/>() ) ).evaluate()(0,0);

	    }
        }
        else
        {
            pod( i,i ) = inner_product( ui,ui );
            podL( i,i ) = inner_product( ui,uj1 );

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


        //std::cout << "\npodL=\n" << podL << "\n";
        std::ofstream fileL("podL.txt");
        if (fileL.is_open())
        {
            fileL << podL;
            fileL.close();
        }
        matrixPodFeel.printMatlab( "podL.m" );

    }

   /* *************************
   // Snapshot POD computation*
   ***************************/
 
   Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > PODeigenSolver;
       PODeigenSolver.compute( pod );


    auto eigenValues = PODeigenSolver.eigenvalues();
    auto eigenVectors = PODeigenSolver.eigenvectors();
    if ( myDb.worldComm().isMasterRank() )
    {
        std::cout << "eigenValues=\n" << eigenValues << "\n";

        std::ofstream file1("EigenValuesPOD.txt");
        if (file1.is_open())
        {
            file1 << eigenValues;
            file1.close();
        }

        std::ofstream file2("EigenVectorsPOD.txt");
        if (file2.is_open())
        {
            file2 << eigenVectors;
            file2.close();
        }

     }


   //Singular value and Sigma
    if ( myDb.worldComm().isMasterRank() )
        std::cout << "Singular values computation \n";
    auto singValues = eigenValues.cwiseSqrt() ;

   /* *************************
   Dynamic modal decomposition*
   ****************************/

   //DMD matrix computation
   if ( myDb.worldComm().isMasterRank() )
        std::cout << "DMD : \n";

  auto nEigenValues = eigenValues.size(); 
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> lhs(nEigenValues,nEigenValues);
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rhs(nEigenValues,nEigenValues);
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Fdmd(nEigenValues,nEigenValues);

   if ( myDb.worldComm().isMasterRank() )
        std::cout << "side matrix computation \n";

  
   lhs = singValues.cwiseInverse().asDiagonal() * eigenVectors ;
   
   rhs = eigenVectors * singValues.cwiseInverse().asDiagonal(); 

   if ( myDb.worldComm().isMasterRank() )
        std::cout << "Modes computation \n";

   Fdmd = lhs * ( podL * rhs ) ;

   // DMD matrix decomposition

   if ( myDb.worldComm().isMasterRank() )
        std::cout<<"DMD EigenSolver\n";
   Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > DMDeigenSolver;
   DMDeigenSolver.compute( Fdmd );

   auto eigenValues_dmd = DMDeigenSolver.eigenvalues();
   auto eigenVectors_dmd = DMDeigenSolver.eigenvectors();
    
   if ( myDb.worldComm().isMasterRank() )
   {
        std::cout << "DMD eigenValues=\n" << eigenValues_dmd << "\n";
        std::ofstream file3("EigenValuesDMD.txt");
        if (file3.is_open())
        {
            file3 << eigenValues_dmd;
            file3.close();
        }

        std::ofstream file4("EigenVectorsDMD.txt");
        if (file4.is_open())
        {
            file4 << eigenVectors_dmd;
            file4.close();
        }
   }



   /*****************************
   //Sparse promoting resolution*
   *****************************/
/*
   // Temporal function
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ones(nEigenValues,nTimeStep);
   ones.setOnes();

   Eigen::Matrix<double,Eigen::Dynamic,1> ind(nTimeStep);
   ind.setLinSpaced(nTimeStep, 0.0, nTimeStep-1);

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Vand(nEigenValues,nTimeStep);
   Vand = eigenValues_dmd.asDiagonal() * ones ;
   Vand.array().pow((ones * ind.asDiagonal()).array());

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> P(nEigenValues,nEigenValues);
   P = (eigenVectors_dmd.adjoint() * eigenVectors_dmd).array() * ((Vand * Vand.adjoint()).conjugate().array());
 
   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> q(nEigenValues); 
   q = (Vand * eigenVectors * singValues.conjugate().asDiagonal() * eigenVectors_dmd).diagonal(0).conjugate();
   auto s = eigenValues.array().sum();


   //ADMM init

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> I(nEigenValues,nEigenValues);
   I.setIdentity();
   
   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> alpha(nEigenValues);
   alpha.setZero();

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> beta(nEigenValues);
   beta.setZero();

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> lambda(nEigenValues);
   lambda.setZero();

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> alphaN(nEigenValues);
   alphaN.setZero();

   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> betaN(nEigenValues);
   betaN.setZero();


   double resPrim = 1e6;
   double resDual = 1e6 ;


   //ADMM parameters
   int gamma = 1 ; 
   // gamma<<1 => maximise quality of least square approximation
   // gamma>>1 => minimize number of non zeros modes

   double rho = 1 ;

   double epsPrim = 1e-6;
   double epsDual = 1e-4 ;

   double epsAbs = 1e-6;
   double epsRel = 1e-4 ;




   // ADMM algo


   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Prho(nEigenValues,nEigenValues);
   Prho = P ; 
   Prho.diagonal() = P.diagonal() + 0.5*rho*I.diagonal() ;
   

   Eigen::LLT< Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> > Pchol_LLT(Prho);
   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> Plow(nEigenValues,nEigenValues); 
   Plow = Pchol_LLT.matrixL();
 
   auto Plow_star = Plow.adjoint(); 

   int ADMMstep = 1 ;
   while ((resPrim < epsPrim) && (resDual<epsDual) )
   {
	   ADMMstep++ ;
	   beta = betaN ; 

	   //alpha minimization step
	   auto u = beta -1/rho*lambda ;

           auto x = Plow.colPivHouseholderQr().solve(rho+0.5*u);
           alphaN = Plow.adjoint().colPivHouseholderQr().solve(x);




   }
*/
   /* *******************
   // Modes projection :*
   *********************/

   auto X = eigenVectors.transpose() * ( singValues.asDiagonal() * eigenVectors_dmd ) ;

   std::ofstream fileL("XEigVec.txt");
   if (fileL.is_open())
   {
      fileL << X;
      fileL.close();
   }


    if ( myDb.worldComm().isMasterRank() )
       std::cout<<"Modes projection and normalization"<<endl;

    tic();
    std::string geoExportType = "static";
    auto e = exporter( _mesh=space->mesh(),_geo=geoExportType);

    e->addRegions();



    int nModes ;
    if ((nOutModes>nEigenValues) || (nOutModes == 0 ))
    {
         nModes = nEigenValues ;
    } else
    {
         nModes = nOutModes ;
    } 

    for ( int i=0;i<nModes;++i )
    {
         
	 auto PHI_real = space->element();
         auto PHI_imag = space->element();

         double f = log(abs(eigenValues_dmd(nEigenValues -1 -i ).imag())) / ( 2*pi*dt ) ;
	 double G = log(abs(eigenValues_dmd(nEigenValues -1 -i ).real())) / dt ;

	 // Mode i computation
	 if ( myDb.worldComm().isMasterRank() )
             std::cout << "Mode " << i << "; eigenValues=" << eigenValues_dmd(nEigenValues -1 - i) << " frequency=" << f << "Hz;  Growth rate=" << G <<"\n";

	 for ( int j=0;j<nTimeStep;++j )
         {
             auto a = X(j,nTimeStep -1 - i);
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


         myDb.save(timeSetIndex[i],"PHI_real",PHI_real );
         myDb.save(timeSetIndex[i],"PHI_imag",PHI_imag );

         // Writing result (ensight as default value)
         e->step( i )-> add(fieldPod, ui ); 
         e->step( i )-> add( "PHI_real", PHI_real );
	 e->step( i )-> add( "PHI_imag", PHI_imag );
	 e->save();
    }
 
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

    auto mesh = myDb.loadMesh( MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED );


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
