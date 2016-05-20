#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include "acusim.h"
//#include "adb.h"
#include "h3dpublic_export.h"
#include <fstream>
#include <cmath>
#include <sstream>
//#include "nrutil.h"

//using Real = double;
typedef double Real;

#define SIGN( a, b ) a?a:b
#define SQR( a ) std::sqrt( a )

using namespace std;
Real** T;
int tlen;
double** a;
double* v;
double** q;
double** b;
double** c;
int Nt;
#define ROTATE( a, i, j, k, l )        \
    g = a[i][j];                       \
    h = a[k][l];                       \
    a[i][j] = g - s * ( h + g * tau ); \
    a[k][l] = h + s * ( g - h * tau );
#define SWAP( g, h )   \
    {                  \
        y = ( g );     \
        ( g ) = ( h ); \
        ( h ) = y;     \
    }
float pythag( float a, float b );

void jacobi( double** a, double* d, double** v, int* nrot );
void tred2( float** a, int n, float d[], float e[] );
void tqli( float d[], float e[], int n, float** z );
void ReportErrorMsg( H3DFileInfo* h3d_file, const char* error )
{
    FILE* errorFile = (FILE*)h3d_file->client_data1;
    if ( !errorFile )
    {
        errorFile = fopen( "export_error_messages", "a" );
        h3d_file->client_data1 = (void*)errorFile;
    }

    fprintf( errorFile, "%s\n", error );
    fflush( errorFile );
}

double scalaire( Real** u, Real** v, int taille, int taille1, int taille2 )
{
    double somme;
    for ( int i( 0 ); i < taille; i++ )
    {
        somme += u[taille1][i] * v[taille2][i];
    }
    return somme;
}

int main( int argc,
          char* argv[] )
{
    AdbHd adbHd;                                  /* an ADB handle		*/
    bool success;                                 /* new update			*/
    int event;                                    /* current event		*/
    int i, i2, i3, i0, itemp, ktemp, iphi, i2phi; /* a running index              */
    int j, k, k2, l, k3;
    int runId, Tstart, Tstop, Tmoy; /* run ID		*/
    int Ntpost;
    Real rTmp;       /* a temporary real		*/
    Real* logEvents; /* log events			*/
    char** logStrs;  /* list of log strings		*/
    int nOutSteps;
    int outStep;
    int* outSteps;
    Real outTime;
    int nNodes, nElmElems, nElmElemNodes;
    int nElms;
    int* elmCnn;
    Real* crd;
    int nOutVars;
    char outVarName, osfName;
    int outVarDim;
    Real* outValues;
    Real* outValues2;
    Real** VarMoyScal;
    bool rc;
    char elmName;
    const char errFilename[] = "export_error_messages";
    const char h3dFilename[] = "API_POD_Pfluct_1024_1050.h3d";
    Real* Pfluct;
    Real* Pfluct2;
    Real* Ufluct;
    Real* Vfluct;
    Real* Wfluct;
    Real* Ufluct2;
    Real* Vfluct2;
    Real* Wfluct2;
    Real* K;
    Real* phitemp;
    float** R;
    float** Rp;
    Real** phi;
    Real** alpha;
    Real** phi_p;
    Real** alpha_p;
    float** tempo;
    int* nrot;
    float* eigenval;
    float* eigenval2;
    float** eigenvect;
    float* eigenval_p;
    float* eigenval2_p;
    float** eigenvect_p;
    Real normephi;
    Real normephi_p;
    Real counttemp;
    //int*	usrIds ;		/* user numbers			*/
    int nOsfs;
    int nQTets;      /* No. qtet elements		*/
    int nQTrias;     /* No. QTria surface elements	*/
    int nQuads;      /* No. Quad surface elements	*/
    int nSrfs;       /* No. Surf sets		*/
    int nTets;       /* No. tet elements		*/
    int nTrias;      /* No. Tria surface elements	*/
    int nWedges;     /* No. wedge elements 		*/
    int nPyramids;   /* No. pyramid elements		*/
    int nBricks;     /* No. brick elements		*/
    int elmType;     /* element type			*/
    int nElemNodes;  /* No. element nodes		*/
    int* usrIds;     /* user numbers			*/
    float h3dCrd[3]; /* H3D Node position   		*/
    int* invMap;     /* inverse map			*/
    int* elmUsrIds;  /* user numbers			*/

    int nOutSurfs;  /* No. output surfaces		*/
    int srfType;    /* surface type (OSF, SBC, ...)	*/
    int srfId;      /* srf Id	       		*/
    int* srfUsrIds; /* user numbers			*/

    // clear old error msg file
    remove( errFilename );

    // open the export file
    H3DFileInfo* h3d_file = Hyper3DExportOpen( h3dFilename, H3D_SINGLEFILE, NULL, ReportErrorMsg );
    if ( !h3d_file ) return 1;

    // define pool names that will be used
    char creating_application[] = "jet pulseup";
    char file_creation_date[] = __DATE__;
    char original_data_file[] = "MESH.DIR";
    char file_comment[] = "export exercise created this file";

    /*---------------------------------------------------------------------------
   * Check the input
   *---------------------------------------------------------------------------

   if ( argc != 4 ) {
   printf( "Usage: %s problem directory run\n", argv[0] ) ;
   exit( 1 ) ;
   }*/

    // create Model block
    unsigned int model_count = 1;
    bool model_tabular = false;
    H3D_TRIBOOL model_adaptive = H3D_BOOL_FALSE;
    H3D_ID model_id = 1;
    rc = Hyper3DModelBegin( h3d_file, model_count );
    if ( !rc ) throw rc;
    rc = Hyper3DModelWrite( h3d_file, "PO_API_post_processing",
                            model_id, model_tabular, model_adaptive );
    if ( !rc ) throw rc;
    rc = Hyper3DModelEnd( h3d_file );
    if ( !rc ) throw rc;

    // when id != H3D_NULL_ID, SetModel must be called
    rc = Hyper3DSetModelToWrite( h3d_file, model_id, model_tabular );
    if ( !rc ) throw rc;

    // create File Info block
    rc = Hyper3DFileInfoBegin( h3d_file, creating_application,
                               file_creation_date );
    if ( !rc ) throw rc;

    rc = Hyper3DFileInfoAddModelFile( h3d_file, original_data_file );
    if ( !rc ) throw rc;
    rc = Hyper3DFileInfoAddResultFile( h3d_file, original_data_file );
    if ( !rc ) throw rc;
    rc = Hyper3DFileInfoAddComment( h3d_file, file_comment );
    if ( !rc ) throw rc;
    rc = Hyper3DFileInfoEnd( h3d_file );
    if ( !rc ) throw rc;

    // create Assemblies
    H3D_ID assm_poolname_id = H3D_NULL_ID;
    rc = Hyper3DAddString( h3d_file, H3D_DEFAULT_ASSMPOOL, &assm_poolname_id );
    if ( !rc ) throw rc;

    unsigned int assm_count = 2;
    rc = Hyper3DAssemblyBegin( h3d_file, assm_count, assm_poolname_id, assm_poolname_id );
    if ( !rc ) throw rc;

    H3D_ID assm_id = 1;
    H3D_ID model_as_parent = 0;
    //rc = Hyper3DAssemblyWrite(h3d_file, "Mesh", assm_id, model_as_parent);if( !rc ) throw rc;

    H3D_ID assm_parent = assm_id; // assm_parent = 1
                                  // assm_id++;                      // id = 2
    rc = Hyper3DAssemblyWrite( h3d_file, "2D", assm_id, model_as_parent );
    if ( !rc ) throw rc;
    assm_id++;
    //assm_parent++;
    rc = Hyper3DAssemblyWrite( h3d_file, "3D", assm_id, model_as_parent );
    if ( !rc ) throw rc;

    rc = Hyper3DAssemblyEnd( h3d_file );
    if ( !rc ) throw rc;

    // create Components
    H3D_ID comp_poolname_id = H3D_NULL_ID;
    rc = Hyper3DAddString( h3d_file, H3D_DEFAULT_COMPPOOL, &comp_poolname_id );
    if ( !rc ) throw rc;
    H3D_ID node_poolname_id = H3D_NULL_ID;
    rc = Hyper3DAddString( h3d_file, H3D_DEFAULT_NODEPOOL, &node_poolname_id );
    if ( !rc ) throw rc;

    unsigned int comp_count = 1;
    H3D_ID comp_id = 100;
    H3D_ID comp_parent_id = 2;

/*---------------------------------------------------------------------------
     * Open the data base and run
     *---------------------------------------------------------------------------
     */
    adbHd	= adbNew( "cda47d", "ACUSIM.DIR", "../", 0 ) ;             //Log name

    if ( adbHd == NULL )
      {
	printf( "Error opening the data base; %s", adbGetError() ) ;
	exit(1) ;
      }

    char* error =adbGetError();

    success = adbOpenRun (adbHd,1);       //which problem


    /* Get nodes*/

    success = adbGetInt0  ( adbHd, "nNodes",&nNodes);
    crd = (Real*) malloc(nNodes*3*sizeof(Real));

    success = adbGetReals0 ( adbHd,"crd", crd, 3 ,nNodes);
    unsigned int max_nodes = nNodes;
    rc = Hyper3DPositionBegin(h3d_file, max_nodes, node_poolname_id);           if( !rc ) throw rc;
    

 // pas de RENUMBER --> C'est la galère  !!!!!!!!!!!
    usrIds = (int*) malloc(nNodes*sizeof(int));
    success = adbGetInts0( adbHd,"usrIds", usrIds, nNodes , 1);
    for( i=0;i<max_nodes;i++)
      {
	float node[]={crd[3*i],crd[3*i+1],crd[3*i+2]};
	rc = Hyper3DPositionWrite(h3d_file,usrIds[i], node, H3D_NULL_ID, H3D_NULL_ID); if( !rc ) throw rc;
      }
    rc = Hyper3DPositionEnd(h3d_file);

     printf("\n ##############################################\n");
  printf("      Liste des surfaces      \n");

 // int nOsfs;         //  N surfaces
  int nOsfSrfNodes ;   //  N noeuds par elements
  int nOsfSrfs;     // N element ds la surface
  int nOsfNodes;    // N noeuds ds la surface
  int nOsfElemNodes;
  int* osfNodes;
  int* osfSrfCnn;
  int ki;
  H3D_ID elem_id = 1;
  H3D_ID elem2D_poolname_id = H3D_NULL_ID;
  rc = Hyper3DAddString(h3d_file, H3D_DEFAULT_ELEMPOOL2D,&elem2D_poolname_id);


  success = adbGetInt0(adbHd,"nOsfs",&nOsfs);
  comp_id=1;
  for(j=0;j<nOsfs;j++)    // Boucle sur les surfaces    VARIABLE J !!!!!!
  {
    success = adbGetStr1(adbHd,"osfName",&osfName,j);
    printf("  Surface %i : %s \n",j,&osfName);

    comp_parent_id = 1;
    rc = Hyper3DComponentBegin(h3d_file, comp_count, comp_poolname_id,
			       assm_poolname_id);                  if( !rc ) throw rc;
    rc = Hyper3DComponentWrite(h3d_file, &osfName, comp_id, node_poolname_id,
			       comp_parent_id);                    if( !rc ) throw rc;
    rc = Hyper3DComponentEnd(h3d_file);   if( !rc ) throw rc;


      success = adbGetInt1(adbHd,"nOsfSrfNodes",&nOsfSrfNodes,j);
//      printf("       Num de noeuds par elements =  %i\n",nOsfSrfNodes);
      
      success = adbGetInt1(adbHd,"nOsfSrfs",&nOsfSrfs,j);
//     printf("       Nombre d'elements dans la surface  %i \n",nOsfSrfs);
      
      success = adbGetInt1(adbHd,"nOsfNodes",&nOsfNodes,j);
//      printf("       Num de noeuds ds la surface =  %i\n",nOsfNodes);

      success = adbGetInt1(adbHd,"nOsfElemNodes",&nOsfElemNodes,j);
//      printf("       num of surface-elem nodes  %i \n",nOsfElemNodes);

      osfNodes=(int*)malloc(nOsfNodes*sizeof(int));
      success = adbGetInts1(adbHd,"osfNodes",osfNodes,j,nOsfNodes,1);  // alloc 320*1

      osfSrfCnn=(int*)malloc(nOsfElemNodes*nOsfSrfs*sizeof(int));
      success = adbGetInts1(adbHd,"osfSrfCnn",osfSrfCnn,j,nOsfSrfNodes,nOsfSrfs); // alloc 4x158

     int kt;
     if (nOsfSrfNodes==3)
     {
      rc = Hyper3DElementBegin(h3d_file, nOsfSrfs, elem2D_poolname_id,
				 H3D_ELEM_CONFIG_TRIA3,comp_id, comp_poolname_id,
				 node_poolname_id);
      if( !rc ) throw rc;
     }
     if (nOsfSrfNodes==4)
     {
      rc = Hyper3DElementBegin(h3d_file, nOsfSrfs, elem2D_poolname_id,
				 H3D_ELEM_CONFIG_QUAD4,comp_id, comp_poolname_id,
				 node_poolname_id);
      if( !rc ) throw rc;
     }
     for(kt=0;kt<nOsfSrfs;kt++)
      {
	//unsigned int elem_count = 1;
        if (nOsfSrfNodes==3)
        {
          unsigned int conns[]={osfSrfCnn[3*kt+0],osfSrfCnn[3*kt+1],osfSrfCnn[3*kt+2]};
          rc = Hyper3DElementWrite(h3d_file,elem_id, conns);
          if( !rc ) throw rc;
         // printf(" conn(%i) = %i %i %i %i \r" , i,elmCnn[4*i],elmCnn[4*i+1],elmCnn[4*i+2],elmCnn[4*i+3]);
         // fflush(stdout);
         }
        if (nOsfSrfNodes==4)
         {
           unsigned int conns[]={osfSrfCnn[4*kt+0],osfSrfCnn[4*kt+1],osfSrfCnn[4*kt+2],osfSrfCnn[4*kt+3]};
           rc = Hyper3DElementWrite(h3d_file,elem_id, conns);
           if( !rc ) throw rc;
           //printf(" conn(%i) = %i %i %i %i \n" , kt,osfSrfCnn[4*kt+0],osfSrfCnn[4*kt+1],osfSrfCnn[4*kt+2],osfSrfCnn[4*kt+3]);
       	  }
	elem_id++;
       }   //END FOR kt
      rc = Hyper3DElementEnd(h3d_file);


   //  }  // end if surface culot
   
           comp_id++;
  }   // end for j
  

  printf("\n ##############################################\n");
  printf("\n");





    /*Create elements*/
    /*  Le paramètre 0 correspond au cnn_id */

    success = adbGetInt0(adbHd,"nElms",&nElms);               //Le nombre de cnn  /!\ nElms = nCNN !!!!!!!!!!!!!
    printf("\n #############################\n");
    printf("  Informations du modele\n");
    printf(" #############################\n");
    printf ("  Nombre de domaines fluides : %i\n", nElms) ;

    success = adbGetStr1(adbHd,"elmName", &elmName,0);           // Le nom du cnn choisi


    //H3D_ID elem_id = 0;

    for (j=0;j<nElms;j++)
    {
    printf("  Chargement du fluide %s\n",&elmName);
    success = adbGetInt1(adbHd,"nElmElems",&nElmElems,j);          //Le nombre d'élems total 3D
    printf("  Le nombre d'elem est %i\n",nElmElems);
    success = adbGetInt1(adbHd,"nElmElemNodes",&nElmElemNodes,j);   // Le nombre de noeuds par elem
    printf("  Nombre de noeuds par element : %i  ( 4 = Tetra    6 = Prisme )\n",nElmElemNodes);

    //Chargement des connectivités
    elmCnn=(int*)malloc(nElmElems*nElmElemNodes*sizeof(int));
    success =adbGetInts1(adbHd,"elmCnn",elmCnn,j,nElmElemNodes,nElmElems);
    unsigned int max_elems=nElmElems;
    H3D_ID elem2D_poolname_id = H3D_NULL_ID;
    rc = Hyper3DAddString(h3d_file, H3D_DEFAULT_ELEMPOOL2D,
			  &elem2D_poolname_id);               if( !rc ) throw rc;
   // comp_id=1;
      //////////////////////
    comp_parent_id = 2;
    rc = Hyper3DComponentBegin(h3d_file, comp_count, comp_poolname_id,
			       assm_poolname_id);                  if( !rc ) throw rc;
    rc = Hyper3DComponentWrite(h3d_file, &elmName, comp_id, node_poolname_id,
			       comp_parent_id);                    if( !rc ) throw rc;
    rc = Hyper3DComponentEnd(h3d_file);   if( !rc ) throw rc;
     ///////////////////////////// 
       
       if(nElmElemNodes==4)   // test de reconnaissance penta ou tetra
        {
        rc = Hyper3DElementBegin(h3d_file, max_elems, elem2D_poolname_id,
				 H3D_ELEM_CONFIG_TETRA4,comp_id, comp_poolname_id,
				 node_poolname_id);                                  if( !rc ) throw rc;
        }
       else
        {
         rc = Hyper3DElementBegin(h3d_file, max_elems, elem2D_poolname_id,
				 H3D_ELEM_CONFIG_PENTA6,comp_id, comp_poolname_id,
				 node_poolname_id);                                  if( !rc ) throw rc;
        }  // end if

    for(i=0;i<max_elems;i++)
      {
	//unsigned int elem_count = 1;
        if(nElmElemNodes==4)
        {
         unsigned int conn[]={elmCnn[4*i],elmCnn[4*i+1],elmCnn[4*i+2],elmCnn[4*i+3]};
         rc = Hyper3DElementWrite(h3d_file,elem_id, conn);                         if( !rc ) throw rc;
         elem_id++;
        }
       else
        {
         unsigned int conn[]={elmCnn[6*i],elmCnn[6*i+1],elmCnn[6*i+2],elmCnn[6*i+3],elmCnn[6*i+4],elmCnn[6*i+5]};
       	 rc = Hyper3DElementWrite(h3d_file,elem_id, conn);     if( !rc ) throw rc;
         elem_id++;
	}
      }   //END FOR i
     
     rc = Hyper3DElementEnd(h3d_file);          if( !rc ) throw rc;
    }   //end j


    success= adbGetInt0(adbHd, "nOutSteps", &nOutSteps);
   /* for(i=0;i<nOutSteps;i++)
      {
	success= adbGetInt1(adbHd,"outStep",&outStep,i);
	printf("Time step %2d\n",outStep);
      }
      */
    outSteps =(int*)malloc(nOutSteps*sizeof(int));
    success =adbGetInts0(adbHd,"outSteps",outSteps,nOutSteps,1);
    printf("  Il y a %i echantillons disponibles \n",nOutSteps);
  printf(" #############################\n\n");
   /* // create Result defaults
    unsigned int  res_count = 1;
    unsigned int asys_count = 1;
    rc = Hyper3DResultBegin(h3d_file, res_count);                               if( !rc ) throw rc;
    rc = Hyper3DResultWrite(h3d_file, "Example1", H3D_SM_BISECT, asys_count);    if( !rc ) throw rc;
    rc = Hyper3DResultAddSystem(h3d_file, elem2D_poolname_id, H3D_DS_GLOBAL,
				H3D_POOL_ELEMENT);                              if( !rc ) throw rc;
    rc = Hyper3DResultEnd(h3d_file);                                            if( !rc ) throw rc;

*/
    success =adbGetInt0(adbHd,"nOutExtVars",&nOutVars);
    // create Subcase (Loadcase)
    int                max_sims = nOutSteps;
    unsigned int       dt_count = nOutVars+27 ;  //nombre de variables +14 pour TKE, Q criterion, Grad V, Vorticity magnitude, Laplacien P, Q/E
    unsigned int datatype_ids[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};

    //initialize
    H3D_ID subcase_id = 1;

    unsigned int      sub_count = 2;
    unsigned int   anim_grp_count = 1;
    unsigned int grp_datatype_ids = 1;
    unsigned int num_dts_per_grp  = 1;
    rc = Hyper3DSimSubcaseBegin(h3d_file, sub_count);                           if( !rc ) throw rc;
    rc = Hyper3DSimSubcaseWrite(h3d_file, "flow_solution", subcase_id, H3D_LINEAR_REAL,
				max_sims, dt_count, datatype_ids, H3D_NODAL_UNKNOWN); if( !rc ) throw rc;

    /*rc = Hyper3DSimSubcaseAnimationGroups(h3d_file, subcase_id, anim_grp_count,
					  &grp_datatype_ids, &num_dts_per_grp, &datatype_ids[1]);  if( !rc ) throw rc;     */

    rc = Hyper3DSimSubcaseEnd(h3d_file);                                        if( !rc ) throw rc;
    // create Simulations
    H3D_SIM_IDX sim_idx;
    rc = Hyper3DSimulationBegin(h3d_file, max_sims, subcase_id);                if( !rc ) throw rc;
    for( sim_idx=0; sim_idx < max_sims; sim_idx++ )
      {
	char name[32];
	snprintf(name, sizeof(name), "Mode = %ld", sim_idx);
	rc = Hyper3DSimulationWrite(h3d_file, sim_idx, name,(float)sim_idx);    if( !rc ) throw rc;
      }
    rc = Hyper3DSimulationEnd(h3d_file);                                        if( !rc ) throw rc;

    subcase_id++;

    // create Subcase2 (Loadcase)- creation du bloc pour les grandeurs moyennes
    rc = Hyper3DSimSubcaseBegin(h3d_file, sub_count);                           if( !rc ) throw rc;
    rc = Hyper3DSimSubcaseWrite(h3d_file, "Moyennes", subcase_id, H3D_LINEAR_REAL,
				1, dt_count, datatype_ids, H3D_NODAL_UNKNOWN); if( !rc ) throw rc;
    /*rc = Hyper3DSimSubcaseAnimationGroups(h3d_file, subcase_id, anim_grp_count,
					  &grp_datatype_ids, &num_dts_per_grp, &datatype_ids[1]);  if( !rc ) throw rc;   */

    rc = Hyper3DSimSubcaseEnd(h3d_file);                                        if( !rc ) throw rc;
    // create Simulationsbis
    rc = Hyper3DSimulationBegin(h3d_file, 1, subcase_id);                if( !rc ) throw rc;
    rc = Hyper3DSimulationWrite(h3d_file,0,"mean",0);              if( !rc ) throw rc;
    rc = Hyper3DSimulationEnd(h3d_file);                                        if( !rc ) throw rc;

    // create result data types
    char edata_type[32];
    H3D_ID                  dt_id = 1;
    unsigned int       pool_count = 0;
    unsigned int      layer_count = 0;
    const char**      layer_names = NULL;    
    bool              has_corners = false;
    H3D_TENSOR_TYPE   tensor_type = (H3D_TENSOR_TYPE)0;     // unused
    float                 poisson = 0.3f;                   // default & unused

    H3D_ID* layername_ids = new H3D_ID[1];
    layername_ids[0] = H3D_NULL_ID;

    rc = Hyper3DDatatypeBegin(h3d_file, dt_count);                              if( !rc ) throw rc;
   /*rc = Hyper3DDatatypeWrite(h3d_file, "Animation", dt_id, H3D_DS_NONE,
			      H3D_DS_UNKNOWN, pool_count);            if( !rc ) throw rc; */

    pool_count = 1;
    dt_id=1;
    for (i=0;i<nOutVars;i++)
    {
         success = adbGetStr1(adbHd,"outExtVarName",&outVarName,i);
         success  = adbGetInt1(adbHd,"outExtVarDim",&outVarDim,i);
         
         if (outVarDim==1)
            {
                 rc = Hyper3DDatatypeWrite(h3d_file, &outVarName, dt_id,
			      H3D_DS_SCALAR, H3D_DS_NODE, pool_count);        if( !rc ) throw rc;
                 rc = Hyper3DDatatypePools(h3d_file, dt_id, node_poolname_id,
			      layer_count, layername_ids, has_corners,
			      tensor_type, poisson);                          if( !rc ) throw rc;
             }
          else
              {
                 rc = Hyper3DDatatypeWrite(h3d_file, &outVarName, dt_id,
			      H3D_DS_VECTOR, H3D_DS_NODE, pool_count);        if( !rc ) throw rc;
                 rc = Hyper3DDatatypePools(h3d_file, dt_id, node_poolname_id,
			      layer_count, layername_ids, has_corners,
			      tensor_type, poisson);                          if( !rc ) throw rc;
             }
    dt_id++;
    }

    rc = Hyper3DDatatypeWrite(h3d_file, "Tke", dt_id,
			      H3D_DS_SCALAR, H3D_DS_NODE, pool_count);        if( !rc ) throw rc;
                 rc = Hyper3DDatatypePools(h3d_file, dt_id, node_poolname_id,
			      layer_count, layername_ids, has_corners,
			      tensor_type, poisson);                          if( !rc ) throw rc;
    dt_id++;

    rc = Hyper3DDatatypeWrite(h3d_file, "Phi", dt_id,
			      H3D_DS_VECTOR, H3D_DS_NODE, pool_count);        if( !rc ) throw rc;
    rc = Hyper3DDatatypePools(h3d_file, dt_id, node_poolname_id,
			      layer_count, layername_ids, has_corners,
			      tensor_type, poisson);                          if( !rc ) throw rc;
    dt_id++;
    rc = Hyper3DDatatypeWrite(h3d_file, "Phi_p", dt_id,
			      H3D_DS_SCALAR, H3D_DS_NODE, pool_count);        if( !rc ) throw rc;
    rc = Hyper3DDatatypePools(h3d_file, dt_id, node_poolname_id,
			      layer_count, layername_ids, has_corners,
			      tensor_type, poisson);                          if( !rc ) throw rc;

    rc = Hyper3DDatatypeEnd(h3d_file);                                          if( !rc ) throw rc;

    /*Get values*/
    unsigned int num_corners = 0;
    unsigned int   num_modes = 0;
    bool             complex = false;
    //subcase_id=1;

//###################################################################################
//###################################################################################
//###################################################################################

  Tstart=nOutSteps-256;       // nombre de modes calculés
  Ntpost = 50;                // Nombre de  modes enregistrés dans l'h3d
  Tstop=nOutSteps;
  Nt=Tstop-Tstart;
  Tmoy = nOutSteps-260;       // Nombre de pas de snapshots à prendre pour la moyenne

//###################################################################################
//###################################################################################
//###################################################################################

    K = new Real[nNodes];
    for(i=0;i<nNodes;i++) K[i]=0;
    R = new float*[Nt+1];
    for(i=0;i<Nt+1;i++)
    {
      R[i]=new float [Nt+1];
      for(j=0;j<Nt+1;j++)R[i][j]=0;
    }
    Rp = new float*[Nt+1];
    for(i=0;i<Nt+1;i++)
    {
      Rp[i]=new float [Nt+1];
      for(j=0;j<Nt+1;j++)Rp[i][j]=0;
    }

  VarMoyScal =new Real*[nOutVars];
  
  for(i=0;i<nOutVars;i++)
  {
    success  = adbGetInt1(adbHd,"outExtVarDim",&outVarDim,i);
    VarMoyScal[i]=new Real[nNodes*outVarDim] ;
  } //end for i
  printf(" #############################\n");
  printf("  Liste des variables :  \n");
  printf(" #############################\n");
  for(j=0;j<nOutVars;j++)    // Boucle sur les grandeurs    VARIABLE J !!!!!!
  {
    success = adbGetStr1(adbHd,"outExtVarName",&outVarName,j);
    printf("  %s - ID = %i \n",&outVarName,j);
  }   // end for j
  printf(" #############################\n");



//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
  printf("\n");
  printf("\n #############################\n");
  printf("  POST TRAITEMENT :  \n");
  printf(" #############################\n");
  for(j=0;j<nOutVars;j++)    // Boucle sur les grandeurs    VARIABLE J !!!!!!  fin ligne~956 !!!!!
  {
    //cout<<j<<endl;
    if ((j==11000) || (j==1))      // On choisi : 0-->vitesse 11-->vorticité 1-->pression
    {
    success = adbGetStr1(adbHd,"outExtVarName",&outVarName,j);
    printf("\n");
    printf("  Processing variable %s  -  dt_id = %i \n",&outVarName,j);
    success  = adbGetInt1(adbHd,"outExtVarDim",&outVarDim,j);
    //cout<<outVarDim<<endl;
    //VarMoyScal[j]=new Real[nNodes*outVarDim] ;

  //################################
  // Creation des GRANDEURS MOYENNES
  //################################

  for (itemp=Tmoy;itemp<Tstop;itemp++)             //boucle sur les pas de temps
   {
    outValues =new Real[nNodes*outVarDim];
   // printf(" paramètres : j=%i / itemp=%i / outVarDim=%i / nNodes=%i \n", j,itemp,&outVarDim,&nNodes);
    success = adbGetReals2(adbHd,"outExtValues",outValues,j,itemp,outVarDim,nNodes);
    counttemp = 100.00*(itemp-Tmoy+1)/(Tstop-Tmoy);
    printf("   Echantillon %i en %s importe : %4.2f%% effectue\r",itemp,&outVarName,counttemp );
    fflush(stdout);
    for(k=0;k<max_nodes;k++)
     {
     for (l=0;l<outVarDim;l++)
      {
       VarMoyScal[j][outVarDim*k+l] = VarMoyScal[j][outVarDim*k+l]+outValues[outVarDim*k+l]/(Tstop-Tmoy);
      }
     }
      delete [] outValues;
      //printf(" valeur moyenne en %s au point 15684 = %f \n",&outVarName,VarMoyScal[j][outVarDim*15684]);
   }          //fin pas de temps
      printf("\n");
   printf("   Moyenne en %s calculee\n", &outVarName);
  // FIN Creation des grandeurs moyennes
   } //END IF

  //##############################
  // Creation des modes de vitesse
  //##############################
  if(j==11000)     // 0 pour vitesse 11 pour vorticité !!!!!!
  {
   for (i=Tstart;i<Tstop;i++)
    {
     Ufluct = new Real[nNodes];
     Vfluct = new Real[nNodes];
     Wfluct = new Real[nNodes];
     //cout<<"Prem boucle time  "<<i<<endl;
     outValues =new Real[nNodes*outVarDim];
     success = adbGetReals2(adbHd,"outExtValues",outValues,j,i,outVarDim,nNodes);
     for(k=0;k<max_nodes;k++)
      {
       // cout<<outVarName<<endl;
       Ufluct[k]=outValues[outVarDim*k]  -VarMoyScal[j][outVarDim*k];
       Vfluct[k]=outValues[outVarDim*k+1]-VarMoyScal[j][outVarDim*k+1];
       Wfluct[k]=outValues[outVarDim*k+2]-VarMoyScal[j][outVarDim*k+2];
       K[k] =K[k]+ (Ufluct[k]*Ufluct[k]+Vfluct[k]*Vfluct[k]+Wfluct[k]*Wfluct[k])/(2*(Tstop-Tstart));
       }    //end node

       //for(i2=Tstart;i2<Tstop;i2++)
       for(i2=i;i2<Tstop;i2++)
        {
         Ufluct2 = new Real[nNodes];
         Vfluct2 = new Real[nNodes];
         Wfluct2 = new Real[nNodes];
         outValues2 =new Real[nNodes*outVarDim];
         success = adbGetReals2(adbHd,"outExtValues",outValues2,j,i2,outVarDim,nNodes);

         //*******************************************
         // Remplissage matrice correlation R(iphi,i2)
         //*******************************************
         for(k2=0;k2<max_nodes;k2++)
          {
             Ufluct2[k2]=outValues2[outVarDim*k2]  -VarMoyScal[j][outVarDim*k2];
             Vfluct2[k2]=outValues2[outVarDim*k2+1]-VarMoyScal[j][outVarDim*k2+1];
             Wfluct2[k2]=outValues2[outVarDim*k2+2]-VarMoyScal[j][outVarDim*k2+2];

             R[i-Tstart+1][i2-Tstart+1]= R[i-Tstart+1][i2-Tstart+1]+Ufluct[k2]*Ufluct2[k2]+Vfluct[k2]*Vfluct2[k2]+Wfluct[k2]*Wfluct2[k2];
          }   // end for k2

          delete [] Ufluct2;
          delete [] Vfluct2;
          delete [] Wfluct2;
          delete [] outValues2;
          //printf(" Correlation %s %i --> %i : R(%i|%i) = %f \n",&outVarName,i,i2,i,i2,R[i-Tstart+1][i2-Tstart+1]);
        }  //fin de la deuxième boucle en temps end for i2

      delete [] Ufluct;
      delete [] Vfluct;
      delete [] Wfluct;
      delete [] outValues;
      printf(" Echantillon %i correle : %f%% effectue\r",i,100.00*(i-Tstart+1)/Nt);
      //printf("-> %4.2f%%\r",100.00*(i-Tstart+1)/Nt);
      fflush(stdout);
     }      //end time  end for i
     printf("\n");
     //**********************
     // La matrice a été remplie moitie, comme elle est symetrique la fonction suivante symetrise
     //**********************
      for (i=1;i<Nt+1;i++)
      {
       for(i2=1;i2<i;i2++)
        {
          R[i][i2]  = R[i2][i];
        }
      }
      //****************
      // Sauvegarde de R
      //****************
      ofstream g2("Rrec.dat");
      {
       for(int s=0;s<Nt+1;s++)
        {
         //printf("Valeur propre %i = %f\n",s,eigenval[s]);
         for(int s2=0;s2<Nt+1;s2++)
         {
          g2<<R[s][s2]<<"  ";                     // NB : alpha[mode][temps]
         }
        g2<<endl;
         }
       }  // end ofstream

      //***************************************
      // Extraction valeurs et vecteurs propres
      //***************************************
      eigenval=new float[Nt+1];
      eigenval2=new float[Nt+1];

      cout<<" Calcul des valeurs propres"<<endl;
      tred2(R,Nt,eigenval,eigenval2);              // see numerical recipes p 480
      cout<<" Valeurs propres extraites"<<endl;
      tqli(eigenval,eigenval2,Nt,R);
      ofstream f("eigenvalues.dat");
      {
      for(int s=1;s<Nt+1;s++)
       {
        //printf("Valeur propre %i = %f\n",Nt-s,eigenval[s]);
        f<<eigenval[s]<<endl;
       }
      }
      
      //**************
      // CALCUL DE PHI
      //**************
      phi = new Real *[Nt];
      alpha = new Real *[Nt];
      for (i=0;i<Nt;i++)
      {
        alpha[i] =new Real[Nt];
        phi[i] = new Real[nNodes*outVarDim];
        //for (i2=0;i2<Nt;i2++) phi[i][i2]=0;
      }

      i0 = Nt-Ntpost;
      for(iphi=i0;iphi<Nt;iphi++)         //boucle sur les modes  : On prends les 25 premiers trop lourd sinon
      {
        printf(" Creation des modes : %4.2f%% effectue \r", 100.0*(iphi+1)/Nt);
        fflush(stdout);
        phitemp = new Real[nNodes*outVarDim];
        normephi=0;

        //printf("    Calcul de phi\n");
        
        for(i2phi=0;i2phi<Nt;i2phi++)               //boucle sur les pas de temps
        {
          outValues =new Real[nNodes*outVarDim];
          success = adbGetReals2(adbHd,"outExtValues",outValues,j,i2phi,outVarDim,nNodes);               //on récupère les données au temps i2phi
          for(k2=0;k2<nNodes;k2++)
          {
            for (l=0;l<outVarDim;l++) phitemp[outVarDim*k2+l]=phitemp[outVarDim*k2+l]+(outValues[outVarDim*k2+l]-VarMoyScal[j][outVarDim*k2+l])*R[i2phi+1][iphi+1];
           // normephi=normephi + phitemp[k2]*phitemp[k2];
          }
         delete [] outValues;
        }  // end time  end for i2phi

        //printf("    Calcul de la norme de phi\n");

        for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          normephi=normephi + phitemp[k2]*phitemp[k2];
         }
        normephi=sqrt(normephi);

        for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          phitemp[k2]=(1./normephi)*phitemp[k2];
         }

        //*****************
        // CALCUL DES ALPHA
        //*****************
        //printf("    Calcul de alpha\n");
        for(i3=0;i3<Nt;i3++)               //boucle sur les pas de temps
        {
          outValues =new Real[nNodes*outVarDim];
          success = adbGetReals2(adbHd,"outExtValues",outValues,j,i3,outVarDim,nNodes);               //on récupère les données au temps i3
          for(k3=0;k3<nNodes;k3++)
          {
            for (l=0;l<outVarDim;l++)  alpha[iphi][i3] = alpha[iphi][i3] + (outValues[outVarDim*k3+l]-VarMoyScal[j][outVarDim*k3+l])*phitemp[outVarDim*k3+l] ;
          }   //endfor k3
          delete [] outValues;
        }    //endfor i3

        //**********************
        // Phi gardes uniquement
        //**********************

        //printf("    Sauvegarde de phi et alpha\n");
        if (iphi >= i0)   // sauvegarde des modes contenus dans Ntpost
        {
          //printf(" iphi = %i  | i0 = %i\n", iphi,i0);
          for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          phi[iphi][k2]= phitemp[k2];
         }
         //printf("    Norme_PHI %i = %f  \n",Nt-iphi,normephi);
        } //end if
        delete [] phitemp;
      }    //endmodes END FOR iphi


      ofstream g("alpha.dat");
      {
      for(int s=i0;s<Nt;s++)
       {
        //printf("Valeur propre %i = %f\n",s,eigenval[s]);
        for(int s2=0;s2<Nt;s2++)
       {
         g<<alpha[s][s2]<<"  ";                     // NB : alpha[mode][temps]
       }
       g<<endl;
       }
      }  // end ofstream
         printf("\n");
   }      // endif




   //###############################
   // Creation des modes de Pression
   //###############################

   if (j==1) // C'est la PRESSION - creation des modes de pression
   {
    for (i=Tstart;i<Tstop;i++)  // boucle sur les snapshots
    {
      Pfluct = new Real[nNodes];
      outValues =new Real[nNodes*outVarDim];
      success = adbGetReals2(adbHd,"outExtValues",outValues,j,i,outVarDim,nNodes);
      for(k=0;k<max_nodes;k++)
      {
       Pfluct[k]=outValues[outVarDim*k]  -VarMoyScal[j][outVarDim*k];
       }    //end node  end for k
      for(i2=i;i2<Tstop;i2++)
        {
         Pfluct2 = new Real[nNodes];
         outValues2 =new Real[nNodes*outVarDim];
         success = adbGetReals2(adbHd,"outExtValues",outValues2,j,i2,outVarDim,nNodes);
         for(k2=0;k2<max_nodes;k2++)
          {
             Pfluct2[k2]=outValues2[outVarDim*k2]  -VarMoyScal[j][outVarDim*k2];
             Rp[i-Tstart+1][i2-Tstart+1]= Rp[i-Tstart+1][i2-Tstart+1] + Pfluct[k2]*Pfluct2[k2] ;
          }
          delete [] Pfluct2;
          delete [] outValues2;
          //printf("   Correlation pression %i --> %i : R(%i|%i) = %f \n",i,i2,i,i2,Rp[i-Tstart+1][i2-Tstart+1]);
         }  //end for i2
      delete [] Pfluct;
      delete [] outValues;
       printf("   Echantillon de pression %i correle : %4.2f%% effectue\r",i,100.00*(i-Tstart+1)/Nt);
    // printf("   Echantillon %i en %s importe : %4.2f%% effectue\r",itemp,&outVarName,counttemp );
    fflush(stdout);
     }      //end time end for i
     printf("\n");
     for (i=1;i<Nt+1;i++)
      {
        for(i2=1;i2<i;i2++)
        {
          Rp[i][i2]  = Rp[i2][i];
        }
      }
      
      //****************
      // Sauvegarde de R
      //****************
      ofstream g2("POD_R_pression.dat");
      {
       for(int s=0;s<Nt+1;s++)
        {
         //printf("   Valeur propre %i = %f\n",s,eigenval[s]);
         for(int s2=0;s2<Nt+1;s2++)
         {
          g2<<Rp[s][s2]<<"  ";                     // NB : alpha[mode][temps]
         }
        g2<<endl;
         }
       }  // end ofstream


      Real normephip;
      eigenval_p=new float[Nt+1];
      eigenval2_p=new float[Nt+1];

      cout<<" Calcul des valeurs propres de pression"<<endl;
      tred2(Rp,Nt,eigenval_p,eigenval2_p);              // see numerical recipes p 480
      cout<<" Valeurs propres extraites"<<endl;
      tqli(eigenval_p,eigenval2_p,Nt,Rp);

      ofstream f("POD_eigenvalues_P.dat");   // ecriture des valeurs propre dans un fichier
      {
      for(int s=1;s<Nt+1;s++)
       {
        //printf("    Valeur propre %i = %f\n",Nt-s,eigenval_p[s]);
        f<<eigenval_p[s]<<endl;
       }   // end for s
      }  // end ofstream

      phi_p = new Real*[Nt];
      alpha_p = new Real *[Nt];
      for (i=0;i<Nt;i++)
      {
        alpha_p[i] =new Real[Nt];
        phi_p[i] = new Real[nNodes*outVarDim];
       // for (i2=0;i2<Nt;i2++) phi_p[i][i2]=0;
      }
      i0 = Nt-Ntpost;
      //CALCUL DE PHI_p
      for(i=i0;i<Nt;i++)         //boucle sur les modes
      {
        printf("   Creation des modes : %4.2f%% effectue \r", 100.0*(i+1)/Nt);
        fflush(stdout);
        phitemp = new Real[nNodes*outVarDim];
        normephi_p=0;
        //printf("      Calcul de phi\n");
        for(i2=0;i2<Nt;i2++)               //boucle sur les pas de temps
        {
          outValues =new Real[nNodes*outVarDim];
          success = adbGetReals2(adbHd,"outExtValues",outValues,j,i2,outVarDim,nNodes);               //on récupère les données au temps i2
          for(k2=0;k2<nNodes;k2++)
          {
            phitemp[outVarDim*k2]=phitemp[outVarDim*k2]+(outValues[outVarDim*k2]-VarMoyScal[j][outVarDim*k2])*Rp[i2+1][i+1];
          }    // end for k2
         delete [] outValues;
        }  // end time  end for i2
        //normephi_p=normephi_p+scalaire(phi_p,phi_p,nNodes*outVarDim,i,i);
        //printf("      Calcul de la norme de phi\n");
        for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          normephi_p=normephi_p + phitemp[k2]*phitemp[k2];
         }

        normephi_p=sqrt(normephi_p);
        for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          phitemp[k2]=(1./normephi_p)*phitemp[k2];
         }
         //printf("    Norme_PHI_p %i = %f  \n",Nt-i,normephi_p);
        //printf("    Calcul de alpha\n");
        for(i3=0;i3<Nt;i3++)               //boucle sur les pas de temps   CALCUL DES ALPHA
        {
          outValues =new Real[nNodes*outVarDim];
          success = adbGetReals2(adbHd,"outExtValues",outValues,j,i3,outVarDim,nNodes);               //on récupère les données au temps i3
          for(k3=0;k3<nNodes;k3++)
          {
            alpha_p[i][i3] = alpha_p[i][i3] + (outValues[outVarDim*k3]-VarMoyScal[j][outVarDim*k3])*phitemp[outVarDim*k3] ;
          }   // end for k3
          delete [] outValues ;
        }    // end for i3
        //printf("    Sauvegarde de phi et alpha\n");
        if (i >= i0)   // sauvegarde des modes contenus dans Ntpost
        {
          //printf(" i = %i  | i0 = %i\n", i,i0);
          for(int k2=0;k2<nNodes*outVarDim;k2++)
         {
          phi_p[i][k2]= phitemp[k2];
         }
         //printf("    Norme_PHI_p %i = %f  \n",Nt-i,normephi_p);
        } //end if
      delete [] phitemp;
      }    //endmodes END FOR i
      printf("\n");
      ofstream h("POD_alpha_p.dat");
      {
      for(int s=i0;s<Nt;s++)
       {
        //printf("Valeur propre %i = %f\n",s,eigenval[s]);
        for(int s2=0;s2<Nt;s2++)
        {
          h<<alpha_p[s][s2]<<"  ";                     // NB : alpha[mode][temps]
        }
       h<<endl;
       }
      }  // end ofstream
   }  // END IF pression


  //#############################
  // create nodal Result Data set
  //#############################
  float p;
  float v[]={0.0f};
  if(outVarDim==3)
  {
  rc = Hyper3DDatasetBegin(h3d_file, max_nodes, 0,
		 subcase_id, H3D_DS_NODE, H3D_DS_VECTOR,
	         num_corners, num_modes,j+1, H3D_DS_NO_LAYER,
                  node_poolname_id, false); if( !rc ) throw rc;
  }
  else if (outVarDim==1)
   {
     rc = Hyper3DDatasetBegin(h3d_file, max_nodes, 0,
		 subcase_id, H3D_DS_NODE, H3D_DS_SCALAR,
	         num_corners, num_modes,j+1, H3D_DS_NO_LAYER,
                  node_poolname_id, complex); if( !rc ) throw rc;
    }
   else  cout<<"fatal error occurred...."<<endl;
   
   for (k=0;k<max_nodes;k++)
   {
    if (outVarDim!=1)
     {
      for (i=0;i<outVarDim;i++)
       {
	v[i]=VarMoyScal[j][outVarDim*k+i];
        }
        rc = Hyper3DDatasetWrite(h3d_file,usrIds[k],v);           if( !rc ) throw rc;
      }
    else
      {
          p=VarMoyScal[j][k];
          rc = Hyper3DDatasetWrite(h3d_file,usrIds[k],&p);           if( !rc ) throw rc;
      }
    }
   rc = Hyper3DDatasetEnd(h3d_file);   if( !rc ) throw rc;
   // end create nodal result data set

 } // end boucle sur les variables  J (Grangeurs)
  printf("\n #############################\n");
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************

   cout <<endl<<"Ecriture des résultats"<<endl<<endl;
   j =  nOutVars; //OK
   float p;
   rc = Hyper3DDatasetBegin(h3d_file, nNodes, 0,
		 subcase_id, H3D_DS_NODE, H3D_DS_SCALAR,
	         num_corners, num_modes,j+1, H3D_DS_NO_LAYER,
                  node_poolname_id, complex); if( !rc ) throw rc;

    for( k=0;k<max_nodes;k++)
	      {
                     p=K[k];
                     rc = Hyper3DDatasetWrite(h3d_file,usrIds[k],&p);           if( !rc ) throw rc;
              }
	    rc = Hyper3DDatasetEnd(h3d_file);   if( !rc ) throw rc;
        delete [] K;


   // On ecrit PHI dans le h3d
    /*float pod[]={0.0f};
    for(sim_idx=0;sim_idx<Ntpost;sim_idx++)
    {
     rc = Hyper3DDatasetBegin(h3d_file, max_nodes,sim_idx,
		 1, H3D_DS_NODE, H3D_DS_VECTOR,
	         num_corners, num_modes,j+2, H3D_DS_NO_LAYER,
                  node_poolname_id, complex); if( !rc ) throw rc;
     for (k=0;k<max_nodes;k++)
     {
        for (i=0;i<3;i++)
         {
       	   pod[i]=phi[Nt-1-sim_idx][3*k+i];
          }
          rc = Hyper3DDatasetWrite(h3d_file,usrIds[k],pod);           if( !rc ) throw rc;
     }
      rc = Hyper3DDatasetEnd(h3d_file);   if( !rc ) throw rc;
    }  // end for sim_idx
      */
      

    // On ecrit PHI_P dans le h3d
    float pod_p;
    for(sim_idx=0;sim_idx<Ntpost;sim_idx++)
    {
     rc = Hyper3DDatasetBegin(h3d_file, max_nodes,sim_idx,
		 1, H3D_DS_NODE, H3D_DS_SCALAR,
	         num_corners, num_modes,j+3, H3D_DS_NO_LAYER,
                  node_poolname_id, complex); if( !rc ) throw rc;
     for (k=0;k<max_nodes;k++)
     {
       pod_p=phi_p[Nt-1-sim_idx][k];
       rc = Hyper3DDatasetWrite(h3d_file,usrIds[k],&pod_p);           if( !rc ) throw rc;
     }
     rc = Hyper3DDatasetEnd(h3d_file);   if( !rc ) throw rc;
    }  // end for sim_idx


  //for (i=0;i<nOutVars;i++) delete VarMoyScal[i];
  //delete [] VarMoyScal;
  FILE* errorFile = (FILE*)h3d_file->client_data1;
  if( errorFile ) fclose(errorFile);
  bool rc2 = Hyper3DExportClose(h3d_file);
  if( !rc ) remove(h3dFilename);
    cout <<endl<<"  FIN du programme"<<endl<<endl;
  return !(rc && rc2) ? 1: 0;

} /* end of main() */

//############
//############
//
// END OF MAIN
//
//############
//############


//Subroutine pour le calcul des valeurs et vecteurs propres
float pythag(float a, float b)
//Computes (a2 + b2)1/2 without destructive underflow or overflow.
{
float absa,absb;
absa=fabs(a);
absb=fabs(b);
if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


void tred2(float **a, int n, float d[], float e[])
/*Householder reduction of a real, symmetric matrix a[1..n][1..n]. On output, a is replaced
by the orthogonal matrix Q effecting the transformation. d[1..n] returns the diagonal elements
of the tridiagonal matrix, and e[1..n] the off-diagonal elements, with e[1]=0. Several
statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
case a contains no useful information on output. Otherwise they are to be included.*/
{
int l,k,j,i;
float scale,hh,h,g,f;
for (i=n;i>=2;i--) {
l=i-1;
h=scale=0.0;
if (l > 1) {
for (k=1;k<=l;k++)
scale += fabs(a[i][k]);
if (scale == 0.0)
e[i]=a[i][l];
else {
for (k=1;k<=l;k++) {
a[i][k] /= scale;
h += a[i][k]*a[i][k];
}
f=a[i][l];
g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
e[i]=scale*g;
h -= f*g;
a[i][l]=f-g;
f=0.0;
for (j=1;j<=l;j++) {
/* Next statement can be omitted if eigenvectors not wanted */
a[j][i]=a[i][j]/h;
g=0.0;
for (k=1;k<=j;k++)
g += a[j][k]*a[i][k];
for (k=j+1;k<=l;k++)
g += a[k][j]*a[i][k];
e[j]=g/h;
f += e[j]*a[i][j];
}
hh=f/(h+h);
for (j=1;j<=l;j++) {
f=a[i][j];
e[j]=g=e[j]-hh*f;
for (k=1;k<=j;k++)
a[j][k] -= (f*e[k]+g*a[i][k]);
}
}
} else
e[i]=a[i][l];
d[i]=h;
}
/* Next statement can be omitted if eigenvectors not wanted */
d[1]=0.0;
e[1]=0.0;
/* Contents of this loop can be omitted if eigenvectors not
wanted except for statement d[i]=a[i][i]; */
for (i=1;i<=n;i++) {
l=i-1;
if (d[i]) {
for (j=1;j<=l;j++) {
g=0.0;
for (k=1;k<=l;k++)
g += a[i][k]*a[k][j];
for (k=1;k<=l;k++)
a[k][j] -= g*a[k][i];
}
}
d[i]=a[i][i];
a[i][i]=1.0;
for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
}
}

/*void tqli(float d[], float e[], int n, float **z)

/*QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric,
tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 §11.2. On
input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns
the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix,
with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix
that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].
{
float pythag(float a, float b);
//cout<<n<<endl;
int m,l,iter,i,k;
float s,r,p,g,f,dd,c,b;
for (i=2;i<=n;i++) e[i-1]=e[i];
e[n]=0.0;
for (l=1;l<=n;l++) {
iter=0;
do {
for (m=l;m<=n-1;m++) {

dd=fabs(d[m])+fabs(d[m+1]);
if ((float)(fabs(e[m])+dd) == dd) break;
}
if (m != l) {
//cout<<"iter="<<iter<<endl;
//if (iter++ == 30) nrerror("Too many iterations in tqli");
g=(d[l+1]-d[l])/(2.0*e[l]);
r=pythag(g,1.0);
g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
s=c=1.0;
p=0.0;
for (i=m-1;i>=l;i--) {
f=s*e[i];
b=c*e[i];
e[i+1]=(r=pythag(f,g));
if (r == 0.0) {
d[i+1] -= p;
e[m]=0.0;
break;
}
s=f/r;
c=g/r;
g=d[i+1]-p;
r=(d[i]-g)*s+2.0*c*b;
d[i+1]=g+(p=s*r);
g=c*r-b;
/* Next loop can be omitted if eigenvectors not wanted
for (k=1;k<=n;k++) {
f=z[k][i+1];
z[k][i+1]=s*z[k][i]+c*f;
z[k][i]=c*z[k][i]-s*f;
}
}
if (r == 0.0 && i >= l) continue;
d[l] -= p;
e[l]=g;
e[m]=0.0;
}
} while (m != l);
}
}*/


void tqli(float d[], float e[], int n, float **z)
{
int m, l, iter, i, k;
float s, r, p, g, f, dd, c, b;
for (i = 2; i <= n; i++)
    e[i-1] = e[i];
e[n] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m]) + fabs(d[m+1]);
          if (fabs(e[m]) + dd == dd) break;
          }
          if (m != l)
             {
             if (iter++ == 30) cout<<"too many iteration"<<endl; ;
             g = (d[l+1] - d[l]) / (2.0 * e[l]);
             r = sqrt((g * g) + 1.0);
             g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--)
                 {
                 f = s * e[i];
                 b = c * e[i];
                 if (fabs(f) >= fabs(g))
                    {
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else
                    {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1] - p;
                 r = (d[i] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                     {
                     f = z[k][i+1];
                     z[k][i+1] = s * z[k][i] + c * f;
                     z[k][i] = c * z[k][i] - s * f;
                     }
                 }
                 d[l] = d[l] - p;
                 e[l] = g;
                 e[m] = 0.0;
             }
          }  while (m != l);
      }
 }
