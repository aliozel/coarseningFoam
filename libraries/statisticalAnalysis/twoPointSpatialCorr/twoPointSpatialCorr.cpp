#include <fstream>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include "math.h"
#include <vector>

// Neighboring algorithm 
#include "ANN.h"

using namespace std;

void read_data
(
    string& filename,
    int& nPAll, 
    double& radii,
    double*& meanVel,
    double**& partAllPos, 
    double**& partAllVel
);

void calc_Rii
(
    int& nBin,
    double& nRp,
    int& nPAll,    
    double& radii,
    double*& meanVel,     
    double**& partAllPos, 
    double**& partAllVel,
    double**& Rii,
    double*& nAve,
    double& aveAlpp
);

void write_data
(
    string& filename,
    int& nBin,
    double& nRp,
    double& radii,
    double**& Rii,
    double*& nAve
);

int main(  int argc, char *argv[] )
{
    //- Filename for dump files
    string filename; 	
    if( argc > 1 )
    	filename = argv[1];
    else
    	filename = ((char*)"data.txt" );
	
    //-	Number of bins
    int nBin;
    if( argc > 2 )
    	nBin = atoi(argv[2]);
    else
    	nBin = 50;    

    //-	Maximum distance fron centering particle
    double nRp;
    if( argc > 3 )
    	nRp = atoi(argv[3]);
    else
    	nRp = 5;  
    //- Domain ave solid
    double aveAlpp; 
    if( argc > 4 )
        aveAlpp = atof(argv[4]);
    else
        aveAlpp = 1.;

    //- Read data	
    int nPAll;
    double radii;
    double* meanVel;
    double** partAllPos;
    double** partAllVel;    
    read_data(filename,nPAll,radii,meanVel,partAllPos,partAllVel);
    
    cout << "nRp = " << nRp << " Max. distance = " << radii*nRp << " Number of bins = " << nBin << " Mean solid vol. frac. " << aveAlpp << endl;	
 
    //- Calculate two-point spatial correlation
    double** Rii;
    double* nAve;
    calc_Rii(nBin,nRp,nPAll,radii,meanVel,partAllPos,partAllVel,Rii,nAve,aveAlpp);
    
    //- Write into the file
    write_data(filename,nBin,nRp,radii,Rii,nAve);
    
    return 0;
}


void read_data
(
    string& filename,
    int& nPAll, 
    double& radii,
    double*& meanVel,
    double**& partAllPos, 
    double**& partAllVel
)
{
    cout<< "Reading particle variables" << endl;
    string HH = string(filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ifstream inputPtr(charFilename);
    cout << "Opening file: " << HH.c_str() << endl;

    // Read data
    string just_read = " ";
    int just_read_int;
    double just_read_scalar;
    
    while( getline( inputPtr, just_read ) )
    {		    
	if( strcmp( just_read.c_str(), "ITEM: NUMBER OF ATOMS" ) == 0 )
	{		        
	    getline( inputPtr, just_read );			
	    nPAll = atoi( just_read.c_str() );
	    break;				
	}		    		    
    }

    cout << "Number of particles = " << nPAll << endl;
    
    //- Define particle vels & pos
    int i;
    double partPosx, partPosy, partPosz;
    double partVelx, partVely, partVelz;   
      
    //- Init double double pointer
    partAllPos = new double*[nPAll];
    partAllVel = new double*[nPAll];
    for(int i=0; i < nPAll; i++)
    {
       partAllPos[i] = new double[3];
       partAllVel[i] = new double[3];
    }
    meanVel = new double[3];
    
    double fieldDummy;
    meanVel[0] = 0.; meanVel[1] = 0.; meanVel[2] = 0.;		

    while( getline( inputPtr, just_read ) )
    {		    		    
	if( strncmp( just_read.c_str(), "ITEM: ATOMS id type", 19 ) == 0 )
	{		       
	    for(int index = 0; index < nPAll; ++index)
	    {			     
		 inputPtr >> fieldDummy
			  >> fieldDummy
			  >> partPosx 
			  >> partPosy 
			  >> partPosz 
			  >> partVelx   
			  >> partVely
			  >> partVelz 
			  >> fieldDummy
			  >> fieldDummy
			  >> fieldDummy
			  >> fieldDummy
			  >> fieldDummy
			  >> fieldDummy				      
			  >> radii;
		 
	 	 /*
		 //- Testing create random numbers
		 partPosx = rand() % 50 + 1;
		 partPosy = rand() % 50 + 1;
		 partPosz = rand() % 50 + 1;
		 partVelx = rand() % 50 + 1;  
		 partVely = rand() % 50 + 1;
		 partVelz = rand() % 50 + 1;
		 radii = 1.;		 
		 //cout << "partVelx = " << partVelx << endl; 		 
		 */
		 
		 partAllPos[index][0] = partPosx;	    
		 partAllPos[index][1] = partPosy;
	         partAllPos[index][2] = partPosz;
	 	 partAllVel[index][0] = partVelx;	    
		 partAllVel[index][1] = partVely;
	         partAllVel[index][2] = partVelz;
	 
		 meanVel[0] += partVelx;
		 meanVel[1] += partVely;
		 meanVel[2] += partVelz;
		 
	     }			
	     break;			
	 }
     } 

     meanVel[0] /= nPAll;
     meanVel[1] /= nPAll;
     meanVel[2] /= nPAll;
     
     cout << "Particle mean velocities = " 
          << meanVel[0] << "," << meanVel[1] << "," << meanVel[2] << endl;
     
     //-Close the file
     inputPtr.close();    
}

void calc_Rii
(
    int& nBin,
    double& nRp,
    int& nPAll,    
    double& radii,
    double*& meanVel,     
    double**& partAllPos, 
    double**& partAllVel,
    double**& Rii,
    double*& nAve,
    double& aveAlpp
)
{
    double delta;
    delta = nRp*radii/( static_cast< double >(nBin) );

    //- Init double double pointer for correlation
    Rii = new double*[nBin];
    for(int i=0; i < nBin; i++)
    {
       Rii[i] = new double[6];
    }    

    //double** nAve;   
    nAve = new double[nBin];
       
    for(int i=0; i < nBin; i++)
    {
        nAve[i] = 0.;
        Rii[i][0] = 0.;  Rii[i][1] = 0.;  Rii[i][2] = 0.;
        Rii[i][3] = 0.;  Rii[i][4] = 0.;  Rii[i][5] = 0.;
    } 
               
    int i,j, iBin;
    double varVelx,varVely,varVelz; 
    varVelx = 0; varVely = 0; varVelz = 0;  

    //- Ann KDtree, max number of particle in fixed radius
    int k = (static_cast<int>(0.64*aveAlpp*nRp*nRp*nRp));
    cout << "Number of particles searched in kDtree " << k << endl; 

    // Dimensions, exact OR approximate  
    int dim =3; double eps = 0;

    // Square of radius		
    double sqRad = ( nRp*radii )*( nRp*radii );

    cout << "Squared radius = " << sqRad << endl;

    // Data points
    ANNpointArray dataPts;
    // Query points
    ANNpoint queryPt;

    ANNidxArray	nnIdx;          	// 	near neighbour indices
    ANNdistArray 	dists;		//	near neighbour distances
    ANNkd_tree* 	kdTree;		//	search structure

    cout << "Created kdTree variables " << endl;

    // Allocate 
    queryPt = annAllocPt(dim);
    //dataPts = annAllocPts(nPAll, dim);
    nnIdx = new ANNidx[k];
    dists = new ANNdist[k];

    cout << "Allocated kdTree variables " << endl;

    cout << "Creating kdTree..." << endl;
    //kdTree = new ANNkd_tree(dataPts, nPts, dim);
    kdTree = new ANNkd_tree(partAllPos, nPAll, dim);
    
    // - Center particle j  
    for( j = 0; j < nPAll; j++ )
    {
      double xtmp = partAllPos[j][0];
      double ytmp = partAllPos[j][1];
      double ztmp = partAllPos[j][2];
    
      double vxj = partAllVel[j][0] - meanVel[0];
      double vyj = partAllVel[j][1] - meanVel[1];
      double vzj = partAllVel[j][2] - meanVel[2];    

      varVelx += vxj * vxj;
      varVely += vyj * vyj;
      varVelz += vzj * vzj;

      queryPt[0] = xtmp;
      queryPt[1] = ytmp;
      queryPt[2] = ztmp;

      //- KDtree algorithm 
      kdTree->annkFRSearch(
                              queryPt,				// query point					
                              sqRad,				// squared radius
                              k,                  		// number of the near neighbours to return
                              nnIdx,				// nearest neighbor array
                              dists,				// dist to near neighbours
                              eps			);
           
     if ( j%100000 == 0 ) cout << "Center particle = " << j << endl; 	 	 

      //- KDtree algorithm particles i
      for( i = 0; i < k; i++ )
      {
	 //cout << "Center particle = " << j << " nnIdx[i] = " << nnIdx[i] << endl;
	 if( nnIdx[i] != ANN_NULL_IDX )
	 {
	   double delx = xtmp - partAllPos[nnIdx[i]][0];
           double dely = ytmp - partAllPos[nnIdx[i]][1];
           double delz = ztmp - partAllPos[nnIdx[i]][2];

           double vxi = partAllVel[nnIdx[i]][0] - meanVel[0];
           double vyi = partAllVel[nnIdx[i]][1] - meanVel[1];
           double vzi = partAllVel[nnIdx[i]][2] - meanVel[2];

	   double dist = sqrt(delx*delx + dely*dely + delz*delz);
	   	
	   iBin = floor(dist/delta);
	   if( iBin > 0 && iBin < nBin && dist > radii )
	   {
	      nAve[iBin]++;	 	
	       Rii[iBin][0] += vxi * vxj;
	       Rii[iBin][1] += vyi * vyj;
	       Rii[iBin][2] += vzi * vzj;
	       Rii[iBin][3] += vyi * vxj;
	       Rii[iBin][4] += vzi * vyj;
	       Rii[iBin][5] += vxi * vzj;
	   }
	   	
 	 }
      }
    } 
    
    //- Normalize Rii
    for( i = 0; i < nBin; i++ )
    {
      if( nAve[i] > 0 ) Rii[i][0] = ( Rii[i][0]/nAve[i] ) / ( varVelx/nPAll );
      if( nAve[i] > 0 ) Rii[i][1] = ( Rii[i][1]/nAve[i] ) / ( varVely/nPAll );
      if( nAve[i] > 0 ) Rii[i][2] = ( Rii[i][2]/nAve[i] ) / ( varVelz/nPAll );
      if( nAve[i] > 0 ) Rii[i][3] = ( Rii[i][3]/nAve[i] ) / ( sqrt(varVelx/nPAll)*sqrt(varVely/nPAll) );
      if( nAve[i] > 0 ) Rii[i][4] = ( Rii[i][4]/nAve[i] ) / ( sqrt(varVely/nPAll)*sqrt(varVelz/nPAll) );
      if( nAve[i] > 0 ) Rii[i][5] = ( Rii[i][5]/nAve[i] ) / ( sqrt(varVelz/nPAll)*sqrt(varVelx/nPAll) );
    }
    
}

void write_data
(
    string& filename,
    int& nBin,
    double& nRp,
    double& radii,
    double**& Rii,
    double*& nAve
)
{
    string HH = string("Rii_"+filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ofstream inputPtr(charFilename);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr << "# |d| " << "\t" << "Rxx" << "\t" << "Ryy" << "\t" << "Rzz " 
                         << "\t" << "Rxy" << "\t" << "Rzy" << "\t" << "Rxz " 
    			 << "\t" << "nAvexx" << "\t" << "nAveyy" << "\t" << "nAvezz " 
    			 << "\t" << "nAvexy" << "\t" << "nAvezy" << "\t" << "nAvexz " << endl;

    double delta;
    delta = nRp*radii/( static_cast< double >(nBin) );
        
    for( int i = 0; i < nBin; i++ )
    {
       inputPtr << (i+0.5)*delta << "\t" 
                <<     Rii[i][0] << "\t" 
		<<     Rii[i][1] << "\t"
		<<     Rii[i][2] << "\t"
                <<     Rii[i][3] << "\t" 
		<<     Rii[i][4] << "\t"
		<<     Rii[i][5] << "\t"
		<<    nAve[i]    << endl;
    }
    
    inputPtr.close();	     	
}
