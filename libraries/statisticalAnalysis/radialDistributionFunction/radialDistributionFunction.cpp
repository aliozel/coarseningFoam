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

#define PI 3.14152

// Neighboring algorithm 
#include "ANN.h"

using namespace std;

void read_data
(
    string& filename,
    int& nP, 
    double& radii,    
    double**& partPos    
);

void calc_Rdf
(
    int& nBin,
    double& rMax,
    int& nP,    
    double& radii,             
    double**& partPos, 
    //double*& rdf,        
    double**& rdf,        
    //double*& nAve,
    double**& nAve,
    double& Lx,
    double& Ly,
    double& Lz,
    double& delta,
    int& nIntPart         
);

void write_data
(
    string& filename,
    int& nBin,
    double& rMax,
    double& dp,
    //double*& rdf,
    double**& rdf,
    //double*& nAve,
    double**& nAve,
    int& nP,
    double& volDomain,
    double& delta,
    int& nIntPart                  
);

int main(  int argc, char *argv[] )
{
    //- Filename for input file
    string filename; 	
    if( argc > 1 )
    {
    	filename = argv[1];
    }
    else
    {
    	cout << "Define input file" << endl; 
	return 0;
    }	
    //- Define variables
    double dp;
    double radii;
    double Lx;
    double Ly;
    double Lz;
    double rMax;
    int nBin;    
    int nP;
    string filenameDatFile;
    
    cout<< "Reading input file" << endl;
    string HH = string(filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ifstream inputPtr(charFilename);
    cout << "Opening file: " << HH.c_str() << endl;
    
    //- 
    string dummyString;
    inputPtr >> dummyString >> filenameDatFile;
    inputPtr >> dummyString >> dp;		
    inputPtr >> dummyString >> nBin;
    inputPtr >> dummyString >> rMax;
    inputPtr >> dummyString >> Lx;
    inputPtr >> dummyString >> Ly;
    inputPtr >> dummyString >> Lz;

    //-Close the file
    inputPtr.close();   

    //- Read data	
    radii = dp/2.;
    double** partPos; 
    read_data(filenameDatFile,nP,radii,partPos);
    
    double volDomain;
    volDomain = Lx*Ly*Lz; 
    double aveAlpp = nP*(PI/6.*dp*dp*dp)/volDomain;  
    cout << "Particle radius = " 	<< radii         << " " 
         << "rMax = "             	<< rMax          << " " 
	 << "Max. distance = "   	<< radii*rMax    << " " 
         << "Number of bins = "  	<< nBin          << " "  
         << "Mean solid vol. frac. = "  << aveAlpp 	 << endl;	
	
 
    //- Variables
    //double* rdf;
    double** rdf;
    //double* nAve;
    double** nAve;
    
    //- 
    double delta;
    delta = rMax*radii/( static_cast< double >(nBin) );   

    //- Calculate g0
    int nIntPart;  
    calc_Rdf(nBin,rMax,nP,radii,partPos,rdf,nAve,Lx,Ly,Lz,delta,nIntPart);
    
    //- Write into the file
    write_data(filenameDatFile,nBin,rMax,radii,rdf,nAve,nP,volDomain,delta,nIntPart);

    return 0;
}


void read_data
(
    string& filename,
    int& nP, 
    double& radii,   
    double**& partPos
)
{
    cout<< "Reading particle variables" << endl;
    string HH = string(filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ifstream inputPtr(charFilename);
    cout << "Opening file: " << HH.c_str() << endl;

    // Read data
    inputPtr >> nP;
    cout << "Number of particles = " << nP << endl;

    //- Init double double pointer
    partPos = new double*[nP];
    for(int i=0; i < nP; i++)
    {
       partPos[i] = new double[3];
    }
    
    //- Define particle vels & pos
    int i;
    double posx, posy, posz; 
      
    for(int index = 0; index < nP; ++index)
    {			     
	 inputPtr >> posx 
		  >> posy 
		  >> posz;

	 partPos[index][0] = posx;	    
	 partPos[index][1] = posy;
	 partPos[index][2] = posz;   	 
     }			

     //-Close the file
     inputPtr.close();    
}

void calc_Rdf
(
    int& nBin,
    double& rMax,
    int& nP,    
    double& radii,             
    double**& partPos, 
    //double*& rdf,        
    double**& rdf,        
    //double*& nAve,
    double**& nAve,
    double& Lx,
    double& Ly,
    double& Lz,
    double& delta,
    int& nIntPart                 
)
{
    //nAve = new double[nBin];
    //rdf = new double[nBin];
       
    /*
    for(int i=0; i < nBin; i++)
    {
        nAve[i] = 0.;
        //rdf[i] = 0.;
    }
    */ 

    nAve = new double*[nBin];
    rdf = new double*[nBin];
    for(int i=0; i < nBin; i++)
    {
       rdf[i] = new double[4];
       nAve[i] = new double[4];       
    }
               
    int i,j;
    int iBin,iBinX,iBinY,iBinZ;
   
    //- Ann KDtree, max number of particle in fixed radius
    int k = (static_cast<int>(nP));
    cout << "Number of particles searched in kDtree = " << k << endl; 

    // Dimensions, exact OR approximate  
    int dim =3; double eps = 0;

    // Square of radius		
    double sqRad = ( rMax*radii )*( rMax*radii );

    cout << "Squared radius = " << sqRad << endl;

    // Data points
    ANNpointArray dataPts;
    // Query points
    ANNpoint queryPt;

    ANNidxArray		nnIdx;         	// 	near neighbour indices
    ANNdistArray 	dists;		//	near neighbour distances
    ANNkd_tree* 	kdTree;		//	search structure

    cout << "Created kdTree variables " << endl;

    // Allocate 
    queryPt = annAllocPt(dim);
    //dataPts = annAllocPts(nP, dim);
    nnIdx = new ANNidx[k];
    dists = new ANNdist[k];

    cout << "Allocated kdTree variables " << endl;

    cout << "Creating kdTree..." << endl;
    //kdTree = new ANNkd_tree(dataPts, nPts, dim);
    kdTree = new ANNkd_tree(partPos, nP, dim);
    
    // - Center particle j  
    nIntPart = 0;
    for( j = 0; j < nP; j++ )
    {
      double xtmp = partPos[j][0];
      double ytmp = partPos[j][1];
      double ztmp = partPos[j][2];
       
      //- Particles away from bcs
      if( (xtmp > rMax*radii) && (xtmp < Lx-rMax*radii)  )
      {
        if( (ytmp > rMax*radii) && (ytmp < Ly-rMax*radii)  )
        {
          if( (ztmp > rMax*radii) && (ztmp < Lz-rMax*radii)  )
          {
	      //- Inerior particles
	      nIntPart++;
	      
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

	      if ( j%(nP/10) == 0 ) cout << "Center particle = " << j << endl; 	 	 

	      //- KDtree algorithm particles i
	      for( i = 0; i < k; i++ )
	      {
		 //cout << "Center particle = " << j << " nnIdx[i] = " << nnIdx[i] << endl;
		 if( nnIdx[i] != ANN_NULL_IDX )
		 {
		   double delx = xtmp - partPos[nnIdx[i]][0];
        	   double dely = ytmp - partPos[nnIdx[i]][1];
        	   double delz = ztmp - partPos[nnIdx[i]][2]; 

		   double dist = sqrt(delx*delx + dely*dely + delz*delz);

		   iBin = floor(dist/delta);
		   iBinX = floor(abs(delx)/delta);
		   iBinY = floor(abs(dely)/delta);
		   iBinZ = floor(abs(delz)/delta);

		   if( iBin > 0 && iBin < nBin && dist >= 0. )
		   {
		      //nAve[iBin]++;
		      //rdf[iBin] += 1.; 
		      nAve[iBin][0]++;
		      rdf[iBin][0] += 1.; 
		   }

		   if( iBinX > 0 && iBinX < nBin && abs(dely) < radii && abs(delz) < radii )
		   {
		      nAve[iBinX][1]++;
		      rdf[iBinX][1] += 1.; 
		   }


		   if( iBinY > 0 && iBinY < nBin && abs(delx) < radii && abs(delz) < radii )
		   {
		      nAve[iBinY][2]++;
		      rdf[iBinY][2] += 1.; 
		   }


		   if( iBinZ > 0 && iBinZ < nBin && abs(delx) < radii && abs(dely) < radii )
		   {
		      nAve[iBinZ][3]++;
		      rdf[iBinZ][3] += 1.; 
		   }
		   		   
 		 }
	      }
	   }
	}
      } 	      
   }
   cout << "Interior particles = " << nIntPart << endl;       
}

void write_data
(
    string& filename,
    int& nBin,
    double& rMax,
    double& radii,
    //double*& rdf,
    double**& rdf,
    //double*& nAve,
    double**& nAve,
    int& nP,
    double& volDomain,
    double& delta,
    int& nIntPart                      
)
{
    string HH = string("rdf_"+filename);
    const char* charFilename3 = HH.c_str();
    ofstream inputPtr3(charFilename3);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr3 << "# |d| " << "\t" << "Rdf" << "\t" << "nAve" << endl;
    
    double numberDensity = nP/volDomain;
    //cout << "Number density = " << numberDensity << endl;   
    
    for( int i = 0; i < nBin; i++ )
    {
       double rOuter = (i+1)*delta;
       double rInner = i*delta;
       inputPtr3 << (rOuter + rInner) / 2. / (2.*radii)     	<< "\t" 
                 << rdf[i][0] / nIntPart / (4.0 / 3.0 * PI * (rOuter*rOuter*rOuter - rInner*rInner*rInner)) / numberDensity << "\t" 
   		 // Number of particles in shell/total number of particles/volume of shell/number density
	 	 << nAve[i][0]  	<< "\t" 
                 // Along x-dir
		 << rdf[i][1] / nIntPart / (4.0 / 3.0 * PI * (rOuter*rOuter*rOuter - rInner*rInner*rInner)) / numberDensity << "\t" 
	 	 << nAve[i][1]  	<< "\t" 		 
                 // Along y-dir
		 << rdf[i][2] / nIntPart / (4.0 / 3.0 * PI * (rOuter*rOuter*rOuter - rInner*rInner*rInner)) / numberDensity << "\t" 
	 	 << nAve[i][2]  	<< "\t" 
                 // Along z-dir
		 << rdf[i][3] / nIntPart / (4.0 / 3.0 * PI * (rOuter*rOuter*rOuter - rInner*rInner*rInner)) / numberDensity << "\t" 
	 	 << nAve[i][3]  		 		 
		 << endl;
    }    

    inputPtr3.close();	
         	
}
