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

class Particle
{
    public:
    
    Particle()
    {}
    
    //Shallow copy
    Particle(const Particle& p)
    {
    	for(int i=0; i<3; i++)
	{
	   x[i] = p.x[i];
	}
        index = p.index;	
	dp = p.dp;
    }	 
    int index;
    double dp; 
    double x[3];
};

class GhostParticle : public Particle
{
    public:
    Particle* origp;
    GhostParticle(Particle* p) : Particle(*p)
    {
    	origp = p;
    }	
};

void read_data
(
    string& filenameDatFile,
    std::vector<Particle*> part, 
    std::vector<GhostParticle>& ghostPart 
);

void calc_rdf
(
    int& nBin,
    double& nDp,
    int& nP,    
    double& dp,            
    double**& partPos, 
    double*& rdf,        
    double*& nAve,
    double& aveAlpp   
);

void write_data
(
    string& filename,
    int& nBin,
    double& nDp,
    double& dp,
    double*& rdf,
    double*& nAve,
    double& aveAlpp,
    int& nP,
    double& volDomain     
);

int main(  int argc, char *argv[] )
{
    //- Filename for dump files
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
    double Lx;
    double Ly;
    double Lz;
    double nDp;
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
    inputPtr >> dummyString >> nDp;
    inputPtr >> dummyString >> Lx;
    inputPtr >> dummyString >> Ly;
    inputPtr >> dummyString >> Lz;

    //-Close the file
    inputPtr.close();   	        

    //- Particles
    std::vector<Particle*> part;
    //- Ghost particles
    std::vector<GhostParticle> ghostPart;
    
    //- Read data	
    read_data(filenameDatFile,part,ghostPart);
    
    double aveAlpp;
    double volP = 1./6.*PI*dp*dp*dp;
    double volDomain = Lx*Ly*Lz;
    aveAlpp=(part.size()*volP)/volDomain;	    
    cout << "Particle diameter = " << dp << " nDp = " << nDp << " Max. distance = " << dp*nDp << " " 
         << "Number of bins = " << nBin << " Mean solid vol. frac. = " << aveAlpp << endl;	
 
    //- Calculate rdf
    double* rdf;
    double* nAve; 
    
    //calc_rdf(nBin,nDp,nP,dp,partPos,rdf,nAve,aveAlpp);
    
    //- Write into the file
    //write_data(filenameDatFile,nBin,nDp,dp,rdf,nAve,aveAlpp,nP,volDomain);
    
    return 0;
}

void read_data
(
    string& filename,
    std::vector<Particle*> part, 
    std::vector<GhostParticle>& ghostPart 
)
{
    cout<< "Reading particle variables" << endl;
    string HH = string(filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ifstream inputPtr(charFilename);
    cout << "Opening file: " << HH.c_str() << endl;

    // Read data
    int nP;
    inputPtr >> nP;
    part.resize(nP);
    cout << "Number of particles = " << nP << endl;
    
    //part = new Particles[nP];
    cout << part.size() << endl;
    
    //- Define particle vels & pos
    int i;
    double posx, posy, posz; 
      
    for(int index = 0; index < nP; ++index)
    {			     
	 inputPtr >> posx 
		  >> posy 
		  >> posz;

	 cout << " index = " << index << endl;
	 part[index]->x[0] = posx;	    
	 part[index]->x[1] = posy;
	 part[index]->x[2] = posz;
	 
     }			
     
     //-Close the file
     inputPtr.close();    
}

void calc_rdf
(
    int& nBin,
    double& nDp,
    int& nP,    
    double& dp,            
    double**& partPos, 
    double*& rdf,        
    double*& nAve,
    double& aveAlpp
)
{
    double delta;
    delta = nDp*dp/( static_cast< double >(nBin) );   

    //double** nAve;   
    nAve = new double[nBin];
    rdf = new double[nBin];
       
    for(int i=0; i < nBin; i++)
    {
        nAve[i] = 0.;
        rdf[i] = 0.;
    } 
               
    int i,j, iBin;
   
    //- Ann KDtree, max number of particle in fixed radius
    int k = (static_cast<int>(nP));
    cout << "Number of particles searched in kDtree = " << k << endl; 

    // Dimensions, exact OR approximate  
    int dim =3; double eps = 0;

    // Square of radius		
    double sqRad = ( nDp*dp )*( nDp*dp );

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
    //dataPts = annAllocPts(nP, dim);
    nnIdx = new ANNidx[k];
    dists = new ANNdist[k];

    cout << "Allocated kdTree variables " << endl;

    cout << "Creating kdTree..." << endl;
    //kdTree = new ANNkd_tree(dataPts, nPts, dim);
    kdTree = new ANNkd_tree(partPos, nP, dim);
    
    // - Center particle j  
    for( j = 0; j < nP; j++ )
    {
      double xtmp = partPos[j][0];
      double ytmp = partPos[j][1];
      double ztmp = partPos[j][2];
       
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
	   //cout << "iBin = " << iBin << endl;
	   //cout << " dp " << dp << endl;
	   
	   
	   if( iBin > 0 && iBin < nBin && dist >= 0. )
	   {
	      nAve[iBin]++;
	      rdf[iBin] += 1.; 
	   }
	   	
 	 }
      }
    } 
    
}

void write_data
(
    string& filename,
    int& nBin,
    double& nDp,
    double& dp,
    double*& rdf,
    double*& nAve,
    double& aveAlpp,
    int& nP,
    double& volDomain      
)
{
    
    string HH = string("rdf_"+filename);
    const char* charFilename3 = HH.c_str();
    ofstream inputPtr3(charFilename3);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr3 << "# |d| " << "\t" << "Rdf" << "\t" << "nAve" << endl;

    double delta;
    delta = nDp*dp/( static_cast< double >(nBin) );
        
    for( int i = 0; i < nBin; i++ )
    {
       double r = (i+0.5)*delta;
       double voldr = 4.*PI*r*r*delta;
       double volS = 4./3.*PI*r*r*r;
       inputPtr3 <<   r      	<< "\t" 
                 <<  rdf[i]*voldr*nP/volDomain << "\t" 
	 	 << nAve[i]  	<< endl;
    }    

    inputPtr3.close();	
         	
}
