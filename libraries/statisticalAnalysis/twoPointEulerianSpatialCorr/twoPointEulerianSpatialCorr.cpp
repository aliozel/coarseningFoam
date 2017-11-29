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
    int& nCell, 
    double& radii,
    double*& meanUf,
    double*& meanUs,    
    double**& cellPos, 
    double**& Uf,
    double**& Us,
    double*& voidfraction  
);

void calc_Rii
(
    int& nBin,
    double& nRp,
    int& nCell,    
    double& radii,
    double*& meanUf,
    double*& meanUs,              
    double**& cellPos, 
    double**& Uf,
    double**& Us,
    double*& voidfraction,    
    double**& Riif,
    double**& Riis,
    double*& rdf,        
    double*& nAve,
    double& aveAlpp
);

void write_data
(
    string& filename,
    int& nBin,
    double& nRp,
    double& radii,
    double**& Riif,
    double**& Riis,
    double*& rdf,
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
    	nRp = atof(argv[3]);
    else
    	nRp = 3;  
    //- Domain ave solid
    double aveAlpp; 
    if( argc > 4 )
        aveAlpp = atof(argv[4]);
    else
        aveAlpp = 0.1;

    //- Particle radius
    double radii; 
    if( argc > 5 )
        radii = atof(argv[5]);
    else
        radii = 37.5e-06;

    //- Read data	
    int nCell;
    double* meanUf;
    double* meanUs;
    double** cellPos;
    double** Uf;
    double** Us;      
    double* voidfraction;   
    read_data(filename,nCell,radii,meanUf,meanUs,cellPos,Uf,Us,voidfraction);
    
    cout << "Particle radius = " << radii << " nRp = " << nRp << " Max. distance = " << radii*nRp << " " 
         << "Number of bins = " << nBin << " Mean solid vol. frac. = " << aveAlpp << endl;	
 
    //- Calculate two-point spatial correlation
    double** Riif;
    double** Riis;
    double* rdf;
    double* nAve;
    calc_Rii(nBin,nRp,nCell,radii,meanUf,meanUs,cellPos,Uf,Us,voidfraction,Riif,Riis,rdf,nAve,aveAlpp);
    
    //- Write into the file
    write_data(filename,nBin,nRp,radii,Riif,Riif,rdf,nAve);
    
    return 0;
}


void read_data
(
    string& filename,
    int& nCell, 
    double& radii,
    double*& meanUf,
    double*& meanUs,    
    double**& cellPos, 
    double**& Uf,
    double**& Us,
    double*& voidfraction
)
{
    cout<< "Reading particle variables" << endl;
    string HH = string(filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ifstream inputPtr(charFilename);
    cout << "Opening file: " << HH.c_str() << endl;

    // Read data
    inputPtr >> nCell;
    cout << "Number of particles = " << nCell << endl;
    
    //- Define particle vels & pos
    int i;
    double posx, posy, posz; 
    double Ufx, Ufy, Ufz; 
    double Usx, Usy, Usz; 
    double voidf; 
      
    //- Init double double pointer
    cellPos = new double*[nCell];
    Uf = new double*[nCell];
    Us = new double*[nCell];
    voidfraction =  new double[nCell];
    
    for(int i=0; i < nCell; i++)
    {
       voidfraction[i] = 0.;
       cellPos[i] = new double[3];
       Uf[i] = new double[3];
       Us[i] = new double[3];
    }
    meanUf = new double[3];
    meanUs = new double[3];
    
    meanUf[0] = 0.; meanUf[1] = 0.; meanUf[2] = 0.;		
    meanUs[0] = 0.; meanUs[1] = 0.; meanUs[2] = 0.;		

    double meanAlpp;
    meanAlpp = 0.;  

    for(int index = 0; index < nCell; ++index)
    {			     
	 inputPtr >> posx 
		  >> posy 
		  >> posz 
		  >> Ufx   
		  >> Ufy
		  >> Ufz 
		  >> Usx   
		  >> Usy
		  >> Usz 
		  >> voidf;

	 /*
	 //- Testing create random numbers
	 posx = rand() % 50 + 1;
	 posy = rand() % 50 + 1;
	 posz = rand() % 50 + 1;
	 Ufx = rand() % 50 + 1;  
	 Ufy = rand() % 50 + 1;
	 Ufz = rand() % 50 + 1;
	 radii = 1.;		 
	 //cout << "Ufx = " << Ufx << endl; 		 
	 */

	 cellPos[index][0] = posx;	    
	 cellPos[index][1] = posy;
	 cellPos[index][2] = posz;
	 
	 Uf[index][0] = Ufx;	    
	 Uf[index][1] = Ufy;
	 Uf[index][2] = Ufz;
	 
	 Us[index][0] = Usx;	    
	 Us[index][1] = Usy;
	 Us[index][2] = Usz;
	 
	 voidfraction[index] = voidf;

	 meanUf[0] += voidf * Ufx;
	 meanUf[1] += voidf * Ufy;
	 meanUf[2] += voidf * Ufz;

	 meanUs[0] += ( 1. - voidf ) * Usx;
	 meanUs[1] += ( 1. - voidf ) * Usy;
	 meanUs[2] += ( 1. - voidf ) * Usz;
	 
	 meanAlpp += ( 1. - voidf );

     }			

     meanUf[0] /= ( 1. - meanAlpp );
     meanUf[1] /= ( 1. - meanAlpp );
     meanUf[2] /= ( 1. - meanAlpp );

     meanUs[0] /= meanAlpp;
     meanUs[1] /= meanAlpp;
     meanUs[2] /= meanAlpp;
          
     cout << "Mean solid volume fraction = " << meanAlpp/nCell << endl;
     cout << "Gas mean velocities = " 
          << meanUf[0] << "," << meanUf[1] << "," << meanUf[2] << endl;
     cout << "Solid mean velocities = " 
          << meanUs[0] << "," << meanUs[1] << "," << meanUs[2] << endl;     
     
     //-Close the file
     inputPtr.close();    
}

void calc_Rii
(
    int& nBin,
    double& nRp,
    int& nCell,    
    double& radii,
    double*& meanUf,
    double*& meanUs,              
    double**& cellPos, 
    double**& Uf,
    double**& Us,
    double*& voidfraction,    
    double**& Riif,
    double**& Riis,
    double*& rdf,        
    double*& nAve,
    double& aveAlpp
)
{
    double delta;
    delta = nRp*radii/( static_cast< double >(nBin) );

    //- Init double double pointer for correlation
    Riif = new double*[nBin];
    Riis = new double*[nBin];
    for(int i=0; i < nBin; i++)
    {
       Riif[i] = new double[6];
       Riis[i] = new double[6];
    }    

    //double** nAve;   
    nAve = new double[nBin];
    rdf = new double[nBin];
       
    for(int i=0; i < nBin; i++)
    {
        nAve[i] = 0.;
        rdf[i] = 0.;
        Riif[i][0] = 0.;  Riif[i][1] = 0.;  Riif[i][2] = 0.;
        Riif[i][3] = 0.;  Riif[i][4] = 0.;  Riif[i][5] = 0.;
        Riis[i][0] = 0.;  Riis[i][1] = 0.;  Riis[i][2] = 0.;
        Riis[i][3] = 0.;  Riis[i][4] = 0.;  Riis[i][5] = 0.;
    } 
               
    int i,j, iBin;
    double varUfx,varUfy,varUfz; 
    double varUsx,varUsy,varUsz; 
    varUfx = 0; varUfy = 0; varUfz = 0;  
    varUsx = 0; varUsy = 0; varUsz = 0; 

    double varalpf;
    varalpf = 0.; 
        
    double varalpp;
    varalpp = 0.; 

    double sumalpp;
    sumalpp = 0.;
    
    //- Ann KDtree, max number of particle in fixed radius
    int k = (static_cast<int>(0.64*aveAlpp*nRp*nRp*nRp));
    cout << "Number of particles searched in kDtree = " << k << endl; 

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
    //dataPts = annAllocPts(nCell, dim);
    nnIdx = new ANNidx[k];
    dists = new ANNdist[k];

    cout << "Allocated kdTree variables " << endl;

    cout << "Creating kdTree..." << endl;
    //kdTree = new ANNkd_tree(dataPts, nPts, dim);
    kdTree = new ANNkd_tree(cellPos, nCell, dim);
    
    // - Center particle j  
    for( j = 0; j < nCell; j++ )
    {
      double xtmp = cellPos[j][0];
      double ytmp = cellPos[j][1];
      double ztmp = cellPos[j][2];
          
      double vfxj = Uf[j][0] - meanUf[0];
      double vfyj = Uf[j][1] - meanUf[1];
      double vfzj = Uf[j][2] - meanUf[2];  
      
      varUfx += vfxj * vfxj;
      varUfy += vfyj * vfyj;
      varUfz += vfzj * vfzj;

      double vsxj = Us[j][0] - meanUs[0];
      double vsyj = Us[j][1] - meanUs[1];
      double vszj = Us[j][2] - meanUs[2];    

      varUsx += vsxj * vsxj;
      varUsy += vsyj * vsyj;
      varUsz += vszj * vszj;
      
      double alpfj = voidfraction[j];
      varalpf += alpfj * alpfj;
      double alppj = 1.-voidfraction[j];
      varalpp += alppj * alppj;
      sumalpp += alppj;
       
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
	   double delx = xtmp - cellPos[nnIdx[i]][0];
           double dely = ytmp - cellPos[nnIdx[i]][1];
           double delz = ztmp - cellPos[nnIdx[i]][2];

           double vfxi = Uf[nnIdx[i]][0] - meanUf[0];
           double vfyi = Uf[nnIdx[i]][1] - meanUf[1];
           double vfzi = Uf[nnIdx[i]][2] - meanUf[2];

           double vsxi = Us[nnIdx[i]][0] - meanUs[0];
           double vsyi = Us[nnIdx[i]][1] - meanUs[1];
           double vszi = Us[nnIdx[i]][2] - meanUs[2];
	   
	   double alpfi = voidfraction[i]; 
	   double alppi = 1.-voidfraction[i]; 

	   double dist = sqrt(delx*delx + dely*dely + delz*delz);
	   	
	   iBin = floor(dist/delta);
	   //cout << "iBin = " << iBin << endl;
	   //cout << " radii " << radii << endl;
	   if( iBin > 0 && iBin < nBin && dist > radii )
	   {
	      nAve[iBin]++;
	      	 	
	       Riif[iBin][0] += alpfi * alpfj * vfxi * vfxj;
	       Riif[iBin][1] += alpfi * alpfj * vfyi * vfyj;
	       Riif[iBin][2] += alpfi * alpfj * vfzi * vfzj;
	       Riif[iBin][3] += alpfi * alpfj * vfyi * vfxj;
	       Riif[iBin][4] += alpfi * alpfj * vfzi * vfyj;
	       Riif[iBin][5] += alpfi * alpfj * vfxi * vfzj;

	       Riis[iBin][0] += alppi * alppj * vsxi * vsxj;
	       Riis[iBin][1] += alppi * alppj * vsyi * vsyj;
	       Riis[iBin][2] += alppi * alppj * vszi * vszj;
	       Riis[iBin][3] += alppi * alppj * vsyi * vsxj;
	       Riis[iBin][4] += alppi * alppj * vszi * vsyj;
	       Riis[iBin][5] += alppi * alppj * vsxi * vszj;
	       
	       rdf[iBin] = alppi * alppj;
	       
	   }
	   	
 	 }
      }
    } 
    
    //- Normalize Rii
    for( i = 0; i < nBin; i++ )
    {
      if( nAve[i] > 0 ) Riif[i][0] = ( Riif[i][0]/nAve[i] ) / ( varalpf/nCell*varUfx/nCell );
      if( nAve[i] > 0 ) Riif[i][1] = ( Riif[i][1]/nAve[i] ) / ( varalpf/nCell*varUfy/nCell );
      if( nAve[i] > 0 ) Riif[i][2] = ( Riif[i][2]/nAve[i] ) / ( varalpf/nCell*varUfz/nCell );
      if( nAve[i] > 0 ) Riif[i][3] = ( Riif[i][3]/nAve[i] ) / ( varalpp/nCell*sqrt(varUfx/nCell)*sqrt(varUfy/nCell) );
      if( nAve[i] > 0 ) Riif[i][4] = ( Riif[i][4]/nAve[i] ) / ( varalpp/nCell*sqrt(varUfy/nCell)*sqrt(varUfz/nCell) );
      if( nAve[i] > 0 ) Riif[i][5] = ( Riif[i][5]/nAve[i] ) / ( varalpp/nCell*sqrt(varUfz/nCell)*sqrt(varUfx/nCell) );

      if( nAve[i] > 0 ) Riis[i][0] = ( Riis[i][0]/nAve[i] ) / ( varalpp/nCell*varUsx/nCell );
      if( nAve[i] > 0 ) Riis[i][1] = ( Riis[i][1]/nAve[i] ) / ( varalpp/nCell*varUsy/nCell );
      if( nAve[i] > 0 ) Riis[i][2] = ( Riis[i][2]/nAve[i] ) / ( varalpp/nCell*varUsz/nCell );
      if( nAve[i] > 0 ) Riis[i][3] = ( Riis[i][3]/nAve[i] ) / ( varalpp/nCell*sqrt(varUsx/nCell)*sqrt(varUsy/nCell) );
      if( nAve[i] > 0 ) Riis[i][4] = ( Riis[i][4]/nAve[i] ) / ( varalpp/nCell*sqrt(varUsy/nCell)*sqrt(varUsz/nCell) );
      if( nAve[i] > 0 ) Riis[i][5] = ( Riis[i][5]/nAve[i] ) / ( varalpp/nCell*sqrt(varUsz/nCell)*sqrt(varUsx/nCell) );

      if( nAve[i] > 0 ) rdf[i] = ( rdf[i]/nAve[i] ) / ( sumalpp/nCell*sumalpp/nCell );
      
    }
    
}

void write_data
(
    string& filename,
    int& nBin,
    double& nRp,
    double& radii,
    double**& Riif,
    double**& Riis,
    double*& rdf,
    double*& nAve
)
{
    double delta;
    delta = nRp*radii/( static_cast< double >(nBin) );
    
    string HH = string("Riif_"+filename);
    const char* charFilename;
    charFilename = HH.c_str();
    ofstream inputPtr(charFilename);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr << "# |d| " << "\t" << "Rfxx" << "\t" << "Rfyy" << "\t" << "Rfzz " 
                         << "\t" << "Rfxy" << "\t" << "Rfzy" << "\t" << "Rfxz " 
    			 << "\t" << "nAve" << endl;
        
    for( int i = 0; i < nBin; i++ )
    {
       inputPtr << (i+0.5)*delta << "\t" 
                <<     Riif[i][0] << "\t" 
		<<     Riif[i][1] << "\t"
		<<     Riif[i][2] << "\t"
                <<     Riif[i][3] << "\t" 
		<<     Riif[i][4] << "\t"
		<<     Riif[i][5] << "\t"
		<<     nAve[i]    << endl;
    }
    
    inputPtr.close();	
    
    HH = string("Riis_"+filename);
    const char* charFilename2 = HH.c_str();
    ofstream inputPtr2(charFilename2);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr2 << "# |d| " << "\t" << "Rsxx" << "\t" << "Rsyy" << "\t" << "Rszz " 
                          << "\t" << "Rsxy" << "\t" << "Rszy" << "\t" << "Rsxz " 
    			  << "\t" << "nAve" << endl;
        
    for( int i = 0; i < nBin; i++ )
    {
       inputPtr2 << (i+0.5)*delta << "\t" 
                 <<     Riis[i][0] << "\t" 
		 <<     Riis[i][1] << "\t"
		 <<     Riis[i][2] << "\t"
                 <<     Riis[i][3] << "\t" 
		 <<     Riis[i][4] << "\t"
		 <<     Riis[i][5] << "\t"
		 <<     nAve[i]    << endl;
    }    

    inputPtr.close();	

    HH = string("rdf_"+filename);
    const char* charFilename3 = HH.c_str();
    ofstream inputPtr3(charFilename3);
    cout << "Opening output file: " << HH.c_str() << endl;

    inputPtr3 << "# |d| " << "\t" << "Rdf" << "\t" << "nAve" << endl;
        
    for( int i = 0; i < nBin; i++ )
    {
       inputPtr3 << (i+0.5)*delta << "\t" 
                 <<     rdf[i]    << "\t" 
	 	 <<     nAve[i]   << endl;
    }    

    inputPtr3.close();	
         	
}
