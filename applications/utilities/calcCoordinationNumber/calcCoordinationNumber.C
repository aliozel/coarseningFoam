/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "IOModel.H"

#include "meshSearch.H"

// Neighboring algorithm 
#include "ANN.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
		
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	

	Foam::timeSelector::addOptions();
	#   include "addRegionOption.H"
	Foam::argList::addBoolOption
	(
	"noWrite",
	"suppress writing results"
	);
	#include "addDictOption.H"

	#include "setRootCase.H"
	#include "createTime.H"
	Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

        // Post-processing dictionary
        #include "postProcessingDict.H"
			
	// Data points
	ANNpointArray dataPts;
	// Query points
	ANNpoint queryPt;

	ANNidxArray	nnIdx;          // 	near neighbour indices
	ANNdistArray 	dists;		//	near neighbour distances
	ANNkd_tree* 	kdTree;		//	search structure

	Pout << tab << "Created kdTree variables " << endl;

	// Allocate 
	//queryPt = annAllocPt(dim);
	//dataPts = annAllocPts(nPts, dim);
	nnIdx = new ANNidx[k];
	dists = new ANNdist[k];

	Pout << tab << "Allocated kdTree variables " << endl;

	//- Parallel computation
	int me;
	me = Pstream::myProcNo();
		
	//- Open cpus boundary file
	std::ifstream inputPtrCPUS;
	inputPtrCPUS.open("../DEM/in.cpusPost");

	int ncpus;
        if (me == 0) {	
	 inputPtrCPUS >> ncpus;
	 Info << "\nOpening ../DEM/in.cpusPost, nCpus = " << ncpus << endl;
        }
        Pstream::scatter(ncpus);
	
	scalarList coords(6*ncpus); 

	if (me == 0) {
	  for(int ii = 0; ii < 6*ncpus; ++ii)
	  {
              inputPtrCPUS >> coords[ii];
	  }
        }
			 				
	//- Global particle velocities
	int nPAll;	
	vector partAllPos(0,0,0);
	scalar partAllRadius(0);
	
	// For cases without in.CPUS file
	if( ncpus == 0 ) ncpus += 1 ;	
		
	//-Time loop
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		// Define coordination number variables
		#include "coordNumberEulVars.H"
		
		//- Local variables
		labelList nP(ncpus,0);	
		List<DynamicList<vector> > positions(ncpus);		
		List<DynamicList<scalar> > radii(ncpus);	

		labelList cellIDs;	
		
		// Read particle data
		#include "readParticleData.H"

		// Allocate searching algorithm variables
	        queryPt = annAllocPt(dim);
	        dataPts = annAllocPts(nP[me], dim);
						
		// Define a squared radius bound (mono-disperse case)
		sqRad = searchRadius * radii[me][0] * searchRadius * radii[me][0];
		
		Info << tab << "Squared-radius = " << sqRad << endl;

    		// Fill dataPoints for kd search
		for( int index = 0; index < nP[me]; index++ )
    		{
			dataPts[index][0] = positions[me][index][0];
			dataPts[index][1] = positions[me][index][1];
			dataPts[index][2] = positions[me][index][2];						
		}
		
		Pout << tab << "Creating kdTree..." << endl;
        	kdTree = new ANNkd_tree(dataPts, nP[me], dim);
		
		
		Info << tab << "Calculating coordination numbers..." << endl;
		    		
		for( int index = 0; index < nP[me]; index++ )
    		{		
			// Center particle cellID
			label cellID = cellIDs[index];
			
			// Number of particles in a cell which particles belong to
			numberOfParticlesInCell[cellID] += 1.;
			
			// Define query
			queryPt[0] = positions[me][index][0];
        		queryPt[1] = positions[me][index][1];
        		queryPt[2] = positions[me][index][2];
		
			// Fixed radius 
			kdTree->annkFRSearch(
                                		queryPt,			// query point					
                                		sqRad,				// squared radius
                                		k,                  		// number of the near neighbours to return
                                		nnIdx,				// nearest neighbor array
                                		dists,				// dist to near neighbours
                                		eps			);
			
			// Particle without collision
			bool withColl(false);
			
			int i = 1;			
			while( i < k )
		        {								
				//if( index < 1 ) Info  << " index = " << index <<  " k = " << i << " nnIdx = " << nnIdx[i] << nl;  
				//if( dists[i] <= 2.*radii[me][0] )
				
				scalar dist(1.e+64);
				if ( nnIdx[i] > -1 )
				       dist =  ( positions[me][index][0] - positions[me][nnIdx[i]][0] )
				              *( positions[me][index][0] - positions[me][nnIdx[i]][0] )
					     + ( positions[me][index][1] - positions[me][nnIdx[i]][1] )
				              *( positions[me][index][1] - positions[me][nnIdx[i]][1] )
					     + ( positions[me][index][2] - positions[me][nnIdx[i]][2] )
				              *( positions[me][index][2] - positions[me][nnIdx[i]][2] )	;
					      
				
				if( dist <= ( 2.*radii[me][0] * 2.*radii[me][0] )  )
				{
					if( index < 1 ) Info  << " colliding index = " << index <<  " k = " << i << " nnIdx = " << nnIdx[i] << endl;  

				 	coordNumberEul[cellID] += 1.;
					withColl = true;										
				}
				i++;
			} 			
			
			// Count particles without collisions
			if(!withColl) numberOfParticlesInCellWithoutColl[cellID] += 1.;
		}
		
		// Normalize coordination number computations and write them into files
		Info << tab << "Writing coordination number variables into the files " << endl;
		forAll(mesh.cells(),cellI)
		{
			coordNumberEul[cellI] /= numberOfParticlesInCell[cellI];
			if( (numberOfParticlesInCell[cellI] - numberOfParticlesInCellWithoutColl[cellI]) > 0 )
			{
				coordNumberEulWithoutZeroContact[cellI] /= ( numberOfParticlesInCell[cellI] - numberOfParticlesInCellWithoutColl[cellI] );  
			}else
			{
				coordNumberEulWithoutZeroContact[cellI] = 0.;
			}
		}
		numberOfParticlesInCell.write();
		numberOfParticlesInCellWithoutColl.write();
		coordNumberEul.write();
		coordNumberEulWithoutZeroContact.write();
		
			
		positions.clear();
		radii.clear();	
	
	}

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
