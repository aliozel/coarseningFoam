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

#include "particleStressBi.H"
#include "calcCollisionalForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	// Initiate output file
	void initOutput
	(	
		const int&          nParticleClass_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_		
	)
	{
	
	    	outputKinFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"kineticStress");
		*outputKinFile_ << "#Time  " << tab << " "
		                << "Vol    " <<        " "  
	                	<< "alpp   " << tab << " "
				<< "Kin_XX " << tab << " "
				<< "Kin_YY " << tab << " "
				<< "Kin_ZZ " << tab << " " 
				<< "Kin_XY " << tab << " " 
				<< "Kin_XZ " << tab << " " 
				<< "Kin_YZ " << tab << " "; 		

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_  << "Kin_XX["<<iPartClass<<"] " << tab << " "
						<< "Kin_YY["<<iPartClass<<"] " << tab << " "
						<< "Kin_ZZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XY["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_YZ["<<iPartClass<<"] " << tab << " "; 
			}
		}
		*outputKinFile_ << endl;
	
		if( calcCollision_ )
		{
    			outputCollFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"collisonalStress");
			*outputCollFile_ << "#Time  " << tab << " "
		                	 << "Vol    " << 	" "  
	                		 << "alpp   " << tab << " "			 
	                		 << "Coll_XX " << tab << " "
					 << "Coll_YY " << tab << " "
					 << "Coll_ZZ " << tab << " " 
					 << "Coll_XY " << tab << " " 
					 << "Coll_XZ " << tab << " " 
					 << "Coll_YZ " << tab << " "; 	
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_  << "Coll_XX["<<iPartClass<<"] " << tab << " "
							 << "Coll_YY["<<iPartClass<<"] " << tab << " "
							 << "Coll_ZZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XY["<<iPartClass<<"] " << tab << " " 
							 << "Coll_XZ["<<iPartClass<<"] " << tab << " " 
							 << "Coll_YZ["<<iPartClass<<"] " << tab << " "; 			
				}
			}
			*outputCollFile_ << endl;	
		}
	}
	
	// Write into output file
	void writeOutput
	(	
		const int&          nParticleClass_,
		const Time& 		   runTime_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		OFstream*&          outputCollFile_,
		const bool&             bidisperse_,
		const bool&          calcCollision_,
		SymmTensor<double>*       sigmaKin_,					
		SymmTensor<double>*	 sigmaColl_,
		const scalar&		 domainVol_,
		const scalar&	   alppInSubVolume_	)
	{
	
		// Write output
		Info << " " << endl;
		Pout << " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"kineticStress" << endl;

		int iPartClass = 0;
		*outputKinFile_ << runTime_.value() 		       << tab << " " 
		                << domainVol_ 			       << tab << " "
				<< alppInSubVolume_		       << tab << " "   	
				<< sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ;								

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_ << sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ; 
			}
		}
		*outputKinFile_ << endl;

		if( calcCollision_ )
		{
    			Pout<< " Writing particle kinetic stresses into the file " << mesh_.time().path()/outputRelativePath_/"collisionalStress" << endl;
			iPartClass = 0;
			*outputCollFile_ << runTime_.value() 		         << tab << " " 	
					 << domainVol_ 			         << tab << " "
					 << alppInSubVolume_		         << tab << " "  
					 << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
					 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
					 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;
			if( bidisperse_ ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        		{
					*outputCollFile_ << sigmaColl_[iPartClass][0]/domainVol_ << tab << " " 		
					  		 << sigmaColl_[iPartClass][1]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][2]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][3]/domainVol_ << tab << " " 
							 << sigmaColl_[iPartClass][4]/domainVol_ << tab << " " 		
							 << sigmaColl_[iPartClass][5]/domainVol_ << tab << " " ;			
				}
			}
			*outputCollFile_ << endl;	
		}	
	}

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

	#include "createFields.H"  

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);
	
	particleCloud.reAllocArraysPost();

	double **positions;
	double **velocities;
	double **omegas;
	double **radii;
	double **cellID;
	
	double **forces;
	double **types;
	
	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.dataExchangeM().allocateArray(omegas,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);
	particleCloud.dataExchangeM().allocateArray(types,0.,1);
				
	//- Read post-processing dicitonary
	#include "postProcessingDict.H"

	// Initiate output file
	OFstream* outputKinFile;
	OFstream* outputCollFile;

    	initOutput(	    nParticleClass,
	                              mesh,
			outputRelativePath,
			     outputKinFile,
			    outputCollFile,
			        bidisperse,
			     calcCollision	);
			     
			     

     
	// Define collision variables
	#include "initCollisionParameters.H" 
						
	if(calcCollision)
	{
		#include "collisionParameters.H"
	}
		
	// Neighboring list parameters				
	int k(0); 

	// Dimensions, exact OR approximate  
	int dim(0); double eps(0);

	// Number of points
	int nPts(0);

	// Data points
	ANNpointArray dataPts;
	// Query points
	ANNpoint queryPt;

	ANNidxArray	nnIdx;			// 	near neighbour indices
	ANNdistArray 	dists;			//	near neighbour distances
	ANNkd_tree* 	kdTree;			//	search structure
	
	// Sub-volume of domain
	scalar domainVol(0);
	// Total volume of domain
	scalar globalDomainVol(0);	
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << " " << endl;
	Pout << "Domain volume[m^3] = " << domainVol << endl;
	globalDomainVol = domainVol;
	if(domainAve&&Pstream::parRun()) reduce(globalDomainVol, sumOp<scalar>());	
	
	

	// Domain Min/Max
	const pointField& pp = mesh.points();	
	
	// Min, max x-coordinates	
	scalar minX = Foam::min(pp & vector(1,0,0));
	scalar maxX = Foam::max(pp & vector(1,0,0));

	// Min, max y-coordinates		
	scalar minY = Foam::min(pp & vector(0,1,0));
	scalar maxY = Foam::max(pp & vector(0,1,0));

	// Min, max z-coordinates		
	scalar minZ = Foam::min(pp & vector(0,0,1));
	scalar maxZ = Foam::max(pp & vector(0,0,1));
	
	Pout << "x["<<minX<<":"<<maxX<<"]"<<endl;
	Pout << "y["<<minY<<":"<<maxY<<"]"<<endl;
	Pout << "z["<<minZ<<":"<<maxZ<<"]"<<endl;
	
	// Results variables
	// Collisional stress tensor
	symmTensor sigmaCollJI(0,0,0,0,0,0);
	symmTensor sigmaColl[nParticleClass+1];
	symmTensor globalSigmaColl[nParticleClass+1];
	
	// Kinetic stress tensor
    	symmTensor sigmaKinI(0,0,0,0,0,0);
	symmTensor sigmaKin[nParticleClass+1]; // nParticleClass+1 --> for mixture velocity
	symmTensor globalSigmaKin[nParticleClass+1]; 
	
	// Mean velocities in sub-volume
	vector meanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	// 
	vector globalMeanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	int npPartClass[nParticleClass+1];
	
	// Mean velocities in each cell
	volVectorField meanVelCell
	(
             IOobject
             (
        	"meanUp",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), vector(0,0,0) )
	);

	// Number of particle in each cell
	scalarField nPCell(mesh.nCells(),scalar(0));

	// Kinetic tensor in each cell
	volSymmTensorField sigmaKinCell
	(
            IOobject
            (
        	"sigmaKin",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
             ),
            mesh,
            dimensionedSymmTensor( "zero", dimensionSet(1,-1,-2,0,0), symmTensor(0,0,0,0,0,0) )
	 );

	// Collisional tensor in each cell
	volSymmTensorField sigmaCollCell
	(
            IOobject
            (
        	"sigmaColl",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
             ),
            mesh,
            dimensionedSymmTensor( "zero", dimensionSet(1,-1,-2,0,0), symmTensor(0,0,0,0,0,0) )
	 );

	// Solid volume fraction in sub-volumes
	scalar alppInSubVolume(0);

        // Number of particle pairs
        scalar numberOfPairs(0);
		
	// Time loop								
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();
	
	        // Read particle data
		#include "readParticleData.H"
									
		// Search for particles touching boundaries
		#include "createGhostParts.H"	
		
		// Define & Init collision variables
		#include "initCollisionVariables.H" 
						
		if(calcCollision)
		{
			// Neighboring list parameters
			//k = 20; 
			k = neighboringListParameters;
			//if(particleCloud.numberOfParticles()<k) k = particleCloud.numberOfParticles();
			if(nPartInSubVolume < k) k = nPartInSubVolume; 

			// Dimensions, exact OR approximate  
			dim=3; eps = 0;

			// Number of points
			//nPts =  particleCloud.numberOfParticles();		       		
			nPts = nPartInSubVolume;  

			// Allocate 
			queryPt = annAllocPt(dim);
			dataPts = annAllocPts(nPts, dim);
			nnIdx = new ANNidx[k];
			dists = new ANNdist[k];

			// Particle collisional stresses
			sigmaCollJI = symmTensor(0,0,0,0,0,0);				
			for(int iPartClass = 0; iPartClass <= nParticleClass; iPartClass++)
	        	{
				sigmaColl[iPartClass] = symmTensor(0,0,0,0,0,0);
				if(domainAve) globalSigmaColl[iPartClass] = symmTensor(0,0,0,0,0,0);
			}	
		}			

		// Particle mass + total mass
		scalar mass_p(0);
		scalar totalmass_p(0);			

		// Initiate
		int iPartClass = 0;
		//npPartClass[iPartClass] = particleCloud.numberOfParticles();
		npPartClass[iPartClass] = nPartInSubVolume;  

		sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
		globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0); 

		if( bidisperse ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        	{
				// Initiate particle number of each class
				npPartClass[iPartClass] = 0;

				// Initiate velocities
				for(int idir=0; idir<3; idir++) 
				{
					meanVel[iPartClass][idir] = 0 ;
					if(domainAve) globalMeanVel[iPartClass][idir] = 0 ;
				}

				// Initiate particle kinetic stresses
				sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
				// Initiate global particle kinetic stresses
				globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
			}
		}		

		// Create particle list
		//for(int index = 0; index <  particleCloud.numberOfParticles(); index++)
		for(int ii = 0; ii <  particlesInSubDomain.size(); ii++)
		{							
			   label index = particlesInSubDomain[ii];	
			//  Cell ID
			// label cellI = particleCloud.cellIDs()[index][0];
			// if(cellI > -1)
			// {
				if(calcCollision)	
				{
					//dataPts[index][0] = particleCloud.position(index).x();		
					//dataPts[index][1] = particleCloud.position(index).y();		
					//dataPts[index][2] = particleCloud.position(index).z();

					dataPts[ii][0] = particleCloud.position(index).x();		
					dataPts[ii][1] = particleCloud.position(index).y();		
					dataPts[ii][2] = particleCloud.position(index).z();						

					for (int dir=0;dir<3;dir++)
					{	
						fcoll[index][dir] = 0;
					 	ftan[index][dir] = 0;		
 					 	fcap[index][dir] = 0;
 						fvisc[index][dir] = 0;
					}
				}

				// Total velocity of particles
				int iPartClass = 0;;	
				//for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
				for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += mass_p * particleCloud.velocity(index)[idir];			
				totalmass_p += mass_p;

                    		// Total velocity of particles in each cell
                    		label cellID;
                   		cellID = particleCloud.cellIDs()[index][0];
                    		if(cellID>-1) 
				{
					nPCell[cellID] += 1.0;
                    			for(int idir=0; idir<3; idir++)
                    			{
                        			meanVelCell[cellID][idir] += particleCloud.velocity(index)[idir];
                    			}

					if( bidisperse )
					{
						iPartClass = particleCloud.type(index);

						for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
						npPartClass[iPartClass]++;	
					}	
				}
			//}
		}	

		// Global domain average
		if (domainAve) 
		{
			for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	        	{
				for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] = meanVel[iPartClass][idir];

				// Parallel computation
				reduce(globalMeanVel[iPartClass], sumOp<vector>());
				// Normalize
				//for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] /= particleCloud.numberOfParticles();
                        	for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] /= nPartInSubVolume;
                	}
		}			

		// Normalize sub-volume velocities
		if(npPartClass[iPartClass]!=0)
		{
			 for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
		}else
		{
			for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] = 0;
            	}

		if( bidisperse )
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        	{
				for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
			}
		}

        	// Normalize velocity at each cell
        	forAll(mesh.cells(),cellI)
        	{
                    if(nPCell[cellI]!=0)
                    {
                	for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir]/= nPCell[cellI];
                    }else
                    {
                	for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir] = 0.;
                    }
        	}

		if(verboseParRun)
		{
			iPartClass = 0;									
			Info << " " << endl;
			Pout << " Particle class = " << iPartClass << endl;
			Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
			Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
			Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

			if (domainAve && Pstream::master())
			{				
				Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
				Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
				Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
			}

			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        		{
					Info << " " << endl;
					Pout << " Particle class = " << iPartClass << endl;
					Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
					Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
					Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

					if (domainAve && Pstream::master())
					{				
						Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
						Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
						Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
					}

				}
			}

		}

		if(calcCollision)
		{
			// Create Tree structure
			kdTree = new ANNkd_tree(dataPts, nPts, dim);
		}

		// Center particle variables
		vector velJ(0,0,0);
		vector posJ(0,0,0);
		vector omeJ(0,0,0);
		scalar radJ(0);
		label typeJ(0);

		// Neighbouring particle variables
		vector velI(0,0,0);
		vector posI(0,0,0);
		vector omeI(0,0,0);
		scalar radI(0);
		label typeI(0);		

		// Other variables
		label cellI;
		label nPContacts(0);
	
		for(int ii = 0; ii <  particlesInSubDomain.size(); ii++)
		{
			label index_j = particlesInSubDomain[ii];

			// Center particle variables
			velJ  = particleCloud.velocity(index_j);
			posJ  = particleCloud.position(index_j);
			omeJ  = particleCloud.omega(index_j);
			radJ  = particleCloud.radius(index_j);
			typeJ = particleCloud.type(index_j) - 1; // Just nclass starts from "0"			

			if((ii%10000)==0) Info << tab << "PartID = " << ii << " ..." <<endl; 

			if(calcCollision)
			{
				//-Search algorithm
				scalar sqRad = 2.*collisionDp*radJ;

				queryPt[0] = posJ[0];
				queryPt[1] = posJ[1];
				queryPt[2] = posJ[2];

				kdTree->annkFRSearch(
			                		queryPt,			// query point					
							sqRad,				// squared radius
							k,				// number of the near neighbours to return
							nnIdx,				// nearest neighbor array
							dists,				// dist to near neighbours
							eps			);					


				//-Real particles loop
				#include "realPartsLoop.H"

				//-Ghost particles loop
				//#include "ghostPartsLoop.H"

			}	

			//-Calculate stresses
			#include "stresses.H"

		}

		//-Number of contacts
		if(calcCollision) Pout << tab << "Number of contacts = " << numberOfPairs << endl;

		//-Output
                #include "output.H"  

	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);
	particleCloud.dataExchangeM().destroy(types,1);
	particleCloud.dataExchangeM().destroy(omegas,3);
	particleCloud.dataExchangeM().destroy(radii,1);
	particleCloud.dataExchangeM().destroy(cellID,1);
	particleCloud.dataExchangeM().destroy(forces,3);
	
	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
