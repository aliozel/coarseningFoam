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

#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "averagingModel.H"
#include "momCoupleModel.H"
#include "forceModel.H"
#include "IOModel.H"
#include "interpolationCellPoint.H"

#include "timeSelector.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "CPCCellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	void EulerianParticleVelocity
	(
		cfdemCloud& sm,			
		const fvMesh& mesh,
		volVectorField&	Up_,
		const label& iType_
	)
	{						
		
		// Local variables
		vector Up(0,0,0);
		vector position(0,0,0); 
		label cellI(-1);
		scalar dist_s(0);

		// Neighbouring cells
		CPCCellToCellStencil neighbourCells(mesh);		
		scalarField               weightScalar(27,scalar(0.0));
		scalarField		  sumWeights(mesh.cells().size(),scalar(0));
		Field <Field <scalar> >   particleWeights(sm.numberOfParticles(),weightScalar);
		
		for(label index = 0; index < sm.numberOfParticles(); index++)
		{			
			if( iType_ == sm.type(index) )
			{			
				cellI = sm.cellIDs()[index][0];
				if(cellI > -1)
				{
					position = sm.position(index);			    
					Up = sm.velocity(index);

					labelList& cellsNeigh = neighbourCells[cellI];
					sumWeights = 0;
					dist_s = 0;

					forAll(cellsNeigh,jj)
					{
					// Find distances between particle and neighbouring cells					
						dist_s = mag(sm.mesh().C()[cellsNeigh[jj]]-position)/pow(sm.mesh().V()[cellsNeigh[jj]],1./3.);

						if(dist_s <= 0.5)
						{		
							particleWeights[index][jj] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
						}
						else if (dist_s > 0.5 && dist_s <= 1.5)
						{		
							particleWeights[index][jj] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
						}	
						else if (dist_s > 1.5 && dist_s <= 2.5)
						{			
							particleWeights[index][jj] =  pow(2.5-dist_s,4)/24.;
						}
						else
						{		
							particleWeights[index][jj] = 0;
						}

						sumWeights[cellI] += particleWeights[index][jj];
					
					}
				
					forAll(cellsNeigh,jj)
					{	
						if ( sumWeights[cellI] != 0 )
						{
							Up_[cellI] 	         +=  Up*particleWeights[index][jj]/sumWeights[cellI];
						}
						else
						{
							Up_[cellI] 		 = vector(0,0,0);
						}
					}
				}
			}	
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
	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **radii;
	double **cellID;
	double **types;
	
	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
        particleCloud.dataExchangeM().allocateArray(types,0.,1);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
        // Read dictionary + dictionaryProps 
        const dictionary dict(particleCloud.couplingProperties());
	
        // Read total number of particles type
        const dictionary particleGranularPropsDict(dict.subDict("particleGranularTemperatureProps"));
	scalar nType(readScalar(particleGranularPropsDict.lookup("nType")));	

        // Debuging
        bool verbose(false);
        // Particle ID for debuging
        labelList exList;
        if(particleGranularPropsDict.found("verbose"))
        {
                verbose = true;
                exList = labelList(particleGranularPropsDict.lookup("exList"));
        }

	for(label iType = 1 ; iType <= nType; iType++)
	{
		// Create Eulerian particle velocity
		volVectorField Up
		(
        	    IOobject
        	    (
        		"Up",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
        	    ),
        	    mesh,
		    dimensionedVector( "zero", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
		);

		// Mean Eulerian particle velocity
		volVectorField meanUp
		(
        	    IOobject
        	    (
        		"meanUp",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
        	    ),
        	    mesh,
		    dimensionedVector( "zero", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
		);

		// Create Eulerian particlre velocity interpolation variable
		interpolationCellPoint<vector> meanUpxpInt(meanUp);

						
		//- Time loop
		label nTime(0);
		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI);
			Foam::Info << " " << endl;
			Foam::Info << "Time = " << runTime.timeName() << Foam::endl;			
			if( runTime.timeName() != 0 )
			{
				
				mesh.readUpdate();

				// Read particle data
				#include "readParticleData.H"

				// Map particle velocity
				EulerianParticleVelocity
				(
					particleCloud,			
					mesh,
					Up,
					iType
				);

				// Calculate mean Eulerian particle velcoity
				forAll(mesh.cells(),ii)
				{
					meanUp[ii] += Up[ii];
				}
				
				nTime++;
			}		

		}
		// Normalize Eulerian particle velocity
		forAll(mesh.cells(),ii)
		{
			meanUp[ii] /= nTime;
		}	
	
		// Write meanUp into file
		meanUp.write();

		// Mean particle velocity variance
		vector meanUpLagXp(0,0,0);
		scalar meanUpUpLag(0);
		
		label cellI(-1);
		vector position(0,0,0);
		
		// Mean particle velocity fluctuation
		scalar velFluc(0);
		
		//- Time loop to calculate granular temperature
		nTime = 0;
		forAll(timeDirs, timeI)
		{
			runTime.setTime(timeDirs[timeI], timeI);						
			if( runTime.timeName() != 0 )
			{			
				#include "readParticleData.H"
				for(int index = 0; index < particleCloud.numberOfParticles(); index++)
				{
					if( iType == particleCloud.type(index) )
					{
						cellI = particleCloud.cellIDs()[index][0];
						if(cellI > -1)
						{
							position = particleCloud.position(index);
							meanUpLagXp = meanUpxpInt.interpolate(position,cellI);
							meanUpUpLag += 1./2. *(   ( particleCloud.velocity(index)[0] - meanUpLagXp[0] ) * ( particleCloud.velocity(index)[0] - meanUpLagXp[0] )
							        		+ ( particleCloud.velocity(index)[1] - meanUpLagXp[1] ) * ( particleCloud.velocity(index)[1] - meanUpLagXp[1] ) 
							        		+ ( particleCloud.velocity(index)[2] - meanUpLagXp[2] ) * ( particleCloud.velocity(index)[2] - meanUpLagXp[2] )	);	
						}
					}
				}
				nTime++;
			}	

		}
		
		// Normalize and calculate fluctuation
		velFluc = Foam::pow(3./2.*meanUpUpLag/(particleCloud.numberOfParticles()*nTime),0.5);
		
		// Write output on screen
		Pout << tab << "Particle type = " << iType << " Velocity fluc. = " << velFluc << endl;
	
	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
