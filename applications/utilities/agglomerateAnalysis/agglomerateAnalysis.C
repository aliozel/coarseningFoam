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

#include "agglomerateCloud.H"

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

	#include "createFields.H" 	 

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **radii;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.get_radii(radii); 
	
	// Post-processing dictionary
	const dictionary dict(particleCloud.couplingProperties());
		
	forAll(timeDirs, timeI)
	{  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;
				
		int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				
	
		if (count > 0)
		{
		   Info<< " " << endl;
		   particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
   		   Info<< tab <<"Reading particle velocities" << endl;
		   Info<< " " << endl;

		   particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
		   Info<< tab <<"Reading particle positions" << endl;
		   Info<< " " << endl;

		   particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
		   Info<< tab <<"Reading particle radius" << endl;		
		   Info<< " " << endl;
		   
		   particleCloud.setPos(positions);
        	   particleCloud.setVel(velocities);
		   
        	   int index = 0;
                   Pout << tab << "index  = " << index << endl;
                   Pout << tab << "rp     = " << particleCloud.radius(index) << endl;
                   Pout << tab << "Vp     = " << particleCloud.velocity(index) << endl;
                   Pout << tab << "Xp     = " << particleCloud.position(index) << endl;
		   Pout << " " << endl;
		   	
		   // Create agglomerate cloud
		   bool locateAgglomerate(false);
		   agglomerateCloud agglomerates(dict,mesh,particleCloud,count,locateAgglomerate);
		   agglomerates.createAgglomerate();
		   Info << "\tNumber of agglomerates = " << agglomerates.numberOfAgglomerates() << endl;	
		}
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
