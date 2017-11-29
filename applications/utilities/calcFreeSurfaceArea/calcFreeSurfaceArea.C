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

#include "IOmanip.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    // Output file
    OFstream* outputFile;
    fileName outpuFilename(mesh.time().path()/"freeSurfaceArea.dat");
    outputFile =  new OFstream(outpuFilename);
    
    // Total volume of domain
    scalar totalVol = gSum(mesh.V());

    *outputFile  << "#Total volume of domain = " << totalVol << endl;   
    *outputFile  << "#Time \t interfacialArea" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
	
	// Read alpha
	IOobject alphaheader
	(
		"alpha.water",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	Info<< " Reading alpha" << endl;
	volScalarField alpha(alphaheader,mesh);	
	
	// Calculate free surface area @ interface
	volScalarField mag_gradAlpha = mag(fvc::grad(alpha));
	scalar totalSurfaceArea = gSum(mag_gradAlpha.internalField()*mesh.V());
		
    	if(Pstream::master()) //Write only if master
    	{
    		
		Info<< " Writing free surface area into the file " << outpuFilename << endl;
		*outputFile	<< alpha.mesh().time().value()  << tab << " " 
                                << totalSurfaceArea 		<< tab << " " 
                		<< endl;
    	}
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
