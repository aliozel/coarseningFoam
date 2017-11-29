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
    argList::validArgs.append("patchName");
    argList::validArgs.append("firstPatchNumber");
    argList::validArgs.append("lastPatchNumber");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    word patchName = args[1];
    const word charfirstPatchNumber = args[2];
    const word charlastPatchNumber = args[3];
    
    char* pEnd;
    int firstPatchNumber = strtod(charfirstPatchNumber.c_str(), &pEnd);  
    int lastPatchNumber  = strtod(charlastPatchNumber.c_str(), &pEnd);  

    Info << "Patch name: " << patchName << " id = ["<< charfirstPatchNumber <<":"<< charlastPatchNumber <<"]" << endl;
    Info << " " << endl;	

    OFstream* outputFile;
    fileName outpuFilename(mesh.time().path()/"wettedAreaPatch"+patchName+charfirstPatchNumber+"-"+charlastPatchNumber);
    outputFile =  new OFstream(outpuFilename);
 
    *outputFile << "#Time " << tab << "totalWettedArea" << tab;
    for(label id=firstPatchNumber; id<=lastPatchNumber; id++)
    {
    	*outputFile  << "wettedArea-" << id << tab;
    }
    *outputFile << endl;	

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
	
	// Read gas density
	IOobject alphaheader
	(
		"alpha.water",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	Info<< " Reading alpha" << endl;
	volScalarField alpha(alphaheader,mesh);	
	
	scalar totalWettedArea(0);
	scalarField wettedArea(lastPatchNumber+1,scalar(0));
	scalar totalSurfaceArea(0);
	
	for(label id=firstPatchNumber; id<=lastPatchNumber; id++)
	{ 
		// Calculate wetted area on the plane
		word patchNumber;
		std::stringstream ss;
		ss << id;
		patchNumber = ss.str();
		
		const label patchI = mesh.boundaryMesh().findPatchID(patchName+patchNumber);	

        	if (patchI < 0)
        	{
        		FatalError
                	<< "Unable to find patch " << patchName << nl
                	<< exit(FatalError);
        	}

		wettedArea[id] = gSum(alpha.boundaryField()[patchI]*mesh.magSf().boundaryField()[patchI]);
		totalWettedArea += wettedArea[id];
		totalSurfaceArea += gSum(mesh.magSf().boundaryField()[patchI]);
	}

    	if(Pstream::master()) //Write only if master
    	{
    		
		Info<< " Writing mass flow rates into the file " << outpuFilename << endl;
		*outputFile	<< alpha.mesh().time().value()  << tab 
				<< totalWettedArea		<< tab;
		for(label id=firstPatchNumber; id<=lastPatchNumber; id++)
		{ 	
                	*outputFile << wettedArea[id]	<< tab << " "; 
                }		
		*outputFile	<< totalSurfaceArea << endl;

    	}
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
