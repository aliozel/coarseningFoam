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
    argList::validArgs.append("Threshold");
    
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    const word patchName = args[1];

    const word charThreshold = args[2];
    
    char* pEnd;
    scalar threshold = strtod(charThreshold.c_str(), &pEnd);

    Info << "Patch name: " << patchName << endl;
    Info << " " << endl;	

/*     // Create output folder	
    fileName outputRelativePath("FreeSurfaceArea");
    if( !isDir(mesh.time().path()/outputRelativePath) )
    {
	mkDir(mesh.time().path()/outputRelativePath );
    } */

    OFstream* outputFile;
    //outputFile =  new OFstream(mesh.time().path()/outputRelativePath/"freeSurfaceArea.dat");
    fileName outpuFilename(mesh.time().path()/"freeWettedSurfaceAreaPatch"+patchName);
    outputFile =  new OFstream(outpuFilename);
    // *outputFile  << "#Time \t FreeSurfaceArea \t TotalWettedArea" << endl;
    
    // Total volume of domain
    scalar totalVol = gSum(mesh.V());

    *outputFile  << "#Total volume of domain = " << totalVol << endl;   
    *outputFile  << "#Time \t totalWettedArea \t freeSurfaceArea \t xmin \t xmax \t ymin \t ymax \t zmin \t zmax \t comx \t comy \t comz" << endl;

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
	
	// Calculate free surface area @ interface
	volScalarField mag_gradAlpha = mag(fvc::grad(alpha));
	scalar totalSurfaceArea = gSum(mag_gradAlpha.internalField()*mesh.V());
    
	// Calculate wetted area on the plane
	//const word patchName = "Bottom";
	const label patchI = mesh.boundaryMesh().findPatchID(patchName);
	
	//min(mesh.C()[0]
	
	scalar min_x = 100.0;
	scalar max_x = -100.0;
	scalar min_y = 100.0;
	scalar max_y = -100.0;
	scalar min_z = 100.0;
	scalar max_z = -100.0;
	vector center_of_mass(0, 0, 0);
	scalar counter = 0;
	
	//loop over all cells to find minimum and maximum cells with liquid
	forAll(mesh.cells(),index)
	{
		if(alpha[index]>threshold)
		{
			if(mesh.C()[index][0]<min_x)
			{
				min_x = mesh.C()[index][0];
			}
			else if(mesh.C()[index][0]>max_x)
			{
				max_x = mesh.C()[index][0];
			}
			
			if(mesh.C()[index][1]<min_y)
			{
				min_y = mesh.C()[index][1];
			}
			else if(mesh.C()[index][1]>max_y)
			{
				max_y = mesh.C()[index][1];
			}
			
			if(mesh.C()[index][2]<min_z)
			{
				min_z = mesh.C()[index][2];
			}
			else if(mesh.C()[index][2]>max_z)
			{
				max_z = mesh.C()[index][2];
			}
		
			center_of_mass = center_of_mass + mesh.C()[index];
			counter = counter+1;
		}
	} 	
	center_of_mass = center_of_mass/counter;
	
	
        if (patchI < 0)
        {
        	FatalError
                << "Unable to find patch " << patchName << nl
                << exit(FatalError);
        }
	
	scalar totalWettedArea = gSum(alpha.boundaryField()[patchI]*mesh.magSf().boundaryField()[patchI]);
	

    	if(Pstream::master()) //Write only if master
    	{
    		
		Info<< " Writing free surface area into the file " << outpuFilename << endl;
		*outputFile	<< alpha.mesh().time().value()  << tab << " " 
                                << totalWettedArea 		<< tab << " " 
				<< totalSurfaceArea 		<< tab << " "
				<< min_x			<< tab << " "
				<< max_x			<< tab << " "
				<< min_y			<< tab << " "
				<< max_y			<< tab << " "
				<< min_z			<< tab << " "
				<< max_z			<< tab << " "
				<< center_of_mass[0]		<< tab << " "
				<< center_of_mass[1]		<< tab << " "
				<< center_of_mass[2]		<< tab << " "
                		<< endl;

    	}
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
