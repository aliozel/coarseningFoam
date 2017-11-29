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

#include "immiscibleIncompressibleTwoPhaseMixture.H"

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
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

    word patchName = args[1];

    const label patchI = mesh.boundaryMesh().findPatchID(patchName);	

    if (patchI < 0)
    {
            FatalError
            << "Unable to find patch " << patchName << nl
            << exit(FatalError);
    }else
    {
    	Info << "Patch name: " << patchName << endl;
    	Info << " " << endl;
    }
    
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
	IOobject
	(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
	),
	mesh
    );

    #include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    #include "readGravitationalAcceleration.H"

    // Create output folder	
    OFstream* outputFile;
    outputFile =  new OFstream(mesh.time().path()/"massFlowRateAtPatch"+patchName+".dat");
    Info<< "\nWriting mass flow rates into the file " << mesh.time().path()/"massFlowRateAtPatch"+patchName+".dat" << endl;
    *outputFile  << "#Time \t liquidFlowRate \t gasFlowRate \t totalPatchArea" << endl;	

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();
	
	// Read color function
	IOobject alphaheader
	(
		"alpha.water",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	Info<< " Reading alpha" << endl;
	volScalarField alpha(alphaheader,mesh);	

	// Read velocity
	IOobject Umheader
	(
		"U",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);

	Info<< " Reading U" << endl;
	volVectorField Um(Umheader,mesh);
	U = Um;	
	
	// Liquid mass flux
	vector totalLiqMassFlux(0,0,0);

	// Gas mass flux
	vector totalGasMassFlux(0,0,0);
	
	// Total patch area
	scalar totalPatchArea(0);
	
	// Velocity aligned with the gravity direction
	volScalarField UHat((U & g)/mag(g));
	
	totalLiqMassFlux = gSum(     alpha.boundaryField()[patchI]*rho1.value()*UHat.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);
	totalGasMassFlux = gSum((1.-alpha.boundaryField()[patchI])*rho2.value()*UHat.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);
	
	totalPatchArea = gSum(mesh.magSf().boundaryField()[patchI]);		

	*outputFile	<< mesh.time().value() 
			<< tab << " " 
                        << totalLiqMassFlux[0] << " " << totalLiqMassFlux[1] << " " << totalLiqMassFlux[2] 
			<< tab << " " 
			<< totalGasMassFlux[0] << " " << totalGasMassFlux[1] << " " << totalGasMassFlux[2] 
			<< tab << "" 
			<< totalPatchArea
                	<< endl;
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
