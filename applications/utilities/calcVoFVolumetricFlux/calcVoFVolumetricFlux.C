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
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"


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

    //#include "createPhi.H"

    //Info<< "Reading transportProperties\n" << endl;
    //immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    //const dimensionedScalar& rho1 = mixture.rho1();
    //const dimensionedScalar& rho2 = mixture.rho2();

   //#include "readGravitationalAcceleration.H"

    // Create output folder	
    OFstream* outputFile;
    outputFile =  new OFstream(mesh.time().path()/"volumetricFlux.dat");
    Info<< "\nWriting volumetric flux into the file " << mesh.time().path()/"volumetricFlux.dat" << endl;
    *outputFile  << "#Time \t liquidVol \t liquidVolFluxX \t liquidVolFluxY \t liquidVolFluxZ \t gasVolFluxX \t gasVolFluxY \t gasVolFluxZ" << endl;	

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

	scalar liquidVol =  fvc::domainIntegrate(alpha).value();
	vector liquidVolFlux = fvc::domainIntegrate(alpha*U).value();
	vector gasVolFlux = fvc::domainIntegrate((1-alpha)*U).value();
	scalar gasVol =  fvc::domainIntegrate(1.-alpha).value();
	scalar totalVol = gSum(mesh.V());
		
	*outputFile	<< mesh.time().value() 
			<< tab << " " 
			<< liquidVol/totalVol
			<< tab << " " 
                        << liquidVolFlux[0]/liquidVol << " " << liquidVolFlux[1]/liquidVol << " " << liquidVolFlux[2]/liquidVol
			<< tab << " " 
			<<    gasVolFlux[0]/gasVol    << " " <<    gasVolFlux[1]/gasVol    << " " <<    gasVolFlux[2]/gasVol
                	<< endl;
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
