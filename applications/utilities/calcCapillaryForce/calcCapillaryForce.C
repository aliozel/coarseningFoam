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

#   include "readGravitationalAcceleration.H"

    // Transport properties dictionary	
    IOdictionary transportPropertiesDict
    (
	IOobject
	(
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
	)
    );
    dimensionedScalar sigma
    (
    	"sigma",
	dimensionSet(1,0,-2,0,0,0,0),	
	transportPropertiesDict.lookup("sigma")
    ); 	

    const dictionary liqDict(transportPropertiesDict.subDict("water"));\
    dimensionedScalar rho
    (
    	"rho",
	dimensionSet(1,-3,0,0,0,0,0),	
	liqDict.lookup("rho")
    );    

    Info << "Surface tension = " << sigma << endl;
    Info << "Liquid density = " << rho << endl;
    Info << "Gravity = " << g << endl;
    Info << " " << endl;	

    OFstream* outputFile;

    fileName outputFilename(mesh.time().path()/"integratedCappilaryForce.dat");
    outputFile =  new OFstream(outputFilename);
 
    *outputFile  << "#Time \t capForceX \t capForceY  \t capForceZ \t liqWeight" << endl;
    Info<< " Writing integrated capillary force into the file " << outputFilename << endl;

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
	
	// Liquid volume
	scalar volLiq = fvc::domainIntegrate(alpha).value();	
	
    	// grad(alpha)
	volVectorField gradAlpha = fvc::grad(alpha);
	
    	// Unit vector at each cell center
	dimensionedScalar nsDelta
	(
    	    "nsDelta",
	    dimensionSet(0,-1,0,0,0,0,0),	
	    scalar(1.e-64)
	);	
	
	
	
	volVectorField ns(gradAlpha/(Foam::mag(gradAlpha) + nsDelta));

    	// Surface curvature	
    	volScalarField kappa(-fvc::div(ns));
	
	// Capillary force
	volVectorField capForce = sigma.value()*kappa*gradAlpha;
	
	// Integrated over domain
	vector intCapForce = fvc::domainIntegrate(capForce).value();
	


    	if(Pstream::master()) //Write only if master
    	{
    		
		*outputFile	<< alpha.mesh().time().value()  << tab << " " 
                                << intCapForce[0] 		<< tab << " " 
				<< intCapForce[1] 		<< tab << " "
				<< intCapForce[2]		<< tab << " "
				<< rho.value()*volLiq*g.value()[2] << tab << " "
                		<< endl;

    	}
        
	Info << " " << endl;
    }

    return 0;
}


// ************************************************************************* //
