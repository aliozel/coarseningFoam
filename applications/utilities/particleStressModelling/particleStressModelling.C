/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 
Application
    Particle stress modelling for parcel approach

Description
    Conditional averaging of particulate phase properties such as; pressure,
    viscosity...

\*---------------------------------------------------------------------------*/


#include "particleStressModelling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

		
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	

	timeSelector::addOptions();
	#   include "addRegionOption.H"
	argList::addBoolOption
	(
	"noWrite",
	"suppress writing results"
	);
	#include "addDictOption.H"

	#include "setRootCase.H"
	#include "createTime.H"
	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

	#include "createFields.H" 	 

    	// Create particle cloud
    	cfdemCloud particleCloud(mesh);

        // Post-processing dictionary
        #include "postProcessingDict.H"
    
	// Create conditional averaging class
	conditionalAve condAve(postProcessingDict,conditionalAveragingDict,
						mesh,nVariable,nAveragingVariable,nTotalCase,conditionalAveraging);
	
	// Create multiple variable conditional averaging class
	multipleVarsConditionalAve multipleVarsCondAve(postProcessingDict,multConditionalAveragingDict,
						mesh,multNVariable,multNAveragingVariable,multNTotalCase,multConditionalAveraging);
			
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Pout << " " << endl;
		Pout << "\nTime = " << runTime.timeName() << endl;

		mesh.readUpdate();

        	// Read gas volume fraction
        	IOobject voidfractionheader
        	(
        	    "voidfraction",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::MUST_READ
        	);

        	Info<< " 	Reading voidfraction" << endl;
        	volScalarField voidfraction(voidfractionheader,mesh);

		// Read Eulerian particle velocity
		IOobject Usheader
		(
		   "Us",
		   runTime.timeName(),
		   mesh,
		   IOobject::MUST_READ	
		);

		Info<< " 	Reading Us" << endl;
		volVectorField Us(Usheader,mesh);

        	// Read particle kinetic stresses
        	IOobject sigmaKinHeader
        	(
        	    "sigmaKin",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::MUST_READ
        	);

        	Info<< " 	Reading sigmaKin" << endl;
        	volSymmTensorField sigmaKin(sigmaKinHeader,mesh);

        	// Read particle collisional stresses
        	IOobject sigmaCollHeader
        	(
        	    "sigmaColl",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::MUST_READ
        	);

        	Info<< " 	Reading sigmaColl" << endl;
        	volSymmTensorField sigmaColl(sigmaCollHeader,mesh);

        	// Particle pressure
        	volScalarField Pp
        	(
        	    IOobject
        	    (
                	"Pp",
                	runTime.timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::AUTO_WRITE
        	    ),
        	    mesh,
        	    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
        	 );

        	Pp = 1./3. * tr( sigmaKin + sigmaColl ) ;
        	// Write into the results folder
        	Pp.write();
        	// Calculate the particulate pressure gradient
        	volVectorField gradPp(fvc::grad(Pp));

        	// Particle shear stress
        	volTensorField sigmap
        	(
        	    IOobject
        	    (
                	"sigmap",
                	runTime.timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::AUTO_WRITE
        	    ),
        	    mesh,
        	    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
        	);

        	sigmap = ( sigmaKin + sigmaColl ) - tensor(I) * Pp;
        	// Write into the results folder
        	sigmap.write();

        	// Particle viscosity
        	volScalarField mup
        	(
        	    IOobject
        	    (
                	"mup",
                	runTime.timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::AUTO_WRITE
        	    ),
        	    mesh,
        	    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
        	);

        	// Particle shear stresses
        	volTensorField S("S",fvc::grad(Us) + fvc::grad(Us)().T());
        	dimensionedScalar SSsmall("zero", dimensionSet(0,0,-2,0,0,0,0), SMALL);
        	mup = ( sigmap && S ) / ( max ( S && S, SSsmall ) );
        	// Limit by zero
        	mup.max(0.);
        	// Write into the results folder
        	mup.write(); 
		
		// Particle viscosity/sqrt(p)
        	dimensionedScalar PpSmall("zero", dimensionSet(1,-1,-2,0,0,0,0), 1.e-06);
        	volScalarField mupSqrtPp
        	(
        	    IOobject
        	    (
                	"mupSqrtPp",
                	runTime.timeName(),
                	mesh,
                	IOobject::NO_READ,
                	IOobject::AUTO_WRITE
        	    ),
        	    mup/sqrt((mag(Pp)+PpSmall)*rhop)/dp
        	); 
		
		//- Dummy word
		word varName("");

		//- Conditional averaging
	        condAve.calc();
		condAve.write(varName);
		
		//- Multi-variable conditional averaging
	        multipleVarsCondAve.calc();
		multipleVarsCondAve.write(varName);

	}        

	Info	<<  " " << endl;
	Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        	<< nl << endl;

	Info<< "End\n" << endl;
	
	

}


// ************************************************************************* //
