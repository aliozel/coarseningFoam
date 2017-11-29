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
#include "turbulentTransportModel.H"

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::validArgs.append("patchName");
    argList::validArgs.append("firstPatchNumber");
    argList::validArgs.append("lastPatchNumber");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"

    word patchName = args[1];
    const word charfirstPatchNumber = args[2];
    const word charlastPatchNumber = args[3];
    
    char* pEnd;
    int firstPatchNumber = strtod(charfirstPatchNumber.c_str(), &pEnd);  
    int lastPatchNumber  = strtod(charlastPatchNumber.c_str(), &pEnd);  

    Info << "Patch name: " << patchName << " id = ["<< charfirstPatchNumber <<":"<< charlastPatchNumber <<"]" << endl;
    Info << " " << endl;	

    OFstream* outputFile;
    fileName outpuFilename(mesh.time().path()/"domainAveForces"+patchName+".dat");
    outputFile =  new OFstream(outpuFilename);

    Info<< "Forces will be written into the file " << outpuFilename << endl;
 
    *outputFile << "#Time "       << tab 
		<< "weightLiquid" << tab
		<< "weightGas"    << tab	
                << "fLiquidSolid" << tab	
                << "fGasSolid"    << tab	
                << "pLiquidSolid" << tab
                << "pGasSolid"    << tab
                << "fLiquidGas"   << tab
                << "gradPLiq"     << tab
                << "gradPGas"     << tab
                << "fLiquidGas"   << tab
                << "fShear"       << endl;

    OFstream* outputFileCap;
    fileName outpuFilenameCap(mesh.time().path()/"capillaryPressures.dat");
    outputFileCap =  new OFstream(outpuFilenameCap);

    Info<< "Capillary pressures will be written into the file" << outpuFilenameCap << endl;

    *outputFileCap << "#Time " << tab 
		   << "pLiq"   << tab
		   << "pGas"   << endl;

    OFstream* outputFileMom;
    fileName outpuFilenameMom(mesh.time().path()/"mixtureMom.dat");
    outputFileMom =  new OFstream(outpuFilenameMom);

    Info<< "Mixture momentum will be written into the file" << outpuFilenameMom << endl;

    *outputFileMom << "#Time " << tab 
		   << "total momentum"   << endl;
		
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;
	
	// Read phase indicator function
	IOobject alpha1header
	(
		"alpha.water",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);
	Info<< " Reading alpha1" << endl;
	volScalarField alpha1m(alpha1header,mesh);
	alpha1 = alpha1m;	
	
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

        // Read pressure
	IOobject pheader
	(
		"p",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ	
	);
	Info<< " Reading p" << endl;
	volScalarField pm(pheader,mesh);
	p = pm;	

	//- Gas-solid interaction force	
	vector fGasSolid(0,0,0);
        //- Liquid-solid interaction force	
	vector fLiqSolid(0,0,0);
        //- Liquid-gas interaction force	
	vector fLiqGas(0,0,0);
	//- Total weight of liquid
	vector weightLiquid(0,0,0);
	//- Total weight of gas
	vector weightGas(0,0,0);
	//- Shear stress due to pressure
	vector pStressLiq(0,0,0);
	vector pStressGas(0,0,0);
	//- Buoyancy force (gradP)
	vector gradPLiq(0,0,0);
	vector gradPGas(0,0,0);
	//- Liquid-gas interaction force (volume approach)
	vector fLiqGasVol(0,0,0); 
	//- Shear stress contribution of force (volume approach)
	vector fShear(0,0,0); 
	
	//- Weight
	weightLiquid = gSum(     alpha1*rho1.value()*mesh.V()*g.value());
	weightGas    = gSum((1.-alpha1)*rho2.value()*mesh.V()*g.value());

        //- Total  momentum
	volScalarField rho = alpha1*rho1.value() + (1.-alpha1)*rho2.value();
	vector totalMom(0,0,0);
	totalMom = gSum(rho.internalField()*U.internalField()*mesh.V()); 
	
	//- Dynamic viscosity
        mixture.correct();		
	volScalarField mu = mixture.mu();
		
	if(firstPatchNumber>0 && lastPatchNumber>0)
	{
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

			fLiqSolid += gSum(    alpha1.boundaryField()[patchI] *mu.boundaryField()[patchI]*
			                         U.boundaryField()[patchI].snGrad()*mesh.magSf().boundaryField()[patchI]);
						 
			fGasSolid += gSum((1.-alpha1.boundaryField()[patchI])*mu.boundaryField()[patchI]*
			                         U.boundaryField()[patchI].snGrad()*mesh.magSf().boundaryField()[patchI]);
						 
			pStressLiq += gSum(    -alpha1.boundaryField()[patchI] *
			                         p.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);
						 			 
			pStressGas += gSum(-(1.-alpha1.boundaryField()[patchI])*
			                         p.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);	

		}
	}else
	{
		const label patchI = mesh.boundaryMesh().findPatchID(patchName);	

        	if (patchI < 0)
        	{
        		//FatalError
                	//<< "Unable to find patch " << patchName << nl
                	//<< exit(FatalError);
        		Info << "Unable to find patch " << patchName << endl;
		}else{

			fLiqSolid += gSum(    alpha1.boundaryField()[patchI] *mu.boundaryField()[patchI]*
			                	 U.boundaryField()[patchI].snGrad()*mesh.magSf().boundaryField()[patchI]);

			fGasSolid += gSum((1.-alpha1.boundaryField()[patchI])*mu.boundaryField()[patchI]*
			                	 U.boundaryField()[patchI].snGrad()*mesh.magSf().boundaryField()[patchI]);

			pStressLiq += gSum(    -alpha1.boundaryField()[patchI] *
			                	 p.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);	

			pStressGas += gSum(-(1.-alpha1.boundaryField()[patchI])*
			                	 p.boundaryField()[patchI]*mesh.Sf().boundaryField()[patchI]);	
		}
	}	
	
	//- Capillary force 
        //- Update curvature
        mixture.correct();		
	//volVectorField capForce( fvc::reconstruct( mixture.surfaceTensionForce() * mesh.magSf() ) );
	//Info << " mixture.sigma() = " << mixture.sigma() << endl;
	surfaceScalarField surfaceTensionF = mixture.sigma() * fvc::interpolate(mixture.K()) * fvc::snGrad(alpha1);
	volVectorField capForce = fvc::reconstruct( surfaceTensionF * mesh.magSf() );
	fLiqGas = gSum(capForce.internalField()*mesh.V());
		
	//- Capillary force with volume approach	
	//- Cell gradient of alpha1
	volVectorField gradalpha1(fvc::grad(alpha1));
	//- Calculate curvature
	dimensionedScalar deltaN("deltaN",1e-8/pow(average(mesh.V()), 1.0/3.0));
	volVectorField nHat(gradalpha1/(mag(gradalpha1) + deltaN));
	volScalarField K(-fvc::div(nHat,"div(phirb,alpha)"));
	volVectorField capForceVol = mixture.sigma()*K*gradalpha1;
	fLiqGasVol = gSum(capForceVol.internalField()*mesh.V());

	//- Bouyancy force	
	volVectorField buoForceLiq = fvc::reconstruct( - fvc::interpolate(alpha1) * fvc::snGrad(p) * mesh.magSf() );
 	volVectorField buoForceGas = fvc::reconstruct( - fvc::interpolate(1.-alpha1) * fvc::snGrad(p) * mesh.magSf() );

        gradPLiq = gSum(buoForceLiq.internalField()*mesh.V()); 
        gradPGas = gSum(buoForceGas.internalField()*mesh.V()); 
	
	//- Shear stresses (volumetric computation)
	volVectorField shearForce = (fvc::div(mu*fvc::grad(U),"div(phirb,alpha)"));
	fShear = gSum(shearForce.internalField()*mesh.V()); 
		
	//- Capillary pressures
	scalar pLiq(0); 
	scalar pGas(0);
		
	pLiq = fvc::domainIntegrate(alpha1*p).value()/fvc::domainIntegrate(alpha1).value();
	pGas = fvc::domainIntegrate((1.-alpha1)*p).value()/fvc::domainIntegrate((1.-alpha1)).value();
	 
	if(Pstream::master()) // Write only if master
    	{    		
		*outputFile	<< alpha1.mesh().time().value()  << tab 
				<< weightLiquid[0] << " " << weightLiquid[1] << " " << weightLiquid[2] << tab		// 4
				<< weightGas[0]    << " " << weightGas[1]    << " " << weightGas[2]    << tab		// 7	
				<< fLiqSolid[0]    << " " << fLiqSolid[1]    << " " << fLiqSolid[2]    << tab		// 10
				<< fGasSolid[0]    << " " << fGasSolid[1]    << " " << fGasSolid[2]    << tab		// 13	
				<< pStressLiq[0]   << " " << pStressLiq[1]   << " " << pStressLiq[2]   << tab		// 16
				<< pStressGas[0]   << " " << pStressGas[1]   << " " << pStressGas[2]   << tab		// 19
				<< fLiqGas[0]      << " " << fLiqGas[1]      << " " << fLiqGas[2]      << tab		// 22
				<< gradPLiq[0]     << " " << gradPLiq[1]     << " " << gradPLiq[2]     << tab		// 25		
				<< gradPGas[0]     << " " << gradPGas[1]     << " " << gradPGas[2]     << tab		// 28		
				<< fLiqGasVol[0]   << " " << fLiqGasVol[1]   << " " << fLiqGasVol[2]   << tab		// 31
				<< fShear[0]       << " " << fShear[1]       << " " << fShear[2]       << endl;		// 34

		*outputFileCap	<< alpha1.mesh().time().value()  << tab 
				<< pLiq << tab 
				<< pGas << endl;
				
		*outputFileMom	<< alpha1.mesh().time().value()  << tab 
				<< totalMom[0] << " " << totalMom[1] << " " << totalMom[2] << endl;

    	}        
    }

    return 0;
}


// ************************************************************************* //
