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

#include "particleKineticStressBi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	// Initiate output file
	void initOutput
	(	
		const int&          nParticleClass_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		const bool&             bidisperse_
	)
	{
	
	    	outputKinFile_ =  new OFstream(mesh_.time().path()/outputRelativePath_/"kineticStress");
		*outputKinFile_ << "#Time  " << tab << " "
		                << "Vol    " <<        " "  
	                	<< "alpp   " << tab << " "
				<< "Kin_XX " << tab << " "
				<< "Kin_YY " << tab << " "
				<< "Kin_ZZ " << tab << " " 
				<< "Kin_XY " << tab << " " 
				<< "Kin_XZ " << tab << " " 
				<< "Kin_YZ " << tab << " "; 		

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_  << "Kin_XX["<<iPartClass<<"] " << tab << " "
						<< "Kin_YY["<<iPartClass<<"] " << tab << " "
						<< "Kin_ZZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XY["<<iPartClass<<"] " << tab << " " 
						<< "Kin_XZ["<<iPartClass<<"] " << tab << " " 
						<< "Kin_YZ["<<iPartClass<<"] " << tab << " "; 
			}
		}
		*outputKinFile_ << endl;
	
	}
	
	// Write into output file
	void writeOutput
	(	
		const int&          nParticleClass_,
		const Time& 		   runTime_,
		const fvMesh&	       	      mesh_,
		fileName&       outputRelativePath_,
		OFstream*&           outputKinFile_,
		const bool&             bidisperse_,
		SymmTensor<double>*       sigmaKin_,					
		const scalar&		 domainVol_,
		const scalar&	   alppInSubVolume_	)
	{
	
		// Write output
		Info << " " << endl;
		Pout << tab << "Writing particle kinetic stresses into the file " <<nl
		     << tab << tab << mesh_.time().path()/outputRelativePath_/"kineticStress" << endl;

		int iPartClass = 0;
		*outputKinFile_ << runTime_.value() 		       << tab << " " 
		                << domainVol_ 			       << tab << " "
				<< alppInSubVolume_		       << tab << " "   	
				<< sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
				<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
				<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ;								

		if( bidisperse_ ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass_; iPartClass++)
	        	{
				*outputKinFile_ << sigmaKin_[iPartClass][0]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][1]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][2]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][3]/domainVol_ << tab << " " 
						<< sigmaKin_[iPartClass][4]/domainVol_ << tab << " " 		
						<< sigmaKin_[iPartClass][5]/domainVol_ << tab << " " ; 
			}
		}
		*outputKinFile_ << endl;	
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
	double **omegas;
	double **radii;
	double **cellID;	
	double **types;
					
	//- Read post-processing dicitonary
	#include "postProcessingDict.H"

	// Initiate output file
	OFstream* outputKinFile;

    	initOutput(	    nParticleClass,
	                              mesh,
			outputRelativePath,
			     outputKinFile,
			        bidisperse	);

	// Initiate output file
	OFstream* outputGlobalKinFile;

	if(domainAve)
	{
    	  initOutput(	      nParticleClass,
	                        	mesh,
	            outputGlobalRelativePath,
			 outputGlobalKinFile,
			          bidisperse	);
	}
	
	// Sub-volume of domain
	scalar domainVol(0);
	// Total volume of domain
	scalar globalDomainVol(0);	
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << " " << endl;
	Pout << "Domain volume[m^3] = " << domainVol << endl;
	globalDomainVol = domainVol;
	if(domainAve&&Pstream::parRun()) reduce(globalDomainVol, sumOp<scalar>());	
	
	// Domain Min/Max
	const pointField& pp = mesh.points();	
	
	// Min, max x-coordinates	
	scalar minX = Foam::min(pp & vector(1,0,0));
	scalar maxX = Foam::max(pp & vector(1,0,0));

	// Min, max y-coordinates		
	scalar minY = Foam::min(pp & vector(0,1,0));
	scalar maxY = Foam::max(pp & vector(0,1,0));

	// Min, max z-coordinates		
	scalar minZ = Foam::min(pp & vector(0,0,1));
	scalar maxZ = Foam::max(pp & vector(0,0,1));
	
	Pout << "x["<<minX<<":"<<maxX<<"]"<<endl;
	Pout << "y["<<minY<<":"<<maxY<<"]"<<endl;
	Pout << "z["<<minZ<<":"<<maxZ<<"]"<<endl;
	
	// Results variables	
	// Kinetic stress tensor
    	symmTensor sigmaKinI(0,0,0,0,0,0);
	symmTensor sigmaKin[nParticleClass+1]; // nParticleClass+1 --> for mixture velocity
	symmTensor globalSigmaKin[nParticleClass+1]; 
	
	// Mean velocities in sub-volume
	vector meanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	// 
	vector globalMeanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	int npPartClass[nParticleClass+1];
	
	// Mean velocities in each cell
	volVectorField meanVelCell
	(
             IOobject
             (
        	"meanUp",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), vector(0,0,0) )
	);

	// Number of particle in each cell
	// scalarField nPCell(mesh.nCells(),scalar(0));
	
	// Number of particles in each cell
	volScalarField nPCell
	(
             IOobject
             (
        	"nP",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
	);

	// Kinetic tensor in each cell
	volSymmTensorField sigmaKinCell
	(
            IOobject
            (
        	"sigmaKin",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
             ),
            mesh,
            dimensionedSymmTensor( "zero", dimensionSet(1,-1,-2,0,0), symmTensor(0,0,0,0,0,0) )
	 );

	// Solid volume fraction in sub-volumes
	scalar alppInSubVolume(0);
	// Solid volume fraction in domain
	scalar globalAlppInSubVolume(0);
	
	// Time loop								
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();
	
	        // Read particle data
		#include "readParticleData.H"

		// Init & calculate mean Vels.
		#include "initCalcMeanVels.H"
                
		// Loop over particles  	   
		for(int ii = 0; ii <  particlesInSubDomain.size(); ii++)
		{
		   //-Calculate stresses
		   #include "stresses.H"
		}
		//-Output
                #include "output.H"  

	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);
	particleCloud.dataExchangeM().destroy(types,1);
	particleCloud.dataExchangeM().destroy(omegas,3);
	particleCloud.dataExchangeM().destroy(radii,1);
	particleCloud.dataExchangeM().destroy(cellID,1);
	
	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
