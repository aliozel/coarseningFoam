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

//#include "boxFilter.H"
#include "constructFilter.H"
#include "convolution.H"
#include "correlationCoeff.H"

#include "conditionalAve.H"
#include "multipleVarsConditionalAve.H"
#include "threeVarsConditionalAve.H"

#include "dragModelPost.H"
#include "locateModelPost.H"
#include "parcelCloud.H"

#include "meshSearch.H"

#include "vectorList.H"

#include "CPCCellToCellStencil.H"


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
    	#include "readGravitationalAcceleration.H"

	//- Parallel computation
	int me;
	me = Pstream::myProcNo();
	
        // Post-processing dictionary
        #include "postProcessingDict.H"
	
	// Create box filters
	//boxFilter lesFilters(mesh,false);
        autoPtr<constructFilter> lesFilters
        (
            constructFilter::New(mesh,postProcessingDict)
        );

	// Define convolution kernel
	convolution convKernel(mesh,lesFilters,readFilteredVars);	
	// Create correlation coefficient class
	correlationCoeff corrCoeff(mesh);

	// Create conditional averaging class
	conditionalAve condAve(postProcessingDict,conditionalAveragingDict,
						mesh,nVariable,nAveragingVariable,nTotalCase,conditionalAveraging);

	// Create multiple variable conditional averaging class
	multipleVarsConditionalAve multipleVarsCondAve(postProcessingDict,multConditionalAveragingDict,
						mesh,multNVariable,multNAveragingVariable,multNTotalCase,multConditionalAveraging);

	// Create three variable conditional averaging class
	threeVarsConditionalAve threeVarsCondAve(postProcessingDict,threeConditionalAveragingDict,
						mesh,threeNVariable,threeNAveragingVariable,threeNTotalCase,threeConditionalAveraging);
						
	// Create particle drag model class
    	autoPtr<dragModelPost> dragF
	(
    		dragModelPost::New(postProcessingDict,mesh)
	);
	Info << "Particle drag model class was created" << endl;
	
	// Create particle resolved drag model class
    	autoPtr<dragModelPost> resolvedDragF
	(
    		dragModelPost::New(postProcessingDict,mesh)
	);	
	Info << "Particle resolved drag model class was created" << endl;

	// Create Eulerian drag model class
    	autoPtr<dragModelPost> dragFEulerian
	(
    		dragModelPost::New(postProcessingDict,mesh)
	);
	Info << "Eulerian drag model class was created" << endl;
	
	// Create Eulerian resolved drag model class
    	autoPtr<dragModelPost> resolvedDragFEulerian
	(
    		dragModelPost::New(postProcessingDict,mesh)
	);	
	Info << "Eulerian resolved drag model class was created" << endl;
	
	autoPtr<dragModelPost> dragFParcel;
	
	//- Open cpus boundary file
	std::ifstream inputPtrCPUS;
	inputPtrCPUS.open("../DEM/in.cpusPost");

	int ncpus;
        if (me == 0) {	
	 inputPtrCPUS >> ncpus;
	 Info << "\nOpening ../DEM/in.cpus, nCpus = " << ncpus << endl;
        }
        Pstream::scatter(ncpus);
	
	scalarList coords(6*ncpus); 

	if (me == 0) {
	  for(int ii = 0; ii < 6*ncpus; ++ii)
	  {
              inputPtrCPUS >> coords[ii];
	  }
        }
			 				
	//- Global particle velocities
	int nPAll;	
	vector partAllVel(0,0,0);
	vector partAllPos(0,0,0);
	scalar partAllRadius(0);
	
	// Forcases without in.CPUS file
	if( ncpus == 0 ) ncpus += 1 ;	

	// Create output folder	
	if( !isDir(mesh.time().path()/outputRelativePath ) )
	{
		mkDir(mesh.time().path()/outputRelativePath );
	}

	/* - AO
	OFstream* outputFile;
        fileName outpuFilename(mesh.time().path()/outputRelativePath/"vdriftBudgetAnalysisZ.dat");
        outputFile =  new OFstream(outpuFilename);
        *outputFile  << "#Time " << "DivAsgs "  << tab << "Bresolved" << tab << "Bsgs"      << tab << "Cresolved" << tab << "Csgs" << tab 
		     		 << "Dresolved" << tab << "Dsgs"      << tab << "Eresolved" << tab << "Esgs"      << tab << "Fresolved" << tab 
		        	 << "Fsgs"      << tab << "Grav" 	  << tab << "Gresolved" << tab << "Gsgs"      << tab << "Hresolved" << tab 
				 << "Hsgs" << tab
				 << endl;   
	AO - */	
	
	//- Dummy var
	List < scalar > listScalar(16);	 
	//- Dummy var
	List < List < scalar > > listScalarScalar(filterWidth.size(),listScalar);
	//- Domain-averaged resolved & sgs variables
	List < List < List < scalar > > > domainAveVars(timeDirs.size(),listScalarScalar);
	
	Info << " timeDirs.size() = " << timeDirs.size() << endl;
	Info << " filterWidth.size() = " << filterWidth.size() << endl;
	
	//Info << " domainAveVars = " << domainAveVars << endl;
				     	
	//-Time loop
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();
	
		//- Local variables
		//int nP;
		labelList nP(ncpus,0);	
		List<DynamicList<vector> > velocities(ncpus);
		List<DynamicList<vector> > positions(ncpus);		
		List<DynamicList<scalar> > radii(ncpus);	

		labelList cellIDs;	

		// Lists for mapping function
		scalarList subCellList(maxCellPerPart,0.);
		labelList subCellIDsList(maxCellPerPart,-1);

		labelList cellsPerParticle;	
		labelListList neighboringCellIDs;

		List<List<scalar> >  particleWeights;
		List<List<scalar> >  particleVolumes;
		
		// Read Eulerian variables
		#include "readEulerianVar.H"
		
		// Read particle data
		#include "readParticleData.H"

		// Map solid volume fraction and calculate weight functions
		#include "mapSolidVolFracWeightFunc.H"
		
		// Map particle velocity
		#include "mapUs.H"
		
		// Calculate particle drag force
		Info << tab << "Calculating particle drag force..." << endl;
		//dragF().setForce(voidfraction,U,nP[me],cellIDs,positions[me],velocities[me],radii[me]);
		dragF().setForce(mappedVoidfraction,U,nP[me],cellIDs,positions[me],velocities[me],radii[me],particleWeights,neighboringCellIDs);
		
		// Calculate particle gas pressure force
		Info << tab << "Calculating particle gas pressure force..." << endl;
		//dragF().setForce(gradPg,nP[me],cellIDs,positions[me],radii[me]);
		dragF().setForce(gradPg,nP[me],cellIDs,positions[me],radii[me],particleWeights,neighboringCellIDs);
				
		// Calculate Eulerian drag force
		Info << tab << "Calculating Eulerian drag force..." << endl;
	        //dragFEulerian().setForceEulerian(voidfraction,U,radii[me],Us);
		dragFEulerian().setForceEulerian(mappedVoidfraction,U,radii[me],mappedUs);
		
		// Calculate Eulerian drag force
		Info << tab << "Calculating Eulerian gas pressure force.." << endl;
	        //dragFEulerian().setForceEulerian(1.-voidfraction,gradPg);
		dragFEulerian().setForceEulerian(1.-mappedVoidfraction,gradPg);

		// Map drag force
		#include "mapDragForce.H"
		 
		// User loop
		label userLoopLabel(0);
		// Multi-var user loop
		label multUserLoopLabel(0);

		if(fluidCoarsening)
		{
		  forAll(filterWidth,fWidth)
		  {
		    //- Create name for filter width			  
		    char charName[100];
		    sprintf(charName, "%dX%dX%d",filterWidth[fWidth],filterWidth[fWidth],filterWidth[fWidth]);
		    word filterWidthName(charName);
	
		    // Fluid coarsening
		    #include "fluidCoarsening.H"
		   
		    //-Eulerian filtering
		    #include "eulerianVars.H"				    

		    // Conditional averaging
		    #include "conditionalAveraging.H"
		    		    		  
		  } 
		}
		
		// Update mappedVoidfraction and weight function
    		forAll(mesh.cells(),index) 
		{
		   mappedVoidfraction[index] = 1.;		        
		   weightField[index] = 0.;	
		}	
                velocities.clear();
		positions.clear();
		radii.clear();	
	
	}

	//-Write outputs 
	#include "output.H"	  
		    

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
