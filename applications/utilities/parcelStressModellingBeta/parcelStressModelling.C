/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 
Application
    Parcel stress modelling 

Description
    Conditional averaging of particulate phase properties such as; pressure,
    viscosity...

\*---------------------------------------------------------------------------*/


#include "parcelStressModelling.H"

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

	//- Parallel computation
	int me;
	me = Pstream::myProcNo();	
   
        // Post-processing dictionary
        #include "postProcessingDict.H"

	// Create box filters
	boxFilter lesFilters(mesh,false);
	// Define convolution kernel
	convolution convKernel(mesh,lesFilters,readFilteredVars);
    
	// Create conditional averaging class
	conditionalAve condAve(postProcessingDict,conditionalAveragingDict,
						mesh,nVariable,nAveragingVariable,nTotalCase,conditionalAveraging);
		
	// Create multiple variable conditional averaging class
	multipleVarsConditionalAve multipleVarsCondAve(postProcessingDict,multConditionalAveragingDict,
						mesh,multNVariable,multNAveragingVariable,multNTotalCase,multConditionalAveraging);

	// Create three variable conditional averaging class
	threeVarsConditionalAve threeVarsCondAve(postProcessingDict,threeConditionalAveragingDict,
						mesh,threeNVariable,threeNAveragingVariable,threeNTotalCase,threeConditionalAveraging);

	// Create conditional averaging class for coarsening
	conditionalAve condAveCoarsening(postProcessingDict,conditionalAveragingDictCoarsening,
						mesh,nVariableCoarsening,nAveragingVariableCoarsening,nTotalCaseCoarsening,conditionalAveragingCoarsening);

	// Create multiple variable conditional averaging class
	multipleVarsConditionalAve multipleVarsCondAveCoarsening(postProcessingDict,multConditionalAveragingDictCoarsening,
						mesh,multNVariableCoarsening,multNAveragingVariableCoarsening,multNTotalCaseCoarsening,multConditionalAveragingCoarsening);

	// Create three variable conditional averaging class
	threeVarsConditionalAve threeVarsCondAveCoarsening(postProcessingDict,threeConditionalAveragingDictCoarsening,
						mesh,threeNVariableCoarsening,threeNAveragingVariableCoarsening,threeNTotalCaseCoarsening,threeConditionalAveragingCoarsening);

	//- Open cpus boundary file
	std::ifstream inputPtrCPUS;
	inputPtrCPUS.open("../DEM/in.cpusPost");

	int ncpus;
        if (me == 0) {	
	 inputPtrCPUS >> ncpus;
	 Info << "\nOpening ../DEM/in.cpusPost, nCpus = " << ncpus << endl;
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
	scalar partAllCoordNumber(0);

	//- Local variables
	//int nP;
	// Forcases without in.CPUS file
	if( ncpus == 0 ) ncpus += 1 ;	
	//labelList nP(ncpus,0);	
	List<DynamicList<vector> > velocities(ncpus);
	List<DynamicList<vector> > positions(ncpus);		
	List<DynamicList<scalar> > radii(ncpus);	
	List<DynamicList<scalar> > coordNumbers(ncpus);	
	
	labelList cellIDs;
	// Lists for mapping function
	scalarList subCellList(maxCellPerPart,0.);
	labelList subCellIDsList(maxCellPerPart,-1);

	labelList cellsPerParticle;	
	labelListList neighboringCellIDs;
	
	List<List<scalar> >  particleWeights;
	List<List<scalar> >  particleVolumes;
						
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Pout << " " << endl;
		Pout << "\nTime = " << runTime.timeName() << endl;

		mesh.readUpdate();

		//- Local variables
		//int nP;
		// Forcases without in.CPUS file
		labelList nP(ncpus,0);	
		List<DynamicList<vector> > velocities(ncpus);
		List<DynamicList<vector> > positions(ncpus);		
		List<DynamicList<scalar> > radii(ncpus);	
		List<DynamicList<scalar> > coordNumbers(ncpus);	

		labelList cellIDs;
		// Lists for mapping function
		scalarList subCellList(maxCellPerPart,0.);
		labelList subCellIDsList(maxCellPerPart,-1);

		labelList cellsPerParticle;	
		labelListList neighboringCellIDs;

		List<List<scalar> >  particleWeights;
		List<List<scalar> >  particleVolumes;
								
		// Read particle data
		#include "readParticleData.H"
		
		// Calculate mean solid velocity for each cell
		#include "calcMeanUpCell.H"		

		// Map solid volume fraction and calculate weight functions
		#include "mapSolidVolFracWeightFunc.H"

		// Map particle velocity
		#include "mapUs.H"

		// Map particle coordination number
		#include "mapCoordNumber.H"
		
		// Read Eulerian variables
		#include "readEulerianVars.H"
		
		// Kinematic stress modeling
		#include "kinematicStressModeling.H"
				
		// Parcel stress modeling
		#include "stressModeling.H" 
		
		// Conditional averaging KT
		#include "conditionalAveraging.H"
        					
		// User loop
		label userLoopLabel(0);
		
		//- Stress coarsening
		if(fluidCoarsening)
		{
		  forAll(filterWidth,fWidth)
		  {
		    //- Create name for filter width			  
		    char charName[100];
		    sprintf(charName, "%dX%dX%d",filterWidth[fWidth],filterWidth[fWidth],filterWidth[fWidth]);
		    word filterWidthName(charName);
			  
		    // Meso-scale coarsening
		    #include "coarsening.H"
		    
		    // Coarse kinetic theory based data
		    #include "coarseningKT.H"
		    
		    // Conditional averaging
		    #include "conditionalAveragingCoarsening.H"
		    	
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
		coordNumbers.clear();
	}        

	Info	<<  " " << endl;
	Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        	<< nl << endl;

	Info<< "End\n" << endl;
	
	

}


// ************************************************************************* //
