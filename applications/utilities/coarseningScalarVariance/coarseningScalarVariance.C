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

#include "boxFilter.H"
#include "convolution.H"
#include "correlationCoeff.H"
#include "conditionalAve.H"
#include "multipleVarsConditionalAve.H"
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

	
        // Post-processing dictionary
        #include "postProcessingDict.H"
		
	// Create conditional averaging class
	conditionalAve condAve(postProcessingDict,conditionalAveragingDict,
						mesh,nVariable,nAveragingVariable,nTotalCase,conditionalAveraging);

	// Create multiple variable conditional averaging class
	multipleVarsConditionalAve multipleVarsCondAve(postProcessingDict,multConditionalAveragingDict,
						mesh,multNVariable,multNAveragingVariable,multNTotalCase,multConditionalAveraging);
	
		
	//-Time loop
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;
		
		if(fluidCoarsening)
		{
		  forAll(filterWidth,fWidth)
		  {
		    //- Create name for filter width			  
		    char charName[100];
		    sprintf(charName, "%dX%dX%d",filterWidth[fWidth],filterWidth[fWidth],filterWidth[fWidth]);
		    word filterWidthName(charName);
			  
		    // Read Eulerian variables
		    #include "readEulerianVar.H"
				    	
		    //- Conditionally averaging		    
		    if(conditionalAveraging)
		    {
			//- Loop over variables
			forAll(variableList,jj)
			{
			   word variableName(variableList[jj]);
			   word varName(filterWidthName);
			   condAve.calc(variableName,varName,averagingVariableList,jj+variableList.size()*fWidth);
			}

	  	    }	
	  	    //- Multi-variable averaging
		    if(multConditionalAveraging) 
		    {
			forAll(multVariableList,jj)
			{
			   word variableName(multVariableList[jj]);
			   word varName(filterWidthName);
			   multipleVarsCondAve.calc(variableName,varName,multAveragingVariableList,jj+multVariableList.size()*fWidth);
			}
		    }
		    		    

		  } 
		}
		
	
	}

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
