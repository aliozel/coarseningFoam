/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Calculate multipleVarsConditionalAve

-----------------------------------------------------------------------------*/

#include "multipleVarsConditionalAve.H"

namespace Foam {
    defineTypeNameAndDebug(multipleVarsConditionalAve,0);
}

Foam::multipleVarsConditionalAve::multipleVarsConditionalAve
(
    const dictionary& dict,
    const dictionary& subDict,
    const objectRegistry& obr,
    const label& nVariable,
    const label& nAveragingVariable,
    const label& nTotalCase,
    const bool& multConditionalAveraging 
)
:
    dict_(dict),
    obr_(obr),
    nVariable_(nVariable),
    nAveragingVariable_(nAveragingVariable),
    nTotalCase_(nTotalCase), 
    //multipleVarsConditionalAveDict_(dict_.subDict("multipleVarsConditionalAve")),
    multipleVarsConditionalAveDict_(subDict),
    //outputPath_(multipleVarsConditionalAveDict_.lookup("outputRelativePath")),	
    outputPath_(),
    variableList(nVariable,""),
    averagingVariableList(nAveragingVariable,""),
    nBin(nAveragingVariable,1),
    minCondVar(nAveragingVariable,scalar(-1.e+12)),
    maxCondVar(nAveragingVariable,scalar( 1.e+12)),
    logBin(nAveragingVariable,false),
    nReal(nTotalCase,scalarField(1,scalar(0))),	
    condAveVar(nTotalCase,scalarField(1,scalar(0))),		
    readFromFile_(false),
    multConditionalAveraging_(multConditionalAveraging)
{
	// Init
	init();	
}

Foam::multipleVarsConditionalAve::~multipleVarsConditionalAve()
{}

void Foam::multipleVarsConditionalAve::init()
{	
	// Create output folder
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	
	if(multConditionalAveraging_)
	{
		//multipleVarsConditionalAveDict_ = dict_.subDict("multiVarsConditionalAve");
		outputPath_ = fileName(multipleVarsConditionalAveDict_.lookup("outputRelativePath"));
		
		if(Pstream::master())
		{
			if(mesh.time().processorCase())
			{
				if(!isDir(mesh.time().path()/".."/outputPath_)) mkDir(mesh.time().path()/".."/outputPath_);	
			}else
			{
				if(!isDir(mesh.time().path()/outputPath_)) mkDir(mesh.time().path()/outputPath_);	
			}
		}

		// Read varibles list and dictionaries
        	variableList = wordList(multipleVarsConditionalAveDict_.lookup("variableList"));
		Info << "Multi-variable conditional averaging" << endl;
        	Info << "The variable(s) \n";
        	forAll(variableList,ii)
        	{
        	    Info << tab << variableList[ii];
        	}
        	Info << "\nwill be conditional averaged by the following variable(s):" << endl;
        	averagingVariableList = wordList(multipleVarsConditionalAveDict_.lookup("averagingVariableList"));
        	forAll(averagingVariableList,ii)
        	{
        	    word& condVariableName = averagingVariableList[ii];
        	    Info << tab << condVariableName;
        	}
        	Info << endl;
		//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );
		label nAnalyses = nTotalCase_; // nAveragingVariable_;
		
		label caseN(0);
        	forAll(averagingVariableList,ii)
        	{
        	    word& currentVar = averagingVariableList[ii];
        	    const dictionary& subDict = multipleVarsConditionalAveDict_.subDict(currentVar);
        	    nBin[ii] = readLabel(subDict.lookup(word("nBin")));

        	    minCondVar[ii] = readScalar(subDict.lookup(word("min")));
        	    maxCondVar[ii] = readScalar(subDict.lookup(word("max")));
		    
		    if( subDict.found(word("logBin")) )
		    {
		    	logBin[ii] = true;
		        Info << "Variable " << currentVar << " will be conditionally averaged in log-scale " << endl; 
		    }
/*         	    for(label nA = 0; nA < nAnalyses; nA++)
		    { 
			    //- Resize variables
                	    condAveVar[caseN].resize(nBin[ii]);
                	    nReal[caseN].resize(nBin[ii]);

			    //- Re-initialise 
			    condAveVar[caseN] = scalarField(nBin[ii],scalar(0));
			    nReal[caseN] = scalarField(nBin[ii],scalar(0));

			    caseN++;
		    }	 */

        	}
				
        	for(label nA = 0; nA < nAnalyses; nA++)
		{ 
			//- Resize variables
                	condAveVar[caseN].resize(nBin[0]*nBin[1]);
                	nReal[caseN].resize(nBin[0]*nBin[1]);

			//- Re-initialise 
			condAveVar[caseN] = scalarField(nBin[0]*nBin[1],scalar(0));
			nReal[caseN] = scalarField(nBin[0]*nBin[1],scalar(0));

			caseN++;
		}	
	
		Info << " " << endl;
	}
	else
	{
		Info << "Multi-variable conditional averaging will NOT be performed" << endl;
	}

}

//- Condition averaging for volScalarField
void Foam::multipleVarsConditionalAve::condAve
(   
	 const fvMesh&         	mesh_,
	 const volScalarField& 	var_,
	 const volScalarField& 	condVarI_,
	 const scalar&     	minCondVarI_,
	 const scalar&     	maxCondVarI_,
	 const label&       	nBinI_,
	 const volScalarField& 	condVarJ_,
	 const scalar&     	minCondVarJ_,
	 const scalar&     	maxCondVarJ_,
	 const label&       	nBinJ_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
)
{ 
	 // Bin interval size
	 scalar binSizeI = ( maxCondVarI_ - minCondVarI_ ) / nBinI_;
	 scalar binSizeJ = ( maxCondVarJ_ - minCondVarJ_ ) / nBinJ_;

	 // Loop over cells
	 forAll(mesh_.cells(),cellI)
	 {
	     label binI = floor(( condVarI_[cellI] - minCondVarI_ ) / binSizeI );
	     label binJ = floor(( condVarJ_[cellI] - minCondVarJ_ ) / binSizeJ );
	     if ( binI > -1 && binI < nBinI_ && binJ > -1 && binJ < nBinJ_ )
	     {
        	 condAveVar_[binI+binJ*nBinI_] += var_[cellI];
        	 nReal_[binI+binJ*nBinI_] += 1.;
	     }
	 }
}

//- Condition averaging for user define scalar field
void Foam::multipleVarsConditionalAve::condAveUserDefine
(   
	 const scalarField& 	var_,
	 const scalarField& 	condVarI_,
	 const scalar&     	minCondVarI_,
	 const scalar&     	maxCondVarI_,
	 const label&       	nBinI_,
	 const scalarField& 	condVarJ_,
	 const scalar&     	minCondVarJ_,
	 const scalar&     	maxCondVarJ_,
	 const label&       	nBinJ_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
)
{
	 // Bin interval size
	 scalar binSizeI = ( maxCondVarI_ - minCondVarI_ ) / nBinI_;
	 scalar binSizeJ = ( maxCondVarJ_ - minCondVarJ_ ) / nBinJ_;

	 // Loop over cells
	 forAll(var_,cellI)
	 {
	     label binI = floor(( condVarI_[cellI] - minCondVarI_ ) / binSizeI );
	     label binJ = floor(( condVarJ_[cellI] - minCondVarJ_ ) / binSizeJ );
	     if ( binI > -1 && binI < nBinI_ && binJ > -1 && binJ < nBinJ_ )
	     {
        	 condAveVar_[binI+binJ*nBinI_] += var_[cellI];
        	 nReal_[binI+binJ*nBinI_] += 1.;
	     }
	 }
}

void Foam::multipleVarsConditionalAve::calc()
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	forAll(variableList,jj)
	{
            word& variableName = variableList[jj];
            volScalarField* ptrVar;
            // Look up for variable
            if (mesh.foundObject<volScalarField>(variableName))
            {
        	ptrVar = new volScalarField(mesh.lookupObject<volScalarField>(variableName));
            }else
            {
        	// Look up for expression
        	const dictionary& subDictVariableName = multipleVarsConditionalAveDict_.subDict(variableName);
        	// Parse expression and write into result folder
        	expressionField var(variableName,mesh,subDictVariableName,false);
        	// Read volume field
        	IOobject varHeader
        	(
                    variableName,
                    //runTime.timeName(),
		    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ
        	);
        	ptrVar = new volScalarField(varHeader,mesh);
            }

            //forAll(averagingVariableList,ii)
            //{
		int ii = 0;
		word& condVariableNameI = averagingVariableList[ii];
        	const dictionary& subDictCondVariableNameI = multipleVarsConditionalAveDict_.subDict(condVariableNameI);
//        	Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameI 
//		     << " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] and "  
//		     << " [" << minCondVar[ii+1] << "," << maxCondVar[ii+1] << "] ... " << endl;  

        	
		// Parse expression and write into result folder
        	expressionField condVarI(condVariableNameI,mesh,subDictCondVariableNameI,false);

        	// Read averaging volume field
        	IOobject condVarHeaderI
        	(
                    condVariableNameI,
                    //runTime.timeName(),
		    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ
        	 );
        	volScalarField condVarFieldI(condVarHeaderI,mesh);

		ii++;
		word& condVariableNameJ = averagingVariableList[ii];
        	const dictionary& subDictCondVariableNameJ = multipleVarsConditionalAveDict_.subDict(condVariableNameJ);
        	//Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameJ 
		//     << " [" << minCondVar[ii+jj] << "," << maxCondVar[ii+jj] << "] ..."  << endl;

        	Info <<  "\nVariable " << variableName << " conditionally averaging by "  
		     << condVariableNameI << " [" << minCondVar[0] << "," << maxCondVar[0] << "] and "  
		     << condVariableNameJ << " [" << minCondVar[1] << "," << maxCondVar[1] << "] ... "   << endl;  

        	// Parse expression and write into result folder
        	expressionField condVarJ(condVariableNameJ,mesh,subDictCondVariableNameJ,false);

        	// Read averaging volume field
        	IOobject condVarHeaderJ
        	(
                    condVariableNameJ,
                    //runTime.timeName(),
		    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ
        	 );
        	volScalarField condVarFieldJ(condVarHeaderJ,mesh);

        	// Conditional averaging subroutine
        	volScalarField varField = *ptrVar;
		
        	//condAve( mesh, varField, condVarField, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[ii], nReal[ii] );
		condAve( mesh, varField, condVarFieldI, minCondVar[0], maxCondVar[0], nBin[0], 
					 condVarFieldJ, minCondVar[1], maxCondVar[1], nBin[1], 
					 condAveVar[jj], nReal[jj] );
            //}
	}
}

void Foam::multipleVarsConditionalAve::calc
(
	word& variableName,
	word& userDefineName,
	wordList& averagingVariableListCoarsening,
	int jj
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* ptrVar;
        // Look up for variable
        if (mesh.foundObject<volScalarField>(variableName))
        {
            ptrVar = new volScalarField(mesh.lookupObject<volScalarField>(variableName));
        }else
        {
            // Look up for expression
            const dictionary& subDictVariableName = multipleVarsConditionalAveDict_.subDict(variableName);
            // Parse expression and write into result folder
            expressionField var(variableName,mesh,subDictVariableName,false);
            // Read volume field
            IOobject varHeader
            (
                variableName,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
            );
            ptrVar = new volScalarField(varHeader,mesh);
        }

        //forAll(averagingVariableList,ii)
        //{
	    int ii = 0;
	    word& condVariableNameI = averagingVariableListCoarsening[ii];
            const dictionary& subDictCondVariableNameI = multipleVarsConditionalAveDict_.subDict(condVariableNameI);
//        	Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameI 
//		     << " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] and "  
//		     << " [" << minCondVar[ii+1] << "," << maxCondVar[ii+1] << "] ... " << endl;  


	    // Parse expression and write into result folder
            expressionField condVarI(condVariableNameI,mesh,subDictCondVariableNameI,false);

            // Read averaging volume field
            IOobject condVarHeaderI
            (
                condVariableNameI,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
            volScalarField condVarFieldI(condVarHeaderI,mesh);

	    ii++;
	    word& condVariableNameJ = averagingVariableListCoarsening[ii];
            const dictionary& subDictCondVariableNameJ = multipleVarsConditionalAveDict_.subDict(condVariableNameJ);
            //Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameJ 
	    //     << " [" << minCondVar[ii+jj] << "," << maxCondVar[ii+jj] << "] ..."  << endl;

            Info <<  "\nVariable " << variableName << " conditionally averaging by "  
		 << condVariableNameI << " [" << minCondVar[0] << "," << maxCondVar[0] << "] and "  
		 << condVariableNameJ << " [" << minCondVar[1] << "," << maxCondVar[1] << "] ... "   << endl;  

            // Parse expression and write into result folder
            expressionField condVarJ(condVariableNameJ,mesh,subDictCondVariableNameJ,false);

            // Read averaging volume field
            IOobject condVarHeaderJ
            (
                condVariableNameJ,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
            volScalarField condVarFieldJ(condVarHeaderJ,mesh);

            // Conditional averaging subroutine
            volScalarField varField = *ptrVar;

            //condAve( mesh, varField, condVarField, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[ii], nReal[ii] );
	    condAve( mesh, varField, condVarFieldI, minCondVar[0], maxCondVar[0], nBin[0], 
				     condVarFieldJ, minCondVar[1], maxCondVar[1], nBin[1], 
				     condAveVar[jj], nReal[jj] );

	    //- Write into file
	    writeOutput( mesh, outputPath_, variableName, 
	    		 averagingVariableListCoarsening[0], minCondVar[0], maxCondVar[0], nBin[0], 
	    		 averagingVariableListCoarsening[1], minCondVar[1], maxCondVar[1], nBin[1],
	    		 condAveVar[jj], nReal[jj], userDefineName );   
        //}

}

void Foam::multipleVarsConditionalAve::calcUserDefine
(
	scalarField& varField,
	autoPtr<dragModelPost>& dragF,
	word& variableName,
	label& userLoopVar
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );			    
	label nAnalyses = nTotalCase_ / nAveragingVariable_;	

	//forAll(averagingVariableList,ii)
        //{
            int ii = 0; 
	    word& condVariableNameI = averagingVariableList[ii];
            const dictionary& subDictCondVariableNameI = multipleVarsConditionalAveDict_.subDict(condVariableNameI);
            Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameI << 
		              " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;

            // Parse expression and write into result folder
            expressionField condVarI(condVariableNameI,mesh,subDictCondVariableNameI,false);
		
            // Read averaging volume field
            IOobject condVarHeaderI
            (
                condVariableNameI,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
	     	     
            volScalarField condVarFieldI(condVarHeaderI,mesh);
	    //- Create interpolation class of conditional variable 
	    scalarField condVarFieldXpI(dragF().nP(),0);
	    interpolationCellPoint<scalar> condVarFieldInterpolatorI(condVarFieldI);
	    label nPartInI(0);
	    for(int index=0; index < dragF().nP(); index++)
	    {
		    label cellI = dragF().cellID(index);
		    vector position = dragF().xp(index); 
		    //- Needed for parallel computation
		    if( cellI > -1 ) 
		    {
		    	//condVarFieldXp[nPartIn] = condVarFieldInterpolator.interpolate(position,cellI);
		    	condVarFieldXpI[nPartInI] = condVarFieldInterpolatorI.interpolate(position,cellI);		    	
			nPartInI++;
		    }
		    
	    }
	    //- Resize condVarFieldXp
	    //condVarFieldXp.resize(nPartIn);

            ii++; 
	    word& condVariableNameJ = averagingVariableList[ii];
            const dictionary& subDictCondVariableNameJ = multipleVarsConditionalAveDict_.subDict(condVariableNameJ);
            Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableNameJ << 
		              " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;

            // Parse expression and write into result folder
            expressionField condVarJ(condVariableNameJ,mesh,subDictCondVariableNameJ,false);
		
            // Read averaging volume field
            IOobject condVarHeaderJ
            (
                condVariableNameJ,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
	     	     
            volScalarField condVarFieldJ(condVarHeaderJ,mesh);
	    //- Create interpolation class of conditional variable 
	    scalarField condVarFieldXpJ(dragF().nP(),0);
	    interpolationCellPoint<scalar> condVarFieldInterpolatorJ(condVarFieldJ);
	    label nPartInJ(0);
	    for(int index=0; index < dragF().nP(); index++)
	    {
		    label cellI = dragF().cellID(index);
		    vector position = dragF().xp(index); 
		    //- Needed for parallel computation
		    if( cellI > -1 ) 
		    {
		    	//condVarFieldXp[nPartIn] = condVarFieldInterpolator.interpolate(position,cellI);
		    	condVarFieldXpJ[nPartInJ] = condVarFieldInterpolatorJ.interpolate(position,cellI);		    	
			nPartInJ++;
		    }
		    
	    }
	
            // Conditional averaging subroutine
/*             condAveUserDefine( varField, condVarFieldXp, minCondVar[ii], maxCondVar[ii], nBin[ii], 
					  condAveVar[userLoopVar+ii*nAnalyses],
				               nReal[userLoopVar+ii*nAnalyses]  );  */	
            condAveUserDefine( varField, condVarFieldXpI, minCondVar[0], maxCondVar[0], nBin[0], 
					 condVarFieldXpJ, minCondVar[1], maxCondVar[1], nBin[1], 
					  condAveVar[userLoopVar],
				               nReal[userLoopVar]  ); 	
	    //- Debugging			       
	    //Info << "In calcCondAve ii " << ii << " userLoopVar " <<  userLoopVar << " caseName " << userLoopVar << endl;			       		          
	    Info << " In writemultVarsCondAve ii " <<  userLoopVar << " caseName " << userLoopVar << endl;

        //}
}

void Foam::multipleVarsConditionalAve::write(word& userDefineName)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);	
        forAll(variableList,jj)
        {
            word& variableName = variableList[jj];
/*             forAll(averagingVariableList,ii)
            {
                word& condVariableName = averagingVariableList[ii];                
                writeOutput( mesh, outputPath_, variableName, condVariableName, minCondVar[jj+ii], maxCondVar[jj+ii], nBin[jj+ii], condAveVar[jj+ii], nReal[jj+ii], userDefineName );            
            } */
            writeOutput( mesh, outputPath_, variableName, 
	    			averagingVariableList[0], minCondVar[0], maxCondVar[0], nBin[0], 
	    			averagingVariableList[1], minCondVar[1], maxCondVar[1], nBin[1],
	    			condAveVar[jj], nReal[jj], userDefineName );            
        } 
}

void Foam::multipleVarsConditionalAve::writeUserDefine
(
	word& variableName,	
	word& userDefineName,
	label& userLoopVar	
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);		    
	//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );			    
	label nAnalyses = nTotalCase_ / nAveragingVariable_;

	//forAll(averagingVariableList,ii)
        //{
            //word& condVariableName = averagingVariableList[ii];
/*             writeOutput( mesh, outputPath_, variableName, condVariableName, minCondVar[ii], maxCondVar[ii], nBin[ii], 
					  condAveVar[userLoopVar+ii*nAnalyses],
				               nReal[userLoopVar+ii*nAnalyses],  
						               userDefineName  ); */
            writeOutput( mesh, outputPath_, variableName, 
	    				averagingVariableList[0], minCondVar[0], maxCondVar[0], nBin[0], 
	    				averagingVariableList[1], minCondVar[1], maxCondVar[1], nBin[1], 
					  condAveVar[userLoopVar],
				               nReal[userLoopVar],  
						  userDefineName  );
     	    Info << "In writemultVarsCondAve ii userLoopVar " <<  userLoopVar << " caseName " << userLoopVar << endl;
	//}
}

void Foam::multipleVarsConditionalAve::writeOutput
(   
	const fvMesh&       mesh_,
	fileName&           outputPath_,
	word&               varName_,
	word&               condVarNameI_,
	const scalar&       minCondVarI_,
	const scalar&       maxCondVarI_,
	const label&        nBinI_,
	word&               condVarNameJ_,
	const scalar&       minCondVarJ_,
	const scalar&       maxCondVarJ_,
	const label&        nBinJ_,
	scalarField&        condAveVar_,
	scalarField&        nReal_,
	word&		    userDefineName_
)
{
	// Create output folder
	fileName outputFolder;
	fileName outputFile;
	if(Pstream::master())
	{

		if(mesh_.time().processorCase())
		{
			outputFolder = mesh_.time().path()/".."/outputPath_/mesh_.time().timeName();
			if(!isDir(outputFolder)) mkDir(outputFolder);
		}else
		{
			outputFolder = mesh_.time().path()/outputPath_/mesh_.time().timeName();
			if(!isDir(outputFolder)) mkDir(outputFolder);
		}
	
		// Create file
		outputFile = "multiVarsCondAve"+varName_+condVarNameI_+condVarNameJ_+userDefineName_;
	}

	//- Total realisation	
	scalarField nTotalEspectValue = condAveVar_; 
	scalarField nTotalRealization = nReal_; 

				
	//- Summation of processors results
	reduce(nTotalEspectValue, sumOp<scalarField>());
	reduce(nTotalRealization, sumOp<scalarField>());
		
	//- Normalize in the following loop (var*nreal and parallel compt var*nreal*nreal)
	// nTotalEspectValue /= nTotalRealization;

	if(Pstream::master())
	{         
		// Bin interval size
		scalar binSizeI = ( maxCondVarI_ - minCondVarI_ ) / nBinI_;
		scalar binSizeJ = ( maxCondVarJ_ - minCondVarJ_ ) / nBinJ_;
		
		//Info << "binSizeI = " << binSizeI << " binSizeJ = " << binSizeJ << endl;
	
		OFstream outputPtr(outputFolder/outputFile);
		// Write header
		outputPtr << "#" << tab << condVarNameI_ << tab << condVarNameJ_ << tab << varName_ << tab << " Nreal" << endl;
		
		// Write into the file
		forAll(nTotalEspectValue,ii)
		{
		    int II = ii%nBinI_;
		    int JJ = floor(ii/nBinI_);
		    
		    //Info << "ii = " << ii << " II = " << binSizeJ << endl;
		    		    
   		    // To avoid dividing  by zero 
		    if ( nTotalRealization[ii]== 0 )
		    {
		      outputPtr   <<  setw(IOstream::defaultPrecision() + 6)
                		  <<  minCondVarI_ + binSizeI*(1./2.+II)  		<<  tab
                		  <<  minCondVarJ_ + binSizeJ*(1./2.+JJ)   		<<  tab				  
                		  <<  nTotalEspectValue[ii] 		   		<<  tab
                		  <<  nTotalRealization[ii]                      	<<  tab
                		  <<  endl;
		    }else
		    {
		      outputPtr   <<  setw(IOstream::defaultPrecision() + 6)
                		  <<  minCondVarI_ + binSizeI*(1./2.+II)  		<<  tab
                		  <<  minCondVarJ_ + binSizeJ*(1./2.+JJ)   		<<  tab				  
                		  <<  nTotalEspectValue[ii] / nTotalRealization[ii]  	<<  tab
                		  <<  nTotalRealization[ii]                      	<<  tab
                		  <<  endl;	    
		    }		  
		}
	
	}
}
// ************************************************************************* //
