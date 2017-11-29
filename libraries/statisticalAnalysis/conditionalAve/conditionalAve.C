/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Calculate conditionalAve

-----------------------------------------------------------------------------*/

#include "conditionalAve.H"

namespace Foam {
    defineTypeNameAndDebug(conditionalAve,0);
}

Foam::conditionalAve::conditionalAve
(
    const dictionary& dict,
    const dictionary& subDict,    
    const objectRegistry& obr,
    const label& nVariable,
    const label& nAveragingVariable,
    const label& nTotalCase,
    const bool& conditionalAveraging  
)
:
    	dict_(dict),
	obr_(obr),
        nVariable_(nVariable),
        nAveragingVariable_(nAveragingVariable),
        nTotalCase_(nTotalCase), 
	//conditionalAveDict_(dict_.subDict("conditionalAve")),
	conditionalAveDict_(subDict),
	//outputPath_(conditionalAveDict_.lookup("outputRelativePath")),	
	outputPath_(),
	variableList(nVariable,""),
	averagingVariableList(nAveragingVariable,""),
	nBin(nAveragingVariable,1),
	minCondVar(nAveragingVariable,scalar(-1.e+12)),
	maxCondVar(nAveragingVariable,scalar( 1.e+12)),
	nReal(nTotalCase,scalarField(1,scalar(0))),	
	condAveVar(nTotalCase,scalarField(1,scalar(0))),		
	readFromFile_(false),
	conditionalAveraging_(conditionalAveraging)
{
	// Init
	init();	
}

Foam::conditionalAve::~conditionalAve()
{}

void Foam::conditionalAve::init()
{	
	// Create output folder
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	if(conditionalAveraging_)
	{
		//conditionalAveDict_ = dict_.subDict("conditionalAve");
		outputPath_ = fileName(conditionalAveDict_.lookup("outputRelativePath"));
		
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

		// Read variables list and dictionaries
        	variableList = wordList(conditionalAveDict_.lookup("variableList"));
        	Info << "\nThe variable(s) \n";
        	forAll(variableList,ii)
        	{
        	    Info << tab << variableList[ii];
        	}
        	Info << "\nwill be conditional averaged by the following variable(s):" << endl;
        	averagingVariableList = wordList(conditionalAveDict_.lookup("averagingVariableList"));
        	forAll(averagingVariableList,ii)
        	{
        	    word& condVariableName = averagingVariableList[ii];
        	    Info << tab << condVariableName;
        	}
        	Info << endl;

		//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );
		label nAnalyses = nTotalCase_ / nAveragingVariable_;
		
		Info << "nAnalyses = " << nAnalyses << endl;

		label caseN(0);
        	forAll(averagingVariableList,ii)
        	{
        	    word& currentVar = averagingVariableList[ii];
        	    const dictionary& subDict = conditionalAveDict_.subDict(currentVar);
        	    nBin[ii] = readLabel(subDict.lookup(word("nBin")));

        	    minCondVar[ii] = readScalar(subDict.lookup(word("min")));
        	    maxCondVar[ii] = readScalar(subDict.lookup(word("max")));

        	    for(label nA = 0; nA < nAnalyses; nA++)
		    { 
			    //- Resize variables
                	    condAveVar[caseN].resize(nBin[ii]);
                	    nReal[caseN].resize(nBin[ii]);

			    //- Re-initialise 
			    condAveVar[caseN] = scalarField(nBin[ii],scalar(0));
			    nReal[caseN] = scalarField(nBin[ii],scalar(0));

			    caseN++;
		    }	
        	}

		Info << " " << endl;
	}
	else
	{
		Info << "Conditional averaging will NOT be performed" << endl;
	}
}

//- Condition averaging for volScalarField
void Foam::conditionalAve::condAve
(   
	 const fvMesh&         	mesh_,
	 const volScalarField& 	var_,
	 const volScalarField& 	condVar_,
	 const scalar&         	minCondVar_,
	 const scalar&         	maxCondVar_,
	 const label&          	nBin_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
)
{
	 // Bin interval size
	 scalar binSize = ( maxCondVar_ - minCondVar_ ) / nBin_;

	 // Loop over cells
	 forAll(mesh_.cells(),cellI)
	 {
	     label binI = floor(( condVar_[cellI] - minCondVar_ ) / binSize );
	     if ( binI > -1 && binI < nBin_ )
	     {
        	 condAveVar_[binI] += var_[cellI];
        	 nReal_[binI] += 1.;
	     }
	 }
}

//- Condition averaging for user define scalar field
void Foam::conditionalAve::condAveUserDefine
(   
	 const scalarField& 	var_,
	 const scalarField& 	condVar_,
	 const scalar&         	minCondVar_,
	 const scalar&         	maxCondVar_,
	 const label&          	nBin_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
)
{
	 // Bin interval size
	 scalar binSize = ( maxCondVar_ - minCondVar_ ) / nBin_;

	 // Loop over cells
	 forAll(var_,cellI)
	 {
	    label binI = floor(( condVar_[cellI] - minCondVar_ ) / binSize );
	    if ( binI > -1 && binI < nBin_ )
	    {
               condAveVar_[binI] += var_[cellI];
               nReal_[binI] += 1.;

	    //if(binI==9) 
	    //Info << " binI = " << binI << " condVar_[cellI] = " << condVar_[cellI] << " condAveVar_[binI] " << condAveVar_[binI]
	    //     << " var_[cellI] = " << var_[cellI] << " nReal_[binI] = " << nReal_[binI] <<endl;

	    }
	 }
}

void Foam::conditionalAve::calc()
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
        	const dictionary& subDictVariableName = conditionalAveDict_.subDict(variableName);
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

            forAll(averagingVariableList,ii)
            {
        	word& condVariableName = averagingVariableList[ii];
        	const dictionary& subDictCondVariableName = conditionalAveDict_.subDict(condVariableName);
        	//Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableName 
		//     << " [" << minCondVar[ii+jj] << "," << maxCondVar[ii+jj] << "] ..."  << endl;
        	Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableName 
		     << " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;	

        	// Parse expression and write into result folder
        	expressionField condVar(condVariableName,mesh,subDictCondVariableName,false);

        	// Read averaging volume field
        	IOobject condVarHeader
        	(
                    condVariableName,
                    //runTime.timeName(),
		    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ
        	 );
        	volScalarField condVarField(condVarHeader,mesh);

        	// Conditional averaging subroutine
        	volScalarField varField = *ptrVar;
        	condAve( mesh, varField, condVarField, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[ii+jj], nReal[ii+jj] );
		Info << " calc routine (all) condAveVar[ii+jj] = " << ii+jj << endl;

            }
	}
}

void Foam::conditionalAve::calc
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
            const dictionary& subDictVariableName = conditionalAveDict_.subDict(variableName);
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

        forAll(averagingVariableListCoarsening,ii)
        {
            word& condVariableName = averagingVariableListCoarsening[ii];
            const dictionary& subDictCondVariableName = conditionalAveDict_.subDict(condVariableName);
            //Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableName 
	    //     << " [" << minCondVar[ii+jj] << "," << maxCondVar[ii+jj] << "] ..."  << endl;
            Info <<  "\nVariable " << variableName << " conditionally averaging by " << condVariableName 
		 << " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;	

            // Parse expression and write into result folder
            expressionField condVar(condVariableName,mesh,subDictCondVariableName,false);

            // Read averaging volume field
            IOobject condVarHeader
            (
                condVariableName,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
            volScalarField condVarField(condVarHeader,mesh);

            // Conditional averaging subroutine
            volScalarField varField = *ptrVar;
            condAve( mesh, varField, condVarField, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[ii+jj], nReal[ii+jj] );
	    Info << " In calc(.,.,.,.), condAveVar[ii+jj] = " << ii+jj << " where ii = " << ii << " jj = " << jj << endl;
	    
	    // Write into file
            writeOutput( mesh, outputPath_, variableName, condVariableName, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[jj+ii], nReal[jj+ii], userDefineName ); 

        }

}

void Foam::conditionalAve::calcUserDefine
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

	forAll(averagingVariableList,ii)
        {
            word& condVariableName = averagingVariableList[ii];
            const dictionary& subDictCondVariableName = conditionalAveDict_.subDict(condVariableName);
            Info <<  "\n-->Variable " << variableName << " conditionally averaging by " << condVariableName << 
		              " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;

            // Parse expression and write into result folder
            expressionField condVar(condVariableName,mesh,subDictCondVariableName,false);
		
            // Read averaging volume field
            IOobject condVarHeader
            (
                condVariableName,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
	     	     
            volScalarField condVarField(condVarHeader,mesh);
	    //- Create interpolation class of conditional variable 
	    scalarField condVarFieldXp(dragF().nP(),0);
	    Info << " var.size = " << dragF().nP() << endl;
	    interpolationCellPoint<scalar> condVarFieldInterpolator(condVarField);
	    label nPartIn(0);
	    for(int index=0; index < dragF().nP(); index++)
	    {
		    label cellI = dragF().cellID(index);
		    vector position = dragF().xp(index); 
		    //- Needed for parallel computation
		    if( cellI > -1 ) 
		    {
		    	condVarFieldXp[nPartIn] = condVarFieldInterpolator.interpolate(position,cellI);
		    	nPartIn++;
		    }
		    
	    }
	    
	    //- Resize condVarFieldXp
	    //condVarFieldXp.resize(nPartIn);
		
            // Conditional averaging subroutine
            condAveUserDefine( varField, condVarFieldXp, minCondVar[ii], maxCondVar[ii], nBin[ii], 
					  condAveVar[userLoopVar+ii*nAnalyses],
				               nReal[userLoopVar+ii*nAnalyses]  ); 	
	    //- Debugging			       
	    Info << "In calcCondAve(calcUserDefine) ii " << ii << " userLoopVar " <<  userLoopVar << " caseName " << userLoopVar+ii*nAnalyses << endl;			       		          

        }
}

void Foam::conditionalAve::calcUserDefine
(
	scalarField& varField,
	autoPtr<dragModelPost>& dragF,
	word& variableName,
	label& userLoopVar,
        List<List<scalar> >& particleWeights_,
        labelListList& neighboringCellIDs_,
	int& maxCellPerPart_
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );			    
	label nAnalyses = nTotalCase_ / nAveragingVariable_;	

	forAll(averagingVariableList,ii)
        {
            word& condVariableName = averagingVariableList[ii];
            const dictionary& subDictCondVariableName = conditionalAveDict_.subDict(condVariableName);
            Info <<  "\n-->Variable " << variableName << " conditionally averaging by " << condVariableName << 
		              " [" << minCondVar[ii] << "," << maxCondVar[ii] << "] ..."  << endl;

            // Parse expression and write into result folder
            expressionField condVar(condVariableName,mesh,subDictCondVariableName,false);
		
            // Read averaging volume field
            IOobject condVarHeader
            (
                condVariableName,
                //runTime.timeName(),
		mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
             );
	     	     
            volScalarField condVarField(condVarHeader,mesh);
	    //- Create interpolation class of conditional variable 
	    scalarField condVarFieldXp(dragF().nP(),0);
	    Info << " var.size = " << dragF().nP() << endl;
	    label nPartIn(0);
	    scalar weightP(0); 
	    for(int index=0; index < dragF().nP(); index++)
	    {
		    label cellI = dragF().cellID(index);
		    //- Needed for parallel computation
		    if( cellI > -1 ) 
		    {
  			// Centering cell
			weightP = particleWeights_[index][0];
			condVarFieldXp[nPartIn] = weightP*condVarField[cellI];

			// Neighboring cells    
			for(int subCell=1 ; subCell < maxCellPerPart_; subCell++)
			{
        		    label cellJ = neighboringCellIDs_[index][subCell];
        		    if (cellJ > -1)
        		    {
        			weightP = particleWeights_[index][subCell];
				condVarFieldXp[nPartIn] += weightP*condVarField[cellI];
        		    }
			}

		    	nPartIn++;
		    }		    
	    }
	    
	    //- Resize condVarFieldXp
	    //condVarFieldXp.resize(nPartIn);
		
            // Conditional averaging subroutine
            condAveUserDefine( varField, condVarFieldXp, minCondVar[ii], maxCondVar[ii], nBin[ii], 
					  condAveVar[userLoopVar+ii*nAnalyses],
				               nReal[userLoopVar+ii*nAnalyses]  ); 	
	    //- Debugging			       
	    Info << "In calcCondAve(calcUserDefine) ii " << ii << " userLoopVar " <<  userLoopVar << " caseName " << userLoopVar+ii*nAnalyses << endl;			       		          

        }
}

void Foam::conditionalAve::write(word& userDefineName)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);	
        forAll(variableList,jj)
        {
            word& variableName = variableList[jj];
            forAll(averagingVariableList,ii)
            {
                word& condVariableName = averagingVariableList[ii];                
                writeOutput( mesh, outputPath_, variableName, condVariableName, minCondVar[ii], maxCondVar[ii], nBin[ii], condAveVar[jj+ii], nReal[jj+ii], userDefineName ); 
		Info << " Output, condAveVar[ii+jj] = " << ii+jj << endl;           
            }
        } 
}

void Foam::conditionalAve::writeUserDefine
(
	word& variableName,	
	word& userDefineName,
	label& userLoopVar	
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);		    
	//label nAnalyses = nTotalCase_ / ( nVariable_ * nAveragingVariable_ );			    
	label nAnalyses = nTotalCase_ / nAveragingVariable_;

	forAll(averagingVariableList,ii)
        {
            word& condVariableName = averagingVariableList[ii];
            writeOutput( mesh, outputPath_, variableName, condVariableName, minCondVar[ii], maxCondVar[ii], nBin[ii], 
					  condAveVar[userLoopVar+ii*nAnalyses],
				               nReal[userLoopVar+ii*nAnalyses],  
						               userDefineName  );
     	    Info << "In writeCondAve(writeUserDefine) ii " << ii << " userLoopVar " <<  userLoopVar << " caseName " << userLoopVar+ii*nAnalyses << endl;
	}
}

void Foam::conditionalAve::writeOutput
(   
	const fvMesh&       mesh_,
	fileName&           outputPath_,
	word&               varName_,
	word&               condVarName_,
	const scalar&       minCondVar_,
	const scalar&       maxCondVar_,
	const label&        nBin_,
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
		outputFile = "condAve"+varName_+condVarName_+userDefineName_;
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
		scalar binSize = ( maxCondVar_ - minCondVar_ ) / nBin_;
	
		OFstream outputPtr(outputFolder/outputFile);
		// Write header
		outputPtr << "#" << tab << condVarName_ << tab << varName_ << tab << " Nreal" << endl;
		
		// Write into the file
		forAll(nTotalEspectValue,ii)
		{
		    // To avoid dividing  by zero
		    if ( nTotalRealization[ii]== 0 )
		    {
		      outputPtr   <<  setw(IOstream::defaultPrecision() + 6)
                		  <<  minCondVar_ + binSize*(1./2.+ii)      		<<  tab
                		  <<  nTotalEspectValue[ii] 		   		<<  tab
                		  <<  nTotalRealization[ii]                      	<<  tab
                		  <<  endl;
		    }else
		    {
		      outputPtr   <<  setw(IOstream::defaultPrecision() + 6)
                		  <<  minCondVar_ + binSize*(1./2.+ii)              	<<  tab
                		  <<  nTotalEspectValue[ii] / nTotalRealization[ii]  	<<  tab
                		  <<  nTotalRealization[ii]                      	<<  tab
                		  <<  endl;	    
		    }		  
		}
	
	}
}
// ************************************************************************* //
