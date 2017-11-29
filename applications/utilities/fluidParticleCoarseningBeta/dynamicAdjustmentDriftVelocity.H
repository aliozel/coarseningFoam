    //- Define test filter size
    //int testFilterWidth(2*((filterWidth[fWidth]-1)/2+1)+1);
    //- Create name for test filter width			  
    char charName[100];
    sprintf(charName, "%dX%dX%d",testFilterWidth,testFilterWidth,testFilterWidth);
    word testFilterWidthName(charName);
    Info << tab << "Test filter size = " << testFilterWidthName << endl; 
    
    //- Creating re-filtered variables at test filter level
    
    // hatBarAlpp
    volScalarField hatBarAlpp
    (
	IOobject
	(
            "hatBarAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	barAlpp
    );  
    hatBarAlpp = convKernel.coarseningSca("hatBarAlpp",testFilterWidth,barAlpp);
    hatBarAlpp.write();
    
    // hatTildeUf
    volVectorField hatTildeUf
    (
	IOobject
	(
            "hatTildeUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	tildeUf
    );  
    hatTildeUf = convKernel.coarseningVecFavre("hatTildeUf",testFilterWidth,barAlpf,tildeUf);
    hatTildeUf.write();
            
    // hatTildeUs
    volVectorField hatTildeUs
    (
	IOobject
	(
            "hatTildeUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	tildeUs
    );  
    hatTildeUs = convKernel.coarseningVecFavre("hatTildeUs",testFilterWidth,barAlpp,tildeUs);
    hatTildeUs.write();

    // hatBaralppTildeUfTildeUs
    volVectorField hatBaralppTildeUfTildeUs
    (
	IOobject
	(
            "hatBaralppTildeUfTildeUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	barAlpp*(tildeUf-tildeUs)
    );  
    hatBaralppTildeUfTildeUs = convKernel.coarseningVec("hatBaralppTildeUfTildeUs",testFilterWidth,hatBaralppTildeUfTildeUs);
    hatBaralppTildeUfTildeUs.write();

    // Dynamic procedure - Lvar
    volVectorField Lvar
    (
	IOobject
	(
            "Lvar"+filterWidthName, //+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	hatBaralppTildeUfTildeUs-hatBarAlpp*(hatTildeUf-hatTildeUs)
    );
    Lvar.write();        

    Info << tab << "Re-filtered Eulerian drag model class was created" << endl;
    // Create Eulerian re-filtered drag model class
    autoPtr<dragModelPost> refilteredDragFEulerian
    (
    	    dragModelPost::New(postProcessingDict,mesh)
    );	
    // Calculate Eulerian drag force
    Info << tab << "Calculating re-filtered Eulerian drag force..." << endl;
    refilteredDragFEulerian().setForceEulerian(hatBarAlpp,hatTildeUf,radii[me],hatTildeUs);
    
    volScalarField filteredTaup
    (
	IOobject
	(
            "filteredTaup", //+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("0", dimensionSet(0, 0, -1, 0, 0), 0)
    );
    forAll(mesh.cells(),cellI)
    {	    
    	filteredTaup[cellI] = refilteredDragFEulerian().taupEulerian(cellI);
    }
    filteredTaup.write();    

    volScalarField filteredTaupFilter
    (
	IOobject
	(
            "filteredTaup"+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	filteredTaup
    );
    filteredTaupFilter.write();		    
    
    // Proposed function form + using swak4Foam
    word& variableName = dynamicAdjustmentProposedFunctionList[0];
    volVectorField* ptrVar;
    // Look up for variable
    if (mesh.foundObject<volVectorField>(variableName))
    {
        ptrVar = new volVectorField(mesh.lookupObject<volVectorField>(variableName));
    }else
    {
        // Look up for expression
        const dictionary& subDictVariableName = dynamicAdjustmentDict.subDict(variableName);
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
        ptrVar = new volVectorField(varHeader,mesh);
    } 

    volVectorField proposedFunction
    (
	IOobject
	(
            variableName+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	*ptrVar
    ); 
    proposedFunction.write();

    volVectorField hatProposedFunction
    (
	IOobject
	(
            "hat"+variableName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	*ptrVar
    ); 
    hatProposedFunction = convKernel.coarseningVec(variableName,testFilterWidth,*ptrVar);
    hatProposedFunction.write();
		
    // Using proposed function form with resolved variables + using swak4Foam
    word& variableNameResolved = dynamicAdjustmentResolvedFunctionList[0];
    volVectorField* ptrVarResolved;
    // Look up for variable
    if (mesh.foundObject<volVectorField>(variableNameResolved))
    {
        ptrVarResolved = new volVectorField(mesh.lookupObject<volVectorField>(variableNameResolved));
    }else
    {
        // Look up for expression
        const dictionary& subDictVariableNameResolved = dynamicAdjustmentDict.subDict(variableNameResolved);
        // Parse expression and write into result folder
        expressionField varResolved(variableNameResolved,mesh,subDictVariableNameResolved,false);
        // Read volume field
        IOobject varHeaderResolved
        (
            variableNameResolved,
            //runTime.timeName(),
	    mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );
        ptrVarResolved = new volVectorField(varHeaderResolved,mesh);
    } 
     
    volVectorField resolvedProposedFunction
    (
	IOobject
	(
            variableNameResolved,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	*ptrVarResolved
    ); 
    resolvedProposedFunction.write();
    
    // Dynamic procedure - Model part
    volVectorField Mvar
    (
	IOobject
	(
            "Mvar"+filterWidthName, //+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	resolvedProposedFunction-hatProposedFunction
    );
    Mvar.write();

    // Dynamic procedure - Model constant Kbetabeta
    volVectorField Kbetabeta
    (
	IOobject
	(
            "Kbetabeta"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("0", dimensionSet(0, 0, -1, 0, 0), vector(0,0,0))
    );
    forAll(mesh.cells(),cellI)
    {	    
       if( Mvar[cellI][0] != 0. ) 
       	Kbetabeta[cellI][0] = Lvar[cellI][0] / Mvar[cellI][0]; 
       else
       	Kbetabeta[cellI][0] = 0.; 
        
       if( Mvar[cellI][1] != 0. ) 
	Kbetabeta[cellI][1] = Lvar[cellI][1] / Mvar[cellI][1]; 
       else
       	Kbetabeta[cellI][1] = 0.; 
       
       if( Mvar[cellI][2] != 0. )  	  
        Kbetabeta[cellI][2] = Lvar[cellI][2] / Mvar[cellI][2];
       else
        Kbetabeta[cellI][2] = 0.;
       	  
    }
    Kbetabeta.write();

    //- Predicted drift velocity
    volVectorField predictedKbetabetaDriftVelocity
    (
	IOobject
	(
            "predictedKbetabetaDriftVelocity"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );
    forAll(mesh.cells(),cellI)
    {
       predictedKbetabetaDriftVelocity[cellI][0] = Kbetabeta[cellI][0] * (*ptrVar)[cellI][0]; 
       predictedKbetabetaDriftVelocity[cellI][1] = Kbetabeta[cellI][1] * (*ptrVar)[cellI][1]; 
       predictedKbetabetaDriftVelocity[cellI][2] = Kbetabeta[cellI][2] * (*ptrVar)[cellI][2];     
    }
    predictedKbetabetaDriftVelocity.write();
    
    scalar LvarXMvarX(0);
    scalar LvarYMvarY(0);
    scalar LvarZMvarZ(0);
    scalar MvarXMvarX(0);
    scalar MvarYMvarY(0);
    scalar MvarZMvarZ(0);    

    forAll(mesh.cells(),cellI)
    {    
       LvarXMvarX += ( Lvar[cellI][0] * Mvar[cellI][0] );
       MvarXMvarX += ( Mvar[cellI][0] * Mvar[cellI][0] );
       LvarYMvarY += ( Lvar[cellI][1] * Mvar[cellI][1] );
       MvarYMvarY += ( Mvar[cellI][1] * Mvar[cellI][1] );
       LvarZMvarZ += ( Lvar[cellI][2] * Mvar[cellI][2] );
       MvarZMvarZ += ( Mvar[cellI][2] * Mvar[cellI][2] );
    }
    
    scalar domainKbetabetaX(LvarXMvarX/MvarXMvarX);
    scalar domainKbetabetaY(LvarYMvarY/MvarYMvarY);
    scalar domainKbetabetaZ(LvarZMvarZ/MvarZMvarZ);
         			     			     
    Info << tab << "Domain-average KbetabetaX = " << domainKbetabetaX << endl;
    Info << tab << "Domain-average KbetabetaY = " << domainKbetabetaY << endl;
    Info << tab << "Domain-average KbetabetaZ = " << domainKbetabetaZ << endl;

    volVectorField predictedDomainKbetabetaDriftVelocity
    (
	IOobject
	(
            "predictedDomainKbetabetaDriftVelocity"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );
    forAll(mesh.cells(),cellI)
    {
       predictedDomainKbetabetaDriftVelocity[cellI][0] = domainKbetabetaX * (*ptrVar)[cellI][0]; 
       predictedDomainKbetabetaDriftVelocity[cellI][1] = domainKbetabetaY * (*ptrVar)[cellI][1]; 
       predictedDomainKbetabetaDriftVelocity[cellI][2] = domainKbetabetaZ * (*ptrVar)[cellI][2]; 
    }
    predictedDomainKbetabetaDriftVelocity.write();

