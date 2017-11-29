
    //- Define test filter size
    int testFilterWidth(2*((filterWidth[fWidth]-1)/2+1)+1);
    //int testFilterWidth(3*((filterWidth[fWidth]-1)/2+1)+1); 
    //- Create name for test filter width			  
    char charName[100];
    sprintf(charName, "%dX%dX%d",testFilterWidth,testFilterWidth,testFilterWidth);
    word testFilterWidthName(charName);
    Info << tab << "Test filter size = " << testFilterWidthName << endl; 
    
    //- Creating re-filtered variables at test filter level
    
    // hatBaralpp
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
    
    // hatDiffBarAlppHatBarAlppSqr
    volScalarField hatDiffBarAlppHatBarAlppSqr
    (
	IOobject
	(
            "hatDiffBarAlppHatBarAlppSqr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	(barAlpp-hatBarAlpp)*(barAlpp-hatBarAlpp)
    );  
    hatDiffBarAlppHatBarAlppSqr = convKernel.coarseningSca("hatDiffBarAlppHatBarAlppSqr",testFilterWidth,(barAlpp-hatBarAlpp)*(barAlpp-hatBarAlpp));
    hatDiffBarAlppHatBarAlppSqr.write();  

    // hatBarAlppBarAlpp
    volScalarField hatBarAlppBarAlpp
    (
	IOobject
	(
            "hatBarAlppBarAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	barAlpp*barAlpp
    );  
    hatBarAlppBarAlpp = convKernel.coarseningSca("hatBarAlppBarAlpp",testFilterWidth,barAlpp*barAlpp);
    hatBarAlppBarAlpp.write();  
    
    // hatHatBaralpp
    volScalarField hatHatBaralpp
    (
	IOobject
	(
            "hatHatBaralpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	hatBarAlpp
    );  
    hatHatBaralpp = convKernel.coarseningSca("hatHatBaralpp",testFilterWidth,hatBarAlpp);
    hatHatBaralpp.write();    
 
    // hatflucAlppFlucAlpp
    volScalarField hatFlucAlppFlucAlpp
    (
	IOobject
	(
            "hatFlucAlppFlucAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	flucAlppFlucAlpp
    );
    hatFlucAlppFlucAlpp = convKernel.coarseningSca("hatflucAlppFlucAlpp",testFilterWidth,flucAlppFlucAlpp);
    hatFlucAlppFlucAlpp.write();
            
    // Dynamic procedure- Leonard stresses
    volScalarField Lvar
    (
	IOobject
	(
            "Lvar"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	hatDiffBarAlppHatBarAlppSqr-hatFlucAlppFlucAlpp
    );
    Lvar.write(); 

    // Dynamic procedure- Leonard stresses
    volScalarField LvarApprox
    (
	IOobject
	(
            "LvarApprox"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	hatBarAlppBarAlpp-hatBarAlpp*hatBarAlpp
    );
    LvarApprox.write(); 
    
    // Proposed function form + using swak4Foam
    word& variableName = dynamicAdjustmentProposedFunctionList[0];
    volScalarField* ptrVar;
    // Look up for variable
    if (mesh.foundObject<volScalarField>(variableName))
    {
        ptrVar = new volScalarField(mesh.lookupObject<volScalarField>(variableName));
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
        ptrVar = new volScalarField(varHeader,mesh);
    } 

    volScalarField hatProposedFunction
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
    hatProposedFunction = convKernel.coarseningSca(variableName,testFilterWidth,*ptrVar);
    hatProposedFunction.write();
		
    // Using proposed function form with resolved variables + using swak4Foam
    word& variableNameResolved = dynamicAdjustmentResolvedFunctionList[0];
    volScalarField* ptrVarResolved;
    // Look up for variable
    if (mesh.foundObject<volScalarField>(variableNameResolved))
    {
        ptrVarResolved = new volScalarField(mesh.lookupObject<volScalarField>(variableNameResolved));
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
        ptrVarResolved = new volScalarField(varHeaderResolved,mesh);
    } 
     
    volScalarField resolvedProposedFunction
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
    volScalarField Mvar
    (
	IOobject
	(
            "Mvar"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	resolvedProposedFunction-hatProposedFunction
    );
    Mvar.write();

    // Dynamic procedure - Model constant Kbetabeta
    volScalarField Kbetabeta
    (
	IOobject
	(
            "Kbetabeta"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	Lvar/(Mvar+scalar(1.e-06))
    );
    Kbetabeta.write();

    // Dynamic procedure - Model constant KbetabetaApprox
    volScalarField KbetabetaApprox
    (
	IOobject
	(
            "KbetabetaApprox"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	LvarApprox/(Mvar+scalar(1.e-06))
    );
    KbetabetaApprox.write();

    //- Predicted scalar variance 
    volScalarField predictedKbetabetaBarFlucAlppFlucAlpp
    (
	IOobject
	(
            "predictedKbetabetaBarFlucAlppFlucAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	Kbetabeta*(*ptrVar)
    );
    predictedKbetabetaBarFlucAlppFlucAlpp.write();

    //- Predicted scalar variance 
    volScalarField predictedKbetabetaApproxPDFScalarVar
    (
	IOobject
	(
            "predictedKbetabetaApproxPDFScalarVar"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	KbetabetaApprox*(*ptrVar)
    );
    predictedKbetabetaApproxPDFScalarVar.write();
    
    scalar domainKbetabeta = fvc::domainIntegrate(Lvar*Mvar).value()/fvc::domainIntegrate(Mvar*Mvar).value();
    Info << tab << "Domain-average Kbetabeta = " << domainKbetabeta << endl;

    scalar domainKbetabetaApprox = fvc::domainIntegrate(LvarApprox*Mvar).value()/fvc::domainIntegrate(Mvar*Mvar).value();
    Info << tab << "Domain-average KbetabetaApprox = " << domainKbetabetaApprox << endl;

    volScalarField predictedDomainKbetabetaBarFlucAlppFlucAlpp
    (
	IOobject
	(
            "predictedDomainKbetabetaBarFlucAlppFlucAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	domainKbetabeta*(*ptrVar)
    );
    predictedDomainKbetabetaBarFlucAlppFlucAlpp.write();

    volScalarField predictedDomainKbetabetaApproxPDFScalarVar
    (
	IOobject
	(
            "predictedDomainKbetabetaApproxPDFScalarVar"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	domainKbetabetaApprox*(*ptrVar)
    );
    predictedDomainKbetabetaApproxPDFScalarVar.write();
    
    //- Calculate ApproxPDFScalarVarAtTestScale
    
    //- hatBarAlppAlpp
    volScalarField hatBarAlppAlpp
    (
	IOobject
	(
            "hatBarAlppAlpp"+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
	),
	barAlppAlpp
    );  
    hatBarAlppAlpp = convKernel.coarseningSca("hatBarAlppAlpp",testFilterWidth,hatBarAlppAlpp); 
    
    //- ApproxPDFScalarVarAtTestScale 
    volScalarField approxPDFScalarVarAtTestScale
    (
	IOobject
	(
            "approxPDFScalarVarAtTestScale"+testFilterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	hatBarAlppAlpp-hatBarAlpp*hatBarAlpp
    );
    approxPDFScalarVarAtTestScale.write();  
