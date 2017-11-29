		
// Filtering gas volume fraction
//volScalarField barAlpf(convKernel.coarseningSca("barAlpf",filterWidth[fWidth],voidfraction));
//barAlpf = convKernel.coarseningSca("barAlpf",filterWidth[fWidth],voidfraction);
barAlpf = convKernel.coarseningSca("barAlpf",filterWidth[fWidth],mappedVoidfraction);
barAlpp = 1.-barAlpf;
//barAlppAlpp = convKernel.coarseningSca("barAlppAlpp",filterWidth[fWidth],(1-voidfraction)*(1-voidfraction));
barAlppAlpp = convKernel.coarseningSca("barAlppAlpp",filterWidth[fWidth],(1-mappedVoidfraction)*(1-mappedVoidfraction));

// Filtering gas velocity
//volVectorField tildeUf(convKernel.coarseningVecFavre("tildeUf",filterWidth[fWidth],voidfraction,U));
//tildeUf = convKernel.coarseningVecFavre("tildeUf",filterWidth[fWidth],voidfraction,U);
//barAlppUf = convKernel.coarseningVec("barAlppUf",filterWidth[fWidth],(1-voidfraction)*U);
tildeUf = convKernel.coarseningVecFavre("tildeUf",filterWidth[fWidth],mappedVoidfraction,U);
barAlppUf = convKernel.coarseningVec("barAlppUf",filterWidth[fWidth],(1-mappedVoidfraction)*U);

// Filtering gas pressure gradient
//volVectorField tildeGradPg(convKernel.coarseningVec("tildeGradPg",filterWidth[fWidth],gradPg));
tildeGradPg = convKernel.coarseningVec("tildeGradPg",filterWidth[fWidth],gradPg);

// Filtering solid velocity
//tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-voidfraction,Us);
tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-mappedVoidfraction,mappedUs);

// Scalar variance
//tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidth[fWidth],(1-voidfraction)*(U-Us));
tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidth[fWidth],(1-mappedVoidfraction)*(U-Us));

// Write gradBarAlpp
volVectorField gradBarAlpp
(
    IOobject
    (
        "gradBarAlpp"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, //MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(barAlpp)
);  
gradBarAlpp.write();

// Create filtered Eulerian force
//Info << " Dimensions coarsening fluid coarsening" << barForceEulerian.dimensions() << endl;	
barForceEulerian = convKernel.coarseningVec("barForceEulerian",filterWidth[fWidth],dragFpEulerian);

// Scalar variance (LES modeling)
flucAlppFlucAlpp = ( ( (1.-mappedVoidfraction) - barAlpp ) * ( (1.-mappedVoidfraction) - barAlpp ) );

// Apply filtering
barFlucAlppFlucAlpp = convKernel.coarseningSca("barFlucAlppFlucAlpp",filterWidth[fWidth],flucAlppFlucAlpp);

//- barAlppFilter
volScalarField barAlppFilter
(
    IOobject
    (
        "barAlpp"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barAlpp
);
barAlppFilter.write();  

//- ApproxPDFScalarVar
volScalarField approxPDFScalarVar
(
    IOobject
    (
        "approxPDFScalarVar"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barAlppAlpp-barAlpp*barAlpp
);
approxPDFScalarVar.write();    

// Dynamic procedure
if(dynamicAdjustment)
{
    //- Define test filter size
    //int testFilterWidth(2*((filterWidth[fWidth]-1)/2+1)+1);
    int testFilterWidth(3*((filterWidth[fWidth]-1)/2+1)+1); 
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
            "Lvar"+testFilterWidthName,
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
            "LvarApprox"+testFilterWidthName,
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
            "Mvar"+testFilterWidthName,
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
            "Kbetabeta"+testFilterWidthName,
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
            "KbetabetaApprox"+testFilterWidthName,
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
            "predictedKbetabetaBarFlucAlppFlucAlpp"+testFilterWidthName,
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
            "predictedKbetabetaApproxPDFScalarVar"+testFilterWidthName,
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
            "predictedDomainKbetabetaBarFlucAlppFlucAlpp"+testFilterWidthName,
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
            "predictedDomainKbetabetaApproxPDFScalarVar"+testFilterWidthName,
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
}
