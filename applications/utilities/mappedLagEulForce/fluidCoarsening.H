		
// Filtering gas volume fraction
//volScalarField barAlpf(convKernel.coarseningSca("barAlpf",filterWidth[fWidth],voidfraction));
barAlpf = convKernel.coarseningSca("barAlpf",filterWidth[fWidth],voidfraction);
barAlpp = 1.-barAlpf;
barAlppAlpp = convKernel.coarseningSca("barAlppAlpp",filterWidth[fWidth],(1-voidfraction)*(1-voidfraction));

// Filtering gas velocity
//volVectorField tildeUf(convKernel.coarseningVecFavre("tildeUf",filterWidth[fWidth],voidfraction,U));
tildeUf = convKernel.coarseningVecFavre("tildeUf",filterWidth[fWidth],voidfraction,U);
barAlppUf = convKernel.coarseningVec("barAlppUf",filterWidth[fWidth],(1-voidfraction)*U);

// Filtering gas pressure gradient
//volVectorField tildeGradPg(convKernel.coarseningVec("tildeGradPg",filterWidth[fWidth],gradPg));
tildeGradPg = convKernel.coarseningVec("tildeGradPg",filterWidth[fWidth],gradPg);

// Filtering solid velocity
tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-voidfraction,Us);

// Scalar variance
tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidth[fWidth],(1-voidfraction)*(U-Us));

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
