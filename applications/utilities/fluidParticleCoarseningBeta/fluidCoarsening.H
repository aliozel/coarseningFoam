		
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

// Filtering alpp*relativeVel
//tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidth[fWidth],(1-voidfraction)*(U-Us));
tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidth[fWidth],(1-mappedVoidfraction)*(U-mappedUs));

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

// Added by Yi to get a new marker of mesoshear rate
volTensorField SMeso("SMeso",fvc::grad(tildeUs) + fvc::grad(tildeUs)().T()-2./3.*tr(fvc::grad(tildeUs))*tensor(I));
volScalarField mesoShearRate
(
    IOobject
    (
        "mesoShearRate", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRate =  sqrt( 2.0 * (SMeso/2.0) && (SMeso/2.0));
mesoShearRate.write();  

//- Vdrift
volVectorField barAlppVdrift
(
    IOobject
    (
        "barAlppVdrift"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tildeAlppUr-barAlpp*(tildeUf-tildeUs)
);
barAlppVdrift.write(); 

//- Delta^2*|Grad(alpp)|^2 
volScalarField delta2NormGradAlpp
(
    IOobject
    (
        "delta2NormGradAlpp"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    magSqr(fvc::grad(barAlpp))
);
delta2NormGradAlpp.internalField() *= (filterWidth[fWidth])*(filterWidth[fWidth])*pow(mesh.V(), 2.0/3.0);
delta2NormGradAlpp.write(); 

//- Delta^2*Grad(alpp)Grad(tildeUf)
volVectorField delta2GradAlppGradTildeUf
(
    IOobject
    (
        "delta2GradAlppGradTildeUf"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   (fvc::grad(barAlpp) & fvc::grad(tildeUf))
);
delta2GradAlppGradTildeUf.internalField() *= (filterWidth[fWidth])*(filterWidth[fWidth])*pow(mesh.V(), 2.0/3.0);
delta2GradAlppGradTildeUf.write(); 


// Reynolds stresses
volScalarField tildeAlppUsxUsx = convKernel.coarseningSca("tildeAlppUsxUsx",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(0));
volScalarField tildeAlppUsyUsy = convKernel.coarseningSca("tildeAlppUsyUsy",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(1)*mappedUs.component(1));
volScalarField tildeAlppUszUsz = convKernel.coarseningSca("tildeAlppUszUsz",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(2)*mappedUs.component(2));

volScalarField tildeAlppUsxUsy = convKernel.coarseningSca("tildeAlppUsxUsy",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(1));
volScalarField tildeAlppUsxUsz = convKernel.coarseningSca("tildeAlppUsxUsz",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(2));
volScalarField tildeAlppUsyUsz = convKernel.coarseningSca("tildeAlppUsyUsz",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(1)*mappedUs.component(2));

volSymmTensorField diffTildeeAlppUsUs
(
    IOobject
    (
        "diffTildeeAlppUsUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
);


forAll(mesh.cells(),cellI)
{
	diffTildeeAlppUsUs[cellI] = symmTensor
					(		
			  		  tildeAlppUsxUsx[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][0] ),
			  		  tildeAlppUsxUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][1] ),
			  		  tildeAlppUsxUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][2] ),
			  		  tildeAlppUsyUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][1] ),
			  		  tildeAlppUsyUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][2] ),
			  		  tildeAlppUszUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][2]*tildeUs[cellI][2] )
					);
}



//- Meso-scale pressure & viscosity
// Particle pressure
volScalarField PpMesoRho
(
    IOobject
    (
        "PpMesoRho",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);
PpMesoRho = 1./3. * tr( diffTildeeAlppUsUs );
// Write into the results folder
PpMesoRho.write();
