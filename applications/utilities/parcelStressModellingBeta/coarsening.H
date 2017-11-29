// Filtering solid volume fraction
//barAlpp = convKernel.coarseningSca("barAlpp",filterWidth[fWidth],1-voidfraction);
barAlpp = convKernel.coarseningSca("barAlpp",filterWidth[fWidth],1-mappedVoidfraction);

//Added by Yile to have a new marker for scalar variance

barAlppAlpp = convKernel.coarseningSca("barAlppAlpp",filterWidth[fWidth],(1-mappedVoidfraction)*(1-mappedVoidfraction));
volScalarField barBarAlpp = convKernel.coarseningSca("barBarAlpp",filterWidth[fWidth]*3,barAlpp);

// Filtering solid velocity
//tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-voidfraction,Us);
tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-mappedVoidfraction,mappedUs);
tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidth[fWidth]*3,barAlpp,tildeUs);

// Filtering gas velocit
tildeUf = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],mappedVoidfraction,U);

// Reynolds stresses
//volScalarField tildeAlppUsxUsx = convKernel.coarseningScaFavre("tildeAlppUsxUsx",filterWidth[fWidth],(1-voidfraction),Us.component(0)*Us.component(0));
//volScalarField tildeAlppUsyUsy = convKernel.coarseningScaFavre("tildeAlppUsyUsy",filterWidth[fWidth],(1-voidfraction),Us.component(1)*Us.component(1));
//volScalarField tildeAlppUszUsz = convKernel.coarseningScaFavre("tildeAlppUszUsz",filterWidth[fWidth],(1-voidfraction),Us.component(2)*Us.component(2));
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
	diffTildeeAlppUsUs[cellI] =   rhoParticle
				    * symmTensor
					(		
			  		  tildeAlppUsxUsx[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][0] ),
			  		  tildeAlppUsxUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][1] ),
			  		  tildeAlppUsxUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][2] ),
			  		  tildeAlppUsyUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][1] ),
			  		  tildeAlppUsyUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][2] ),
			  		  tildeAlppUszUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][2]*tildeUs[cellI][2] )
					);
}


diffTildeeAlppUsUs.write();

//- Meso-scale pressure & viscosity
// Particle pressure
volScalarField PpMeso
(
    IOobject
    (
        "PpMeso",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);
PpMeso = 1./3. * tr( diffTildeeAlppUsUs );
// Write into the results folder
PpMeso.write();

volScalarField PpMesoXX
(
    IOobject
    (
        "PpMesoXX",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

PpMesoXX=diffTildeeAlppUsUs.component(symmTensor::XX);
PpMesoXX.write();

volScalarField PpMesoYY
(
    IOobject
    (
        "PpMesoYY",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );
PpMesoYY=diffTildeeAlppUsUs.component(symmTensor::YY);
PpMesoYY.write();

volScalarField PpMesoZZ
(
    IOobject
    (
        "PpMesoZZ",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );
PpMesoZZ=diffTildeeAlppUsUs.component(symmTensor::ZZ);
PpMesoZZ.write();

volScalarField PpMesoXY
(
    IOobject
    (
        "PpMesoXY",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );
PpMesoXY=diffTildeeAlppUsUs.component(symmTensor::XY);
PpMesoXY.write();

volScalarField PpMesoXZ
(
    IOobject
    (
        "PpMesoXZ",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );
PpMesoXZ=diffTildeeAlppUsUs.component(symmTensor::XZ);
PpMesoXZ.write();

volScalarField PpMesoYZ
(
    IOobject
    (
        "PpMesoYZ",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );
PpMesoYZ=diffTildeeAlppUsUs.component(symmTensor::YZ);
PpMesoYZ.write();


volTensorField sigmapMeso
(
    IOobject
    (
        "sigmapMeso", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
);

sigmapMeso =  tensor(I)*PpMeso - diffTildeeAlppUsUs; //tensor(I) * Pp - ( tildeSigmaKin - sigmaColl );
// Write into the results folder
sigmapMeso.write();

volScalarField mupMeso
(
    IOobject
    (
        "mupMeso", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

volScalarField mupMesoIncomp
(
    IOobject
    (
        "mupMesoIncomp", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

volTensorField SMeso("SMeso",fvc::grad(tildeUs) + fvc::grad(tildeUs)().T()-2./3.*tr(fvc::grad(tildeUs))*tensor(I));
volTensorField SMesoIncomp("SMesoIncomp",fvc::grad(tildeUs) + fvc::grad(tildeUs)().T());

dimensionedScalar SSsmallMeso("zero", dimensionSet(0,0,-2,0,0,0,0), SMALL);
dimensionedScalar SsmallMeso("zero", dimensionSet(0,0,-1,0,0,0,0), SMALL);

//mup = ( sigmap && S ) / ( max ( S && S, SSsmall ) );
mupMeso = ( sigmapMeso && SMeso ) / ( max ( SMeso && SMeso, SSsmallMeso ) );
mupMesoIncomp = ( sigmapMeso && SMesoIncomp ) / ( max ( SMesoIncomp && SMesoIncomp, SSsmallMeso ) );
// Limit by zero
//mup.max(0.);
// Write into the results folder
mupMeso.write(); 
mupMesoIncomp.write(); 

// INCORRECT MesoTemperature (keep here but it is incorrect)
volScalarField mesoTemp
(
    IOobject
    (
        "mesoTemp",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1./2.*tr((mappedUs-tildeUs)*(mappedUs-tildeUs))
);
// Write into the results folder
mesoTemp.write();

// CORRECT MesoTemperature  (use this one)
//tildeUsPrime = convKernel.coarseningVecFavre("tildeUsPrime",filterWidth[fWidth],1.-mappedVoidfraction,mappedUs-tildeUs);
volScalarField mesoTempCORR
(
    IOobject
    (
        "mesoTempCORR",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,2,-2,0,0), scalar(0) )
);
// Write into the results folder
mesoTempCORR = convKernel.coarseningScaFavre("mesoTempCORR",filterWidth[fWidth],1.-mappedVoidfraction, 1./2.*(mappedUs-tildeUs) & (mappedUs-tildeUs) );
mesoTempCORR.write();


volScalarField mesoTempOzel
(
    IOobject
    (
        "mesoTempOzel",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);
// Write into the results folder
mesoTempOzel = PpMeso/rhoParticle/max(barAlpp, SMALL)/2.0;
mesoTempOzel.write();

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

volScalarField mesoShearRateIncomp
(
    IOobject
    (
        "mesoShearRateIncomp", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRateIncomp =  sqrt( 2.0 * (SMesoIncomp/2.0) && (SMesoIncomp/2.0));
mesoShearRateIncomp.write();


volScalarField mesoShearRateZeroed
(
    IOobject
    (
        "mesoShearRateZeroed", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRateZeroed =  sqrt( SMeso.component(symmTensor::XY)*SMeso.component(symmTensor::XY)+SMeso.component(symmTensor::XZ)*SMeso.component(symmTensor::XZ)+SMeso.component(symmTensor::YZ)*SMeso.component(symmTensor::YZ));
mesoShearRateZeroed.write();

volScalarField mesoShearRateXY
(
    IOobject
    (
        "mesoShearRateXY", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRateXY =  sqrt( SMeso.component(symmTensor::XY)*SMeso.component(symmTensor::XY));
mesoShearRateXY.write();

volScalarField mesoShearRateXZ
(
    IOobject
    (
        "mesoShearRateXZ", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRateXZ =  sqrt( SMeso.component(symmTensor::XZ)*SMeso.component(symmTensor::XZ));
mesoShearRateXZ.write();

volScalarField mesoShearRateYZ
(
    IOobject
    (
        "mesoShearRateYZ", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

mesoShearRateYZ =  sqrt( SMeso.component(symmTensor::YZ)*SMeso.component(symmTensor::YZ));
mesoShearRateYZ.write();


volScalarField mupMesoXY
(
    IOobject
    (
        "mupMesoXY", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mupMesoXY =  sqrt( PpMesoXY*PpMesoXY)/max(mesoShearRateXY,SsmallMeso);
mupMesoXY.write();

volScalarField mupMesoXZ
(
    IOobject
    (
        "mupMesoXZ", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mupMesoXZ =  sqrt( PpMesoXZ*PpMesoXZ)/max(mesoShearRateXZ,SsmallMeso);
mupMesoXZ.write();

volScalarField mupMesoYZ
(
    IOobject
    (
        "mupMesoYZ", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mupMesoYZ =  sqrt( PpMesoYZ*PpMesoYZ)/max(mesoShearRateYZ,SsmallMeso);
mupMesoYZ.write();





volScalarField mesoShearStress
(
    IOobject
    (
        "mesoShearStress", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

mesoShearStress =  sqrt( 1.0/2.0 * (sigmapMeso) && (sigmapMeso));
mesoShearStress.write();

volScalarField mesoShearStressZeroed
(
    IOobject
    (
        "mesoShearStressZeroed", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

mesoShearStressZeroed =  sqrt( 1.0/2.0 * (2.0*PpMesoXY*PpMesoXY+2.0*PpMesoXZ*PpMesoXZ+2.0*PpMesoYZ*PpMesoYZ));
mesoShearStressZeroed.write();

volScalarField barsigmagradus
(
    IOobject
    (
        "barsigmagradus", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-3,0,0), scalar(0) )
);

barsigmagradus = -diffTildeeAlppUsUs&&fvc::grad(tildeUs);

/*
volScalarField tildeUsUsXX = convKernel.coarseningSca("tildeUsUsXX",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(0));
volScalarField tildeUsUsYY = convKernel.coarseningSca("tildeUsUsYY",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(1)*mappedUs.component(1));
volScalarField tildeUsUsZZ = convKernel.coarseningSca("tildeUsUsZZ",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(2)*mappedUs.component(2));
volScalarField tildeUsUsXY = convKernel.coarseningSca("tildeUsUsXY",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(1));
volScalarField tildeUsUsXZ = convKernel.coarseningSca("tildeUsUsXZ",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(0)*mappedUs.component(2));
volScalarField tildeUsUsYZ = convKernel.coarseningSca("tildeUsUsYZ",filterWidth[fWidth],(1-mappedVoidfraction)*mappedUs.component(1)*mappedUs.component(2));
*/
volScalarField tildeUsUsXX = convKernel.coarseningScaFavre("tildeUsUsXX",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(0));
volScalarField tildeUsUsYY = convKernel.coarseningScaFavre("tildeUsUsYY",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(1)*mappedUs.component(1));
volScalarField tildeUsUsZZ = convKernel.coarseningScaFavre("tildeUsUsZZ",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(2)*mappedUs.component(2));
volScalarField tildeUsUsXY = convKernel.coarseningScaFavre("tildeUsUsXY",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(1));
volScalarField tildeUsUsXZ = convKernel.coarseningScaFavre("tildeUsUsXZ",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(2));
volScalarField tildeUsUsYZ = convKernel.coarseningScaFavre("tildeUsUsYZ",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(1)*mappedUs.component(2));

forAll(mesh.cells(),cellI)
{
	tildeUsUs[cellI] =   	symmTensor(		
			  		  tildeUsUsXX[cellI],
			  		  tildeUsUsXY[cellI],
			  		  tildeUsUsXZ[cellI],
			  		  tildeUsUsYY[cellI],
			  		  tildeUsUsYZ[cellI],
			  		  tildeUsUsZZ[cellI] 
					);
}

volScalarField tildetUstUsXX = convKernel.coarseningScaFavre("tildetUstUsXX",filterWidth[fWidth]*3,barAlpp,tildeUs.component(0)*tildeUs.component(0));
volScalarField tildetUstUsYY = convKernel.coarseningScaFavre("tildetUstUsYY",filterWidth[fWidth]*3,barAlpp,tildeUs.component(1)*tildeUs.component(1));
volScalarField tildetUstUsZZ = convKernel.coarseningScaFavre("tildetUstUsZZ",filterWidth[fWidth]*3,barAlpp,tildeUs.component(2)*tildeUs.component(2));
volScalarField tildetUstUsXY = convKernel.coarseningScaFavre("tildetUstUsXY",filterWidth[fWidth]*3,barAlpp,tildeUs.component(0)*tildeUs.component(1));
volScalarField tildetUstUsXZ = convKernel.coarseningScaFavre("tildetUstUsXZ",filterWidth[fWidth]*3,barAlpp,tildeUs.component(0)*tildeUs.component(2));
volScalarField tildetUstUsYZ = convKernel.coarseningScaFavre("tildetUstUsYZ",filterWidth[fWidth]*3,barAlpp,tildeUs.component(1)*tildeUs.component(2));

forAll(mesh.cells(),cellI)
{
	tildetUstUs[cellI] =   	symmTensor(		
			  		  tildetUstUsXX[cellI],
			  		  tildetUstUsXY[cellI],
			  		  tildetUstUsXZ[cellI],
			  		  tildetUstUsYY[cellI],
			  		  tildetUstUsYZ[cellI],
			  		  tildetUstUsZZ[cellI] 
					);
}


volSymmTensorField filterSGSTensor
(
    IOobject
    (
        "filterSGSTensor"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
);

filterSGSTensor = tildeUsUs - symm(tildeUs*tildeUs);
filterSGSTensor.write();

volSymmTensorField testSGSTensor
(
    IOobject
    (
        "testSGSTensor"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
);

testSGSTensor =  tildetUstUs - symm(tildetildeUs*tildetildeUs);
testSGSTensor.write();


volScalarField filterSGSTensorZZ
(
    IOobject
    (
        "filterSGSTensorZZ"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,2,-2,0,0), scalar(0) )
);

filterSGSTensorZZ =  filterSGSTensor.component(symmTensor::ZZ);
filterSGSTensorZZ.write();

volScalarField testSGSTensorZZ
(
    IOobject
    (
        "testSGSTensorZZ"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,2,-2,0,0), scalar(0) )
);

testSGSTensorZZ =  testSGSTensor.component(symmTensor::ZZ);
testSGSTensorZZ.write();


volScalarField tildetUstUsHatStarXX = convKernel.coarseningSca("tildetUstUsHatStarXX",filterWidth[fWidth]*3,tildeUs.component(0)*tildeUs.component(0));
volScalarField tildetUstUsHatStarYY = convKernel.coarseningSca("tildetUstUsHatStarYY",filterWidth[fWidth]*3,tildeUs.component(1)*tildeUs.component(1));
volScalarField tildetUstUsHatStarZZ = convKernel.coarseningSca("tildetUstUsHatStarZZ",filterWidth[fWidth]*3,tildeUs.component(2)*tildeUs.component(2));
volScalarField tildetUstUsHatStarXY = convKernel.coarseningSca("tildetUstUsHatStarXY",filterWidth[fWidth]*3,tildeUs.component(0)*tildeUs.component(1));
volScalarField tildetUstUsHatStarXZ = convKernel.coarseningSca("tildetUstUsHatStarXZ",filterWidth[fWidth]*3,tildeUs.component(0)*tildeUs.component(2));
volScalarField tildetUstUsHatStarYZ = convKernel.coarseningSca("tildetUstUsHatStarYZ",filterWidth[fWidth]*3,tildeUs.component(1)*tildeUs.component(2));

volVectorField tildetildeUsHatStar = convKernel.coarseningVec("tildetildeUsHatStar",filterWidth[fWidth]*3,tildeUs);

forAll(mesh.cells(),cellI)
{
	tildetUstUsHatStar[cellI] =   	symmTensor(		
			  		  tildetUstUsHatStarXX[cellI],
			  		  tildetUstUsHatStarXY[cellI],
			  		  tildetUstUsHatStarXZ[cellI],
			  		  tildetUstUsHatStarYY[cellI],
			  		  tildetUstUsHatStarYZ[cellI],
			  		  tildetUstUsHatStarZZ[cellI] 
					);
}

volSymmTensorField testSGSTensorHatStar
(
    IOobject
    (
        "testSGSTensorHatStar"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
);

testSGSTensorHatStar =  tildetUstUsHatStar - symm(tildetildeUsHatStar*tildetildeUsHatStar);
testSGSTensorHatStar.write();


// Avoid dividing by alppp
volScalarField barBarAlpptUsxtUsx = convKernel.coarseningSca("barBarAlpptUsxtUsx",filterWidth[fWidth]*3,barAlpp*tildeUs.component(0)*tildeUs.component(0));
volScalarField barBarAlpptUsytUsy = convKernel.coarseningSca("barBarAlpptUsytUsy",filterWidth[fWidth]*3,barAlpp*tildeUs.component(1)*tildeUs.component(1));
volScalarField barBarAlpptUsztUsz = convKernel.coarseningSca("barBarAlpptUsztUsz",filterWidth[fWidth]*3,barAlpp*tildeUs.component(2)*tildeUs.component(2));

volScalarField barBarAlpptUsxtUsy = convKernel.coarseningSca("barBarAlpptUsxtUsy",filterWidth[fWidth]*3,barAlpp*tildeUs.component(0)*tildeUs.component(1));
volScalarField barBarAlpptUsxtUsz = convKernel.coarseningSca("barBarAlpptUsxtUsz",filterWidth[fWidth]*3,barAlpp*tildeUs.component(0)*tildeUs.component(2));
volScalarField barBarAlpptUsytUsz = convKernel.coarseningSca("barBarAlpptUsytUsz",filterWidth[fWidth]*3,barAlpp*tildeUs.component(1)*tildeUs.component(2));

volSymmTensorField diffbarBarAlpptUstUs
(
    IOobject
    (
        "diffbarBarAlpptUstUs",
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
	diffbarBarAlpptUstUs[cellI] =   rhoParticle
				    * symmTensor
					(		
			  		  barBarAlpptUsxtUsx[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][0]*tildetildeUs[cellI][0] ),
			  		  barBarAlpptUsxtUsy[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][0]*tildetildeUs[cellI][1] ),
			  		  barBarAlpptUsxtUsz[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][0]*tildetildeUs[cellI][2] ),
			  		  barBarAlpptUsytUsy[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][1]*tildetildeUs[cellI][1] ),
			  		  barBarAlpptUsytUsz[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][1]*tildetildeUs[cellI][2] ),
			  		  barBarAlpptUsztUsz[cellI] - barBarAlpp[cellI]*( tildetildeUs[cellI][2]*tildetildeUs[cellI][2] )
					);
}
diffbarBarAlpptUstUs.write();
// Define correlation coefficient kernel
//correlationCoeff corr(mesh,false);
//Info << corr.calc(correlationXX.component(symmTensor::XZ),correlationXX(component::XX)) << endl;
