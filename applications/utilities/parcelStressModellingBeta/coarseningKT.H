// Filtered total Particle Stress
volScalarField tildePartxx = convKernel.coarseningSca("tildePartxx",filterWidth[fWidth],totalParticleStress.component(symmTensor::XX));
volScalarField tildePartyy = convKernel.coarseningSca("tildePartyy",filterWidth[fWidth],totalParticleStress.component(symmTensor::YY));
volScalarField tildePartzz = convKernel.coarseningSca("tildePartzz",filterWidth[fWidth],totalParticleStress.component(symmTensor::ZZ));

volScalarField tildePartxy = convKernel.coarseningSca("tildePartxy",filterWidth[fWidth],totalParticleStress.component(symmTensor::XY));
volScalarField tildePartxz = convKernel.coarseningSca("tildePartxz",filterWidth[fWidth],totalParticleStress.component(symmTensor::XZ));
volScalarField tildePartyz = convKernel.coarseningSca("tildePartyz",filterWidth[fWidth],totalParticleStress.component(symmTensor::YZ));

// Total stress
volSymmTensorField filteredTotalParticleStress
(
    IOobject
    (
	"filteredTotalParticleStress",
	mesh.time().timeName(),
	mesh,
	IOobject::NO_READ, //READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
);

forAll(mesh.cells(),cellI)
{
	filteredTotalParticleStress[cellI] = symmTensor
				     	     (		
						tildePartxx[cellI],
			  	       		tildePartxy[cellI],
			  	       		tildePartxz[cellI],
			  	       		tildePartyy[cellI],
			  	       		tildePartyz[cellI],						 
			  	       		tildePartzz[cellI]
					     );
}

//  Filtered KT pressure
volScalarField filteredPpKT
(
    IOobject
    (
        "filteredPpKT",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);
filteredPpKT = 1./3. * ( tildePartxx + tildePartyy + tildePartzz );
// Write into the results folder
filteredPpKT.write();

//- Filtered KT shear stress
volTensorField filteredSigmapKT
(
    IOobject
    (
        "filteredSigmapKT", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
);
//filteredSigmapKT =  filteredTotalParticleStress - tensor(I)*filteredPpKT; 
// See stressModelling.H Line 145
filteredSigmapKT =  tensor(I)*filteredPpKT - filteredTotalParticleStress; 
// Write into the results folder
filteredSigmapKT.write();

//- Filtered KT viscosity
volScalarField filteredMupKT
(
    IOobject
    (
        "filteredMupKT", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

volTensorField filteredSKT("filteredSKT",fvc::grad(tildeUs) + fvc::grad(tildeUs)().T()-2./3.*tr(fvc::grad(tildeUs))*tensor(I));

dimensionedScalar SSsmallSKT("zero", dimensionSet(0,0,-2,0,0,0,0), SMALL);


filteredMupKT = ( filteredSigmapKT && filteredSKT ) / ( max ( filteredSKT && filteredSKT, SSsmallSKT ) );
// Write into the results folder
filteredMupKT.write(); 



//- Filtered granularTemperature
volScalarField filteredGranularTemp
(
    IOobject
    (
        "filteredGranularTemp", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    granularTemp
);
//filteredGranularTemp = convKernel.coarseningSca("filteredGranularTemp",filterWidth[fWidth],granularTemp);
filteredGranularTemp = convKernel.coarseningScaFavre("filteredGranularTemp",filterWidth[fWidth],1.-mappedVoidfraction,granularTemp);
filteredGranularTemp.write();

//- Filtered Sigma Gradus
volScalarField filteredSigmagradus
(
    IOobject
    (
        "filteredSigmagradus", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sigmagradus
);
//filteredSigmagradus = convKernel.coarseningSca("filteredSigmagradus",filterWidth[fWidth],sigmagradus);
filteredSigmagradus = -filteredTotalParticleStress && fvc::grad(tildeUs);
filteredSigmagradus.write();


