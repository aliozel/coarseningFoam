  // Calculate Eulerian resolved drag force
Info << "\n" << endl;
Info << tab << "Calculating Eulerian resolved drag force..." << endl;
resolvedDragFEulerian().setForceEulerian(barAlpf,tildeUf,radii[me],tildeUs);

// Calculate normalizedEulerianDragCoeff
volVectorField normalizedEulerianDragCoeff
(
    IOobject
    (
        "normalizedEulerianDragCoeff"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedEulerianDragCoeffNum
(
    IOobject
    (
        "normalizedEulerianDragCoeffNum"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedEulerianDragCoeffDom
(
    IOobject
    (
        "normalizedEulerianDragCoeffDom"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


forAll(mesh.cells(),index)
 {	
      for(label iDir=0; iDir<3; iDir++)
      {		
	   {
	    if( resolvedDragFEulerian().FpEulerian(index)[iDir] != 0 ) 
	       normalizedEulerianDragCoeff[index][iDir] =  barForceEulerian[index][iDir] 
				  			 / resolvedDragFEulerian().FpEulerian(index)[iDir];	
	       normalizedEulerianDragCoeffNum[index][iDir] = barForceEulerian[index][iDir];
	       normalizedEulerianDragCoeffDom[index][iDir] = resolvedDragFEulerian().FpEulerian(index)[iDir];	
	   }
      }
 }
 normalizedEulerianDragCoeff.write();
 normalizedEulerianDragCoeffNum.write();	
 normalizedEulerianDragCoeffDom.write();	


volVectorField normalizedAlppAlpfEulerianDragCoeff
(
    IOobject
    (
        "normalizedAlppAlpfEulerianDragCoeff"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedAlppAlpfEulerianDragCoeffNum
(
    IOobject
    (
        "normalizedAlppAlpfEulerianDragCoeffNum"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedAlppAlpfEulerianDragCoeffDom
(
    IOobject
    (
        "normalizedAlppAlpfEulerianDragCoeffDom"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


forAll(mesh.cells(),index)
{	
      for(label iDir=0; iDir<3; iDir++)
      {		
	    // Ratio between filtered drag coefficient and resolved drag coeff.
	   //if( resolvedDragFEulerian().FpEulerian(index)[iDir] > SMALL ) 
	   {			      

	    if( resolvedDragFEulerian().FpEulerian(index)[iDir] != 0 ) 

	       normalizedAlppAlpfEulerianDragCoeff[index][iDir] =  barAlppAlpfForceEulerian[index][iDir] 
		 / resolvedDragFEulerian().FpEulerian(index)[iDir]/(barAlpp[index]/barAlpf[index]);	

	       normalizedAlppAlpfEulerianDragCoeffNum[index][iDir] = barAlppAlpfForceEulerian[index][iDir];
	       normalizedAlppAlpfEulerianDragCoeffDom[index][iDir] = barAlpp[index]/barAlpf[index]*resolvedDragFEulerian().FpEulerian(index)[iDir];	


	   }
      }
 }

normalizedAlppAlpfEulerianDragCoeff.write();
normalizedAlppAlpfEulerianDragCoeffNum.write();	
normalizedAlppAlpfEulerianDragCoeffDom.write();	



// Particle density
scalar rhoParticle = resolvedDragFEulerian().rhoParticle(); 
dimensionedScalar rhoParticleDim("0",dimensionSet(1,-3,0,0,0),rhoParticle); 
Info << tab << "rhoParticle = " << resolvedDragFEulerian().rhoParticle() << endl; 
//Info << tab << "rho = " << rho << endl;
// Gas kinematic viscosity
const volScalarField nuf = resolvedDragFEulerian().nuf(); 
//Info << tab << "nuf = [" << min(nuf).value() <<":"<<max(nuf).value()<<"]" <<endl;
 
//barAlppUfUf
volScalarField barAlppUfxUfx = convKernel.coarseningSca("barAlppUfxUfx",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(0)*U.component(0));
volScalarField barAlppUfyUfy = convKernel.coarseningSca("barAlppUfyUfy",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(1)*U.component(1));
volScalarField barAlppUfzUfz = convKernel.coarseningSca("barALppUfzUfz",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(2)*U.component(2));

volScalarField barAlppUfxUfy = convKernel.coarseningSca("barAlppUfxUfy",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(0)*U.component(1));
volScalarField barAlppUfxUfz = convKernel.coarseningSca("barAlppUfxUfz",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(0)*U.component(2));
volScalarField barAlppUfyUfz = convKernel.coarseningSca("barAlppUfyUfz",filterWidth[fWidth],(1-mappedVoidfraction)*U.component(1)*U.component(2));

//Asgs term - diffbarAlppUfUf
volSymmTensorField diffbarAlppUfUf
(
    IOobject
    (
     "diffbarAlppUfUf",
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
        diffbarAlppUfUf[cellI] = symmTensor
	  (
	   barAlppUfxUfx[cellI] - barAlppUf[cellI][0]*tildeUf[cellI][0],
	   barAlppUfxUfy[cellI] - barAlppUf[cellI][0]*tildeUf[cellI][1],
	   barAlppUfxUfz[cellI] - barAlppUf[cellI][0]*tildeUf[cellI][2],
	   barAlppUfyUfy[cellI] - barAlppUf[cellI][1]*tildeUf[cellI][1],
	   barAlppUfyUfz[cellI] - barAlppUf[cellI][1]*tildeUf[cellI][2],
	   barAlppUfzUfz[cellI] - barAlppUf[cellI][2]*tildeUf[cellI][2]
	  );
}

// Divergence of Asgs term
volVectorField divAsgs = fvc::div(diffbarAlppUfUf);

volVectorField magDivAsgs
(
    IOobject
    (
        "magDivAsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magDivAsgs.internalField() = cmptMag(divAsgs.internalField());

//Info << "\nDivergence of Assg = " << fvc::domainIntegrate(divAsgs).value() << endl;


//resolved Bsgs term
volVectorField tildeUfdivtildeUf = tildeUf*fvc::div(tildeUf);

//Bsgs microscale part - barUfdivUf
volVectorField barUfdivUf = convKernel.coarseningVec("barUfdivUf",filterWidth[fWidth],U*fvc::div(U));

//Bsgs
volVectorField Bsgs = barUfdivUf-tildeUfdivtildeUf;

volVectorField magBsgs
(
    IOobject
    (
        "magBsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magBsgs.internalField() = cmptMag(Bsgs.internalField());


//Csgs resolved part
volVectorField tildeUfDivbarAlpftildeUf = tildeUf*fvc::div(barAlpf*tildeUf);

//Csgs microscale part - barUfDivAlpfUf
volVectorField barUfDivAlpfUf = convKernel.coarseningVec("barUfDivAlpfUf",filterWidth[fWidth],U*fvc::div(mappedVoidfraction*U));

//Csgs
volVectorField Csgs = barUfDivAlpfUf-tildeUfDivbarAlpftildeUf;

volVectorField magCsgs
(
    IOobject
    (
        "magCsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magCsgs.internalField() = cmptMag(Csgs.internalField());




//Dsgs resolved part
volVectorField tildeUfDivbarAlpptildeUs = tildeUf*fvc::div(barAlpp*tildeUs);

//Dsgs microscale part - barUfDivAlppUs
volVectorField barUfDivAlppUs = convKernel.coarseningVec("barUfDivAlppUs",filterWidth[fWidth],U*fvc::div((1-mappedVoidfraction)*mappedUs));

//Dsgs 
volVectorField Dsgs = barUfDivAlppUs - tildeUfDivbarAlpptildeUs;

volVectorField magDsgs
(
    IOobject
    (
        "magDsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magDsgs.internalField() = cmptMag(Dsgs.internalField());



//Esgs resolved part
volVectorField barAlpptildeGradPg = barAlpp*tildeGradPg;

//Esgs microscale part - barAlppGradPg
volVectorField barAlppGradPg = convKernel.coarseningVec("barAlppGradPg",filterWidth[fWidth],(1-mappedVoidfraction)*gradPg);

//Esgs
volVectorField Esgs = barAlppGradPg - barAlpptildeGradPg;								

volVectorField magEsgs
(
    IOobject
    (
        "magEsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magEsgs.internalField() = cmptMag(Esgs.internalField());

						  

//Fsgs resolved part
volTensorField gradtildeUf = fvc::grad(tildeUf);
volTensorField tildeTau = nuf*(gradtildeUf+gradtildeUf.T()-2.0/3.0*fvc::div(tildeUf)*tensor(I));
volVectorField barAlppDivTildeTau = barAlpp*fvc::div(tildeTau);

//Fsgs microscale part
volTensorField gradU = fvc::grad(U);
volTensorField Tau = nuf*(gradU+gradU.T()-2.0/3.0*fvc::div(U)*tensor(I));
volVectorField barAlppDivTau = convKernel.coarseningVec("barAlppDivTau",filterWidth[fWidth],(1-mappedVoidfraction)*fvc::div(Tau));

//Fsgs
volVectorField Fsgs = barAlppDivTau - barAlppDivTildeTau;

volVectorField magFsgs
(
    IOobject
    (
        "magFsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magFsgs.internalField() = cmptMag(Fsgs.internalField());



//Gsgs resolved part
//volVectorField resolvedDragF(U-mappedUs);
volVectorField resolvedDragF
(
    IOobject
    (
        "resolvedDragF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimensionSet(1, -2, -2, 0, 0), vector(0,0,0))
);
forAll(mesh.cells(),cellI)
{
  resolvedDragF[cellI] = -resolvedDragFEulerian().FpEulerian(cellI);
}
resolvedDragF.write();

volVectorField AlppoverAlpftildeF = barAlpp/barAlpf*resolvedDragF/rho;

//Gsgs microscale part
volVectorField barForceEulerian = convKernel.coarseningVec("barForceEulerian",filterWidth[fWidth],-(1-mappedVoidfraction)*dragFpEulerian/mappedVoidfraction);
volVectorField barDragF
(
    IOobject
    (
        "barDragF"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barForceEulerian/rho
); 
barDragF.write();
//volVectorField barDragF = barForceEulerian/rho;

//Gsgs
volVectorField Gsgs = barDragF - AlppoverAlpftildeF;

volVectorField magGsgs
(
    IOobject
    (
        "magGsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);


magGsgs.internalField() = cmptMag(Gsgs.internalField());



//Hsgs - diffbarAlppUfUf2
volSymmTensorField diffbarAlppUfUf2
(
    IOobject
    (
     "diffbarAlppUfUf2",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
 );

forAll(mesh.cells(),cellI)
{
        diffbarAlppUfUf2[cellI] = symmTensor
	  (
	   barAlppUfxUfx[cellI] - barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][0],
	   barAlppUfxUfy[cellI] - barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][1],
	   barAlppUfxUfz[cellI] - barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][2],
	   barAlppUfyUfy[cellI] - barAlpp[cellI]*tildeUf[cellI][1]*tildeUf[cellI][1],
	   barAlppUfyUfz[cellI] - barAlpp[cellI]*tildeUf[cellI][1]*tildeUf[cellI][2],
	   barAlppUfzUfz[cellI] - barAlpp[cellI]*tildeUf[cellI][2]*tildeUf[cellI][2]
	  );
}

// Divergence of Hsgs term
volVectorField divHsgs = fvc::div(diffbarAlppUfUf2);
//Info << "\nDivergence of Assg = " << fvc::domainIntegrate(divAsgs).value() << endl;


volVectorField magDivHsgs
(
    IOobject
    (
        "magDivHsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magDivHsgs.internalField() = cmptMag(divHsgs.internalField());


//UfUf
volScalarField barUfxUfx = convKernel.coarseningSca("barUfxUfx",filterWidth[fWidth],U.component(0)*U.component(0));
volScalarField barUfyUfy = convKernel.coarseningSca("barUfyUfy",filterWidth[fWidth],U.component(1)*U.component(1));
volScalarField barUfzUfz = convKernel.coarseningSca("barUfzUfz",filterWidth[fWidth],U.component(2)*U.component(2));

volScalarField barUfxUfy = convKernel.coarseningSca("barUfxUfy",filterWidth[fWidth],U.component(0)*U.component(1));
volScalarField barUfxUfz = convKernel.coarseningSca("barUfxUfz",filterWidth[fWidth],U.component(0)*U.component(2));
volScalarField barUfyUfz = convKernel.coarseningSca("barUfyUfz",filterWidth[fWidth],U.component(1)*U.component(2));

//Isgs term - diffUfUf
volSymmTensorField diffUfUf
(
    IOobject
    (
     "diffUfUf",
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
        diffUfUf[cellI] = symmTensor
	  (
	   barUfxUfx[cellI] - tildeUf[cellI][0]*tildeUf[cellI][0],
	   barUfxUfy[cellI] - tildeUf[cellI][0]*tildeUf[cellI][1],
	   barUfxUfz[cellI] - tildeUf[cellI][0]*tildeUf[cellI][2],
	   barUfyUfy[cellI] - tildeUf[cellI][1]*tildeUf[cellI][1],
	   barUfyUfz[cellI] - tildeUf[cellI][1]*tildeUf[cellI][2],
	   barUfzUfz[cellI] - tildeUf[cellI][2]*tildeUf[cellI][2]
	  );
}

// Divergence of Isgs term
volVectorField divIsgs = fvc::div(diffUfUf);


volVectorField magDivIsgs
(
    IOobject
    (
        "magDivIsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magDivIsgs.internalField() = cmptMag(divIsgs.internalField());


//Jsgs resolved part
volVectorField barAlpftildeGradPg = barAlpf*tildeGradPg;

//Jsgs microscale part - barAlpfGradPg
volVectorField barAlpfGradPg = convKernel.coarseningVec("barAlpfGradPg",filterWidth[fWidth],mappedVoidfraction*gradPg);

//Jsgs
volVectorField Jsgs = barAlpfGradPg - barAlpftildeGradPg;

//Jsgs term
volVectorField AlppJsgs = barAlpp/barAlpf*Jsgs;


volVectorField magAlppJsgs
(
    IOobject
    (
        "magAlppJsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magAlppJsgs.internalField() = cmptMag(AlppJsgs.internalField());


//Ksgs resolved part
volVectorField barAlpfDivTildeTau = barAlpf*fvc::div(tildeTau);

//Ksgs microscale part
volVectorField barAlpfDivTau = convKernel.coarseningVec("barAlpfDivTau",filterWidth[fWidth],mappedVoidfraction*fvc::div(Tau));

//Ksgs
volVectorField Ksgs = barAlpfDivTau - barAlpfDivTildeTau;

//Ksgs term
volVectorField AlppKsgs = barAlpp/barAlpf*Ksgs;

volVectorField magAlppKsgs
(
    IOobject
    (
        "magAlppKsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magAlppKsgs.internalField() = cmptMag(AlppKsgs.internalField());


														  

//barAlpfUfUf
volScalarField barAlpfUfxUfx = convKernel.coarseningSca("barAlpfUfxUfx",filterWidth[fWidth],mappedVoidfraction*U.component(0)*U.component(0));
volScalarField barAlpfUfyUfy = convKernel.coarseningSca("barAlpfUfyUfy",filterWidth[fWidth],mappedVoidfraction*U.component(1)*U.component(1));
volScalarField barAlpfUfzUfz = convKernel.coarseningSca("barALpfUfzUfz",filterWidth[fWidth],mappedVoidfraction*U.component(2)*U.component(2));

volScalarField barAlpfUfxUfy = convKernel.coarseningSca("barAlpfUfxUfy",filterWidth[fWidth],mappedVoidfraction*U.component(0)*U.component(1));
volScalarField barAlpfUfxUfz = convKernel.coarseningSca("barAlpfUfxUfz",filterWidth[fWidth],mappedVoidfraction*U.component(0)*U.component(2));
volScalarField barAlpfUfyUfz = convKernel.coarseningSca("barAlpfUfyUfz",filterWidth[fWidth],mappedVoidfraction*U.component(1)*U.component(2));

//Lsgs term - diffbarAlpfUfUf
volSymmTensorField diffbarAlpfUfUf
(
    IOobject
    (
     "diffbarAlpfUfUf",
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
        diffbarAlpfUfUf[cellI] = symmTensor
	  (
	   barAlpfUfxUfx[cellI] - barAlpf[cellI]*tildeUf[cellI][0]*tildeUf[cellI][0],
	   barAlpfUfxUfy[cellI] - barAlpf[cellI]*tildeUf[cellI][0]*tildeUf[cellI][1],
	   barAlpfUfxUfz[cellI] - barAlpf[cellI]*tildeUf[cellI][0]*tildeUf[cellI][2],
	   barAlpfUfyUfy[cellI] - barAlpf[cellI]*tildeUf[cellI][1]*tildeUf[cellI][1],
	   barAlpfUfyUfz[cellI] - barAlpf[cellI]*tildeUf[cellI][1]*tildeUf[cellI][2],
	   barAlpfUfzUfz[cellI] - barAlpf[cellI]*tildeUf[cellI][2]*tildeUf[cellI][2]
	  );
}

// Divergence of Lsgs term
volVectorField divLsgs = fvc::div(diffbarAlpfUfUf)/barAlpf;


volVectorField magDivLsgs
(
    IOobject
    (
        "magDivLsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magDivLsgs.internalField() = cmptMag(divLsgs.internalField());


//Msgs resolved part
//resolvedDragF

//Msgs microscale part
volVectorField barF = convKernel.coarseningVec("barF",filterWidth[fWidth],-dragFpEulerian);
volVectorField barFoverrho = barF/rho;

//Msgs
volVectorField Msgs = barFoverrho - resolvedDragF/rho;

//Msgs term
volVectorField AlppMsgs = barAlpp/barAlpf*Msgs;


volVectorField magAlppMsgs
(
    IOobject
    (
        "magAlppMsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magAlppMsgs.internalField() = cmptMag(AlppMsgs.internalField());


//LHS diffbarAlppUfUfbarAlppUfUf
volSymmTensorField diffbarAlppUfUfbarAlppUfUf
(
    IOobject
    (
     "diffbarAlppUfUfbarAlppUfUf",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
 );

forAll(mesh.cells(),cellI)
{
        diffbarAlppUfUfbarAlppUfUf[cellI] = symmTensor
	  (
	   barAlppUf[cellI][0]*tildeUf[cellI][0]-barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][0],
	   barAlppUf[cellI][0]*tildeUf[cellI][1]-barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][1],
	   barAlppUf[cellI][0]*tildeUf[cellI][2]-barAlpp[cellI]*tildeUf[cellI][0]*tildeUf[cellI][2],
	   barAlppUf[cellI][1]*tildeUf[cellI][1]-barAlpp[cellI]*tildeUf[cellI][1]*tildeUf[cellI][1],
	   barAlppUf[cellI][1]*tildeUf[cellI][2]-barAlpp[cellI]*tildeUf[cellI][1]*tildeUf[cellI][2],
	   barAlppUf[cellI][2]*tildeUf[cellI][2]-barAlpp[cellI]*tildeUf[cellI][2]*tildeUf[cellI][2]
	  );
}


//Xsgs microscale part
//volVectorField barAlpfGradPg = convKernel.coarseningVec("barAlpfGradPg",filterWidth[fWidth],mappedVoidfraction*gradPg/rho);

//Xsgs resolved part
volVectorField barGradPg = convKernel.coarseningVec("barGradPg",filterWidth[fWidth],gradPg);
volVectorField barAlpfbarGradPg = barAlpf*barGradPg;

//Xsgs
//volVectorField Xsgs = barAlpfGradPg - barAlpfbarGradPg;
//write Xsgs
volVectorField Xsgs
(
    IOobject
    (
        "Xsgs"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barAlpfGradPg - barAlpfbarGradPg
); 
Xsgs.write();



//1/AlpfXsgs
//volVectorField AlpfXsgs = Xsgs/barAlpf;
// Write gradBarAlpp
volVectorField AlpfXsgs
(
    IOobject
    (
        "AlpfXsgs"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Xsgs/barAlpf
); 
AlpfXsgs.write();

volVectorField magAlpfXsgs
(
    IOobject
    (
        "magAlpfXsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magAlpfXsgs.internalField() = cmptMag(AlpfXsgs.internalField());





//Ysgs microscale part
//same as Gsgs microscale - barDragF

//Ysgs resolved part
//volVectorField barF = convKernel.coarseningVec("barF",filterWidth[fWidth],-dragFpEulerian);
volVectorField barAlppbarAlpfbarFrho = barAlpp/barAlpf*barF/rho;

//Ysgs
//volVectorField Ysgs = barDragF - barAlppbarAlpfbarFrho;

volVectorField Ysgs
(
    IOobject
    (
        "Ysgs"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barDragF - barAlppbarAlpfbarFrho
); 
Ysgs.write();




volVectorField magYsgs
(
    IOobject
    (
        "magYsgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magYsgs.internalField() = cmptMag(Ysgs.internalField());




//divLHS
volVectorField divLHS = fvc::div(diffbarAlppUfUfbarAlppUfUf);


volVectorField magDivLHS
(
    IOobject
    (
        "magDivLHS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

magDivLHS.internalField() = cmptMag(divLHS.internalField());



//barAlppG
volVectorField barAlppG = barAlpp*g;

//Addedsourceterm
dimensionedScalar aveRho = (fvc::domainIntegrate(barAlpf*rho).value() + fvc::domainIntegrate(barAlpp*rhoParticleDim).value())/gSum(mesh.V());
//Info << " aveRho = " << aveRho <<endl;
volVectorField alppbarrho = -barAlpp/rho*aveRho*g;

//gradient of barAlpp
/*volVectorField gradBarAlpp
(
    IOobject
    (
        "gradBarAlpp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(barAlpp)
);

gradBarAlpp.write();
*/
/*
volVectorField Vdrift
(
    IOobject
    (
        "Vdrift"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barAlppVdrift/barAlpp
);
Vdrift.write(); 

volVectorField magVdrift
(
    IOobject
    (
        "magVdrift",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), vector(0,0,0) )
);

magVdrift.internalField() = cmptMag(Vdrift.internalField());
*/											  

/* - AO
// Write results into a file
*outputFile	<< runTime.timeName() << tab 
		<< fvc::domainIntegrate(divAsgs).value()[2]  					<< tab  
                << fvc::domainIntegrate(Bsgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(Csgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(Dsgs).value()[2]                		 	<< tab                
		<< fvc::domainIntegrate(Esgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(Fsgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(Gsgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(AlppJsgs).value()[2]                 			<< tab
                << fvc::domainIntegrate(AlppKsgs).value()[2]                 			<< tab
		<< fvc::domainIntegrate(AlppMsgs).value()[2]					<< tab
                << fvc::domainIntegrate(divLsgs).value()[2]					<< tab
                << fvc::domainIntegrate(divIsgs).value()[2]					<< tab
                << fvc::domainIntegrate(divHsgs).value()[2]					<< tab
                << fvc::domainIntegrate(divLHS).value()[2]					<< tab
                << fvc::domainIntegrate(AlpfXsgs).value()[2]					<< tab
                << fvc::domainIntegrate(Ysgs).value()[2]					<< tab
<< endl;
AO - */ 

/*
domainAveVars[timeI][fWidth][0] = fvc::domainIntegrate(divAsgs).value()[2];

domainAveVars[timeI][fWidth][1] = fvc::domainIntegrate(Bsgs).value()[2];          
domainAveVars[timeI][fWidth][2] = fvc::domainIntegrate(Csgs).value()[2];          
domainAveVars[timeI][fWidth][3] = fvc::domainIntegrate(Dsgs).value()[2];                		 	        
domainAveVars[timeI][fWidth][4] = fvc::domainIntegrate(Esgs).value()[2];          
domainAveVars[timeI][fWidth][5] = fvc::domainIntegrate(Fsgs).value()[2];          
domainAveVars[timeI][fWidth][6] = fvc::domainIntegrate(Gsgs).value()[2];          
domainAveVars[timeI][fWidth][7] = fvc::domainIntegrate(AlppJsgs).value()[2];          
domainAveVars[timeI][fWidth][8] = fvc::domainIntegrate(AlppKsgs).value()[2];          
domainAveVars[timeI][fWidth][9] = fvc::domainIntegrate(AlppMsgs).value()[2];					
domainAveVars[timeI][fWidth][10] = fvc::domainIntegrate(divLsgs).value()[2];					
domainAveVars[timeI][fWidth][11] = fvc::domainIntegrate(divIsgs).value()[2];					
domainAveVars[timeI][fWidth][12] = fvc::domainIntegrate(divHsgs).value()[2];					
domainAveVars[timeI][fWidth][13] = fvc::domainIntegrate(divLHS).value()[2];					
domainAveVars[timeI][fWidth][14] = fvc::domainIntegrate(AlpfXsgs).value()[2];					
domainAveVars[timeI][fWidth][15] = fvc::domainIntegrate(Ysgs).value()[2];					
*/


domainAveVars[timeI][fWidth][0] = fvc::domainIntegrate(magDivAsgs).value()[2];
domainAveVars[timeI][fWidth][1] = fvc::domainIntegrate(magBsgs).value()[2];          
domainAveVars[timeI][fWidth][2] = fvc::domainIntegrate(magCsgs).value()[2];          
domainAveVars[timeI][fWidth][3] = fvc::domainIntegrate(magDsgs).value()[2];                		 	         
domainAveVars[timeI][fWidth][4] = fvc::domainIntegrate(magEsgs).value()[2];          
domainAveVars[timeI][fWidth][5] = fvc::domainIntegrate(magFsgs).value()[2];          
domainAveVars[timeI][fWidth][6] = fvc::domainIntegrate(magGsgs).value()[2];          
domainAveVars[timeI][fWidth][7] = fvc::domainIntegrate(magAlppJsgs).value()[2];          
domainAveVars[timeI][fWidth][8] = fvc::domainIntegrate(magAlppKsgs).value()[2];          
domainAveVars[timeI][fWidth][9] = fvc::domainIntegrate(magAlppMsgs).value()[2];					
domainAveVars[timeI][fWidth][10] = fvc::domainIntegrate(magDivLsgs).value()[2];					
domainAveVars[timeI][fWidth][11] = fvc::domainIntegrate(magDivIsgs).value()[2];					
domainAveVars[timeI][fWidth][12] = fvc::domainIntegrate(magDivHsgs).value()[2];					
domainAveVars[timeI][fWidth][13] = fvc::domainIntegrate(magDivLHS).value()[2];					
domainAveVars[timeI][fWidth][14] = fvc::domainIntegrate(magAlpfXsgs).value()[2];					
domainAveVars[timeI][fWidth][15] = fvc::domainIntegrate(magYsgs).value()[2];					

