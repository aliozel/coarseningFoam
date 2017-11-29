 volVectorField predictedDriftVelocity
 (
     IOobject
     (
         "predictedADMDriftVelocity"+filterWidthName,
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
 );

 volVectorField computedDriftVelocity
 (
     IOobject
     (
         "computedDriftVelocity"+filterWidthName,
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
 );	     
 // Calculate computed scalar variance
 //volVectorField barAlpp = convKernel.coarseningSca("barAlpp",filterWidthADM,alpp));
 volVectorField tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidthADM,alpp,Us);
 volVectorField tildeUf = convKernel.coarseningVecFavre("tildeUf",filterWidthADM,(1.-alpp),U);
 volVectorField tildeAlppUr = convKernel.coarseningVec("tildeAlppUr",filterWidthADM,alpp*(U-Us));
 computedDriftVelocity = tildeAlppUr - barAlpp*(tildeUf-tildeUs);
 computedDriftVelocity.write();

 // nth order ADM
 if(ADMOrder==0)
 {
 	volVectorField UsStar = tildeUs;
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);		
 }
 else if(ADMOrder==1)
 {
 	volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
 	volScalarField phiStar =     barAlpp 
	                         + ( barAlpp - barbarAlpp );
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
 }
 else if(ADMOrder==2)
 {
 	volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
	volScalarField barbarbarAlpp = convKernel.coarseningSca("barbarbarAlpp",filterWidthADM,barbarAlpp);	
 	volScalarField phiStar =     barAlpp 
	                         + ( barAlpp -    barbarAlpp ) 
				 + ( barAlpp - 2.*barbarAlpp + barbarbarAlpp );
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
 }
 else if(ADMOrder==3)
 {
 	volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
	volScalarField barbarbarAlpp = convKernel.coarseningSca("barbarbarAlpp",filterWidthADM,barbarAlpp);	
	volScalarField barbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarAlpp",filterWidthADM,barbarbarAlpp);	
 	volScalarField phiStar =     barAlpp 
	                         + ( barAlpp -    barbarAlpp ) 
				 + ( barAlpp - 2.*barbarAlpp +    barbarbarAlpp ) 
				 + ( barAlpp - 3.*barbarAlpp + 3.*barbarbarAlpp - barbarbarbarAlpp ) ;
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
 }
 else if(ADMOrder==4)
 {
 	volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
	volScalarField barbarbarAlpp = convKernel.coarseningSca("barbarbarAlpp",filterWidthADM,barbarAlpp);	
	volScalarField barbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarAlpp",filterWidthADM,barbarbarAlpp);	
	volScalarField barbarbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarbarAlpp",filterWidthADM,barbarbarbarAlpp);	
 	volScalarField phiStar =     barAlpp 
	                         + ( barAlpp -    barbarAlpp ) 
				 + ( barAlpp - 2.*barbarAlpp +    barbarbarAlpp ) 
				 + ( barAlpp - 3.*barbarAlpp + 3.*barbarbarAlpp -    barbarbarbarAlpp ) 
	                         + ( barAlpp - 4.*barbarAlpp + 6.*barbarbarAlpp - 4.*barbarbarbarAlpp + barbarbarbarbarAlpp )     ;  
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
 }
 else if(ADMOrder==5)
 {
 	volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
	volScalarField barbarbarAlpp = convKernel.coarseningSca("barbarbarAlpp",filterWidthADM,barbarAlpp);	
	volScalarField barbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarAlpp",filterWidthADM,barbarbarAlpp);	
	volScalarField barbarbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarbarAlpp",filterWidthADM,barbarbarbarAlpp);	
	volScalarField barbarbarbarbarbarAlpp = convKernel.coarseningSca("barbarbarbarbarbarAlpp",filterWidthADM,barbarbarbarbarAlpp);	
 	volScalarField phiStar =     barAlpp 
	                         + ( barAlpp -    barbarAlpp ) 
				 + ( barAlpp - 2.*barbarAlpp +     barbarbarAlpp ) 
				 + ( barAlpp - 3.*barbarAlpp +  3.*barbarbarAlpp -     barbarbarbarAlpp ) 
	                         + ( barAlpp - 4.*barbarAlpp +  6.*barbarbarAlpp -  4.*barbarbarbarAlpp +    barbarbarbarbarAlpp )     
	                         + ( barAlpp - 5.*barbarAlpp + 10.*barbarbarAlpp - 10.*barbarbarbarAlpp + 5.*barbarbarbarbarAlpp - barbarbarbarbarbarAlpp)     ;  
	barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
 	volScalarField phiStarsqr = phiStar*phiStar;
	barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
 }

  // Predicted scalar variance
  predictedDriftVelocity = barphiStarsqr - barAlpp*(tildeUf-tildeUs);
  predictedDriftVelocity.write();	
*/
