volScalarField predictedScalarVariance
(
    IOobject
    (
        "predictedADMScalarVariance"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);

volScalarField computedScalarVariance
(
    IOobject
    (
        "computedScalarVariance"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);

// Calculate computed scalar variance
volScalarField barAlpp = convKernel.coarseningSca("barAlpp",filterWidthADM,alpp);		
volScalarField barAlppAlpp = convKernel.coarseningSca("barAlppAlpp",filterWidthADM,alpp*alpp);		
computedScalarVariance = barAlppAlpp - barAlpp*barAlpp;
computedScalarVariance.write();

volScalarField barphiStar
(
    IOobject
    (
        "barphiStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);

volScalarField barphiStarsqr
(
    IOobject
    (
        "barphiStarsqr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);

//- Drift velocity
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

volVectorField tildeUsStar
(
    IOobject
    (
        "tildeUsStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
);

volVectorField tildeUfStar
(
    IOobject
    (
        "tildeUfStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
);

volVectorField bardriftVelStar
(
    IOobject
    (
        "bardriftVelStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
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
       volScalarField phiStar = barAlpp;
       barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
       volScalarField phiStarsqr = phiStar*phiStar;
       barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
       
       volVectorField UsStar = tildeUs;
       volVectorField UfStar = tildeUf;
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar);              		
}
else if(ADMOrder==1)
{
       volScalarField barbarAlpp = convKernel.coarseningSca("barbarAlpp",filterWidthADM,barAlpp);	
       volScalarField phiStar =     barAlpp 
	                        + ( barAlpp - barbarAlpp );
       barphiStar = convKernel.coarseningSca("barphiStar",filterWidthADM,phiStar);		
       volScalarField phiStarsqr = phiStar*phiStar;
       barphiStarsqr = convKernel.coarseningSca("barphiStarsqr",filterWidthADM,phiStarsqr);
       
       volVectorField tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidthADM,barAlpp,tildeUs);	
       volVectorField tildetildeUf = convKernel.coarseningVecFavre("tildetildeUf",filterWidthADM,(1.-barAlpp),tildeUf);
       volVectorField UsStar = 	    tildeUs
       	                        + ( tildeUs - tildetildeUs );
       volVectorField UfStar = 	    tildeUf
       	                        + ( tildeUf - tildetildeUf );
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar);              		       
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
       
       volVectorField tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidthADM,barAlpp,tildeUs);	
       volVectorField tildetildeUf = convKernel.coarseningVecFavre("tildetildeUf",filterWidthADM,(1.-barAlpp),tildeUf);
       volVectorField tildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildeUs",filterWidthADM,barbarAlpp,tildetildeUs);	
       volVectorField tildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildeUf",filterWidthADM,(1.-barbarAlpp),tildetildeUf);       
       volVectorField UsStar = 	    tildeUs
       	                        + ( tildeUs - tildetildeUs )
				+ ( tildeUs - 2.*tildetildeUs + tildetildetildeUs );				
       volVectorField UfStar = 	    tildeUf
       	                        + ( tildeUf - tildetildeUf )
				+ ( tildeUf - 2.*tildetildeUf + tildetildetildeUf );				
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar);                 
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
       
       volVectorField tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidthADM,barAlpp,tildeUs);	
       volVectorField tildetildeUf = convKernel.coarseningVecFavre("tildetildeUf",filterWidthADM,(1.-barAlpp),tildeUf);
       volVectorField tildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildeUs",filterWidthADM,barbarAlpp,tildetildeUs);	
       volVectorField tildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildeUf",filterWidthADM,(1.-barbarAlpp),tildetildeUf);       
       volVectorField tildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildeUs",filterWidthADM,barbarbarAlpp,tildetildetildeUs);	
       volVectorField tildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildeUf",filterWidthADM,(1.-barbarbarAlpp),tildetildetildeUf);  
       volVectorField UsStar = 	    tildeUs
       	                        + ( tildeUs -    tildetildeUs )
				+ ( tildeUs - 2.*tildetildeUs +    tildetildetildeUs )				
				+ ( tildeUs - 3.*tildetildeUs + 3.*tildetildetildeUs - tildetildetildetildeUs ) ;
       volVectorField UfStar = 	    tildeUf
       	                        + ( tildeUf -    tildetildeUf )
				+ ( tildeUf - 2.*tildetildeUf +    tildetildetildeUf )				
				+ ( tildeUf - 3.*tildetildeUf + 3.*tildetildetildeUf - tildetildetildetildeUf ) ;
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar); 
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

       volVectorField tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidthADM,barAlpp,tildeUs);	
       volVectorField tildetildeUf = convKernel.coarseningVecFavre("tildetildeUf",filterWidthADM,(1.-barAlpp),tildeUf);
       volVectorField tildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildeUs",filterWidthADM,barbarAlpp,tildetildeUs);	
       volVectorField tildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildeUf",filterWidthADM,(1.-barbarAlpp),tildetildeUf);       
       volVectorField tildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildeUs",filterWidthADM,barbarbarAlpp,tildetildetildeUs);	
       volVectorField tildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildeUf",filterWidthADM,(1.-barbarbarAlpp),tildetildetildeUf);  
       volVectorField tildetildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildetildeUs",filterWidthADM,barbarbarbarAlpp,tildetildetildetildeUs);	
       volVectorField tildetildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildetildeUf",filterWidthADM,(1.-barbarbarbarAlpp),tildetildetildetildeUf);  
       volVectorField UsStar = 	    tildeUs
       	                        + ( tildeUs -    tildetildeUs )
				+ ( tildeUs - 2.*tildetildeUs +    tildetildetildeUs )				
				+ ( tildeUs - 3.*tildetildeUs + 3.*tildetildetildeUs -    tildetildetildetildeUs ) 
	                        + ( tildeUs - 4.*tildetildeUs + 6.*tildetildetildeUs - 4.*tildetildetildetildeUs + tildetildetildetildetildeUs )     ;  
       volVectorField UfStar = 	    tildeUf
       	                        + ( tildeUf -    tildetildeUf )
				+ ( tildeUf - 2.*tildetildeUf +    tildetildetildeUf )				
				+ ( tildeUf - 3.*tildetildeUf + 3.*tildetildetildeUf -    tildetildetildetildeUf ) 
	                        + ( tildeUf - 4.*tildetildeUf + 6.*tildetildetildeUf - 4.*tildetildetildetildeUf + tildetildetildetildetildeUf )     ;  
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar); 
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

       volVectorField tildetildeUs = convKernel.coarseningVecFavre("tildetildeUs",filterWidthADM,barAlpp,tildeUs);	
       volVectorField tildetildeUf = convKernel.coarseningVecFavre("tildetildeUf",filterWidthADM,(1.-barAlpp),tildeUf);
       volVectorField tildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildeUs",filterWidthADM,barbarAlpp,tildetildeUs);	
       volVectorField tildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildeUf",filterWidthADM,(1.-barbarAlpp),tildetildeUf);       
       volVectorField tildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildeUs",filterWidthADM,barbarbarAlpp,tildetildetildeUs);	
       volVectorField tildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildeUf",filterWidthADM,(1.-barbarbarAlpp),tildetildetildeUf);  
       volVectorField tildetildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildetildeUs",filterWidthADM,barbarbarbarAlpp,tildetildetildetildeUs);	
       volVectorField tildetildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildetildeUf",filterWidthADM,(1.-barbarbarbarAlpp),tildetildetildetildeUf);  
       volVectorField tildetildetildetildetildetildeUs = convKernel.coarseningVecFavre("tildetildetildetildetildetildeUs",filterWidthADM,barbarbarbarAlpp,tildetildetildetildetildeUs);	
       volVectorField tildetildetildetildetildetildeUf = convKernel.coarseningVecFavre("tildetildetildetildetildetildeUf",filterWidthADM,(1.-barbarbarbarAlpp),tildetildetildetildetildeUf); 
       volVectorField UsStar = 	    tildeUs
       	                        + ( tildeUs -    tildetildeUs )
				+ ( tildeUs - 2.*tildetildeUs +     tildetildetildeUs )				
				+ ( tildeUs - 3.*tildetildeUs +  3.*tildetildetildeUs -     tildetildetildetildeUs ) 
	                        + ( tildeUs - 4.*tildetildeUs +  6.*tildetildetildeUs -  4.*tildetildetildetildeUs +    tildetildetildetildetildeUs )      
	                        + ( tildeUs - 5.*tildetildeUs + 10.*tildetildetildeUs - 10.*tildetildetildetildeUs + 5.*tildetildetildetildetildeUs - tildetildetildetildetildetildeUs)     ;  
       volVectorField UfStar = 	    tildeUf
       	                        + ( tildeUf -    tildetildeUf )
				+ ( tildeUf - 2.*tildetildeUf +     tildetildetildeUf )				
				+ ( tildeUf - 3.*tildetildeUf +  3.*tildetildetildeUf -     tildetildetildetildeUf ) 
	                        + ( tildeUf - 4.*tildetildeUf +  6.*tildetildetildeUf -  4.*tildetildetildetildeUf +    tildetildetildetildetildeUf )       
	                        + ( tildeUf - 5.*tildetildeUf + 10.*tildetildetildeUf - 10.*tildetildetildetildeUf + 5.*tildetildetildetildetildeUf - tildetildetildetildetildetildeUf)     ;  
       volVectorField driftVelStar = phiStar*(UfStar-UsStar);
       bardriftVelStar = convKernel.coarseningVec("bardriftVelStar",filterWidthADM,driftVelStar); 
}

// Predicted scalar variance
predictedScalarVariance = barphiStarsqr - barphiStar*barphiStar;
predictedScalarVariance.write();	


// Predicted drift velocity
predictedDriftVelocity = bardriftVelStar - barAlpp*(tildeUf-tildeUs);
predictedDriftVelocity.write();	
