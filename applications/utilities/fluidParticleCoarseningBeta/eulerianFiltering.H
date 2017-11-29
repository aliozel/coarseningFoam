 // Calculate Eulerian resolved drag force
 Info << tab << "Calculating Eulerian resolved drag force..." << endl;
 resolvedDragFEulerian().setForceEulerian(barAlpf,tildeUf,radii[me],tildeUs);

 // Calculate Eulerian resolved gas pressure force
 Info << tab << "Calculating Eulerian resolved gas pressure force..." << endl;
 //resolvedDragFEulerian().setForceEulerian(barAlpf,tildeGradPg);
 resolvedDragFEulerian().setForceEulerian(1.-barAlpf,tildeGradPg);	

 vectorField uSlipEulerian(mesh.cells().size(),vector(0,0,0));				  
 vectorField gijEulerian(mesh.cells().size(),vector(0,0,0));		    
 vectorField alppVdiEulerian(mesh.cells().size(),vector(0,0,0));		    			  		

 forAll(mesh.cells(),index)
 {	
      for(label iDir=0; iDir<3; iDir++)
      {		
	    // Ratio between filtered drag coefficient and resolved drag coeff.
	   //if( resolvedDragFEulerian().FpEulerian(index)[iDir] > SMALL ) 
	   {			      
	     /*
	     if( index < 5 )
	     {
	       //Info << " Mesh = " << dragFEulerian().cellID(index) << " alpp = " << (1.-voidfraction[index])
	       Info << " Mesh = " << dragFEulerian().cellID(index) << " alpp = " << (1.-mappedVoidfraction[index])
		    //<< " resolvedDragF().gradPgFpEulerian = " << -resolvedDragFEulerian().gradPgFpEulerian(index)[iDir] 
		    //<< " Ksl*(U-Us) = " << Ksl[index] * ( U[index][iDir] - Us[index][iDir] )		       
		    //<< " dragF().gradPgFpEulerian = " << dragFEulerian().gradPgFpEulerian(index)[iDir]
		    << " resolvedDragF().FpEulerian = " << resolvedDragFEulerian().FpEulerian(index)[iDir]
		    //<< " dragF().FpEulerian = " << dragFEulerian().FpEulerian(index)[iDir] 
		    << " barForceEulerian = " << barForceEulerian[index][iDir]
		    << " dragFpEulerian = " << dragFpEulerian[index][iDir] 
		    << endl;					    

	       //     << " Particle force = " << dragF().Fp(particleIDinParcel) << endl;
	       //Info << tab << " resolvedDragF().Fp(particleIDinParcel) = " << resolvedDragF().Fp(particleIDinParcel) << endl;
	       //Info << tab << "                dragFParcel().Fp(index) = " << dragFParcel().Fp(index) << endl;
	     }
	     */			      

	    if( resolvedDragFEulerian().FpEulerian(index)[iDir] != 0 ) 
	       //normalizedEulerianDragCoeff[index][iDir] =  ( barForceEulerian[index][iDir] / mesh.V()[index] )
	       //					    / resolvedDragFEulerian().FpEulerian(index)[iDir];	

	       normalizedEulerianDragCoeff[index][iDir] =  barForceEulerian[index][iDir] 
				  			 / resolvedDragFEulerian().FpEulerian(index)[iDir];	

	       normalizedEulerianDragCoeffNum[index][iDir] = barForceEulerian[index][iDir];
	       normalizedEulerianDragCoeffDom[index][iDir] = resolvedDragFEulerian().FpEulerian(index)[iDir];	

	    /*
	    if( index < 5 )
	    {
	       //Info << " alpp = " << (1.-voidfraction[index])
	       Info << " alpp = " << (1.-mappedVoidfraction[index])
		    << " normalizedEulerianDragCoeff = " << normalizedEulerianDragCoeff[index][iDir] << endl;
	    } 
	    */


	   }

	    // Lagrangian relative velocity 		
	    uSlipEulerian[index][iDir] = tildeUf[index][iDir] - tildeUs[index][iDir];

	    // Drag correction gij
	    if ( ( tildeUf[index][iDir] - tildeUs[index][iDir] ) > SMALL &&  barAlpp[index] > SMALL )
	    gijEulerian[index][iDir] =  ( barAlppUf[index][iDir] - ( barAlpp[index] * tildeUf[index][iDir] ) )
			               /(   barAlpp[index] * ( tildeUf[index][iDir] - tildeUs[index][iDir] ) );

	    if( gijEulerian[index][iDir] >  0 ) gijEulerian[index][iDir] =  0.;	
	    if( gijEulerian[index][iDir] < -1 ) gijEulerian[index][iDir] = -1.;

	    // Drift velocity
	    alppVdiEulerian[index][iDir] = barAlppUf[index][iDir] - ( barAlpp[index] * tildeUf[index][iDir] );

      }		    	
 }

 normalizedEulerianDragCoeff.write();
 normalizedEulerianDragCoeffNum.write();	
 normalizedEulerianDragCoeffDom.write();	

 //- Relaxation time @filter size & @testfiltersize
 volScalarField taup
 (
     IOobject
     (
         "taup",
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
     taup[cellI] = resolvedDragFEulerian().taupEulerian(cellI);
 }
 taup.write();
 //- Writing taup with filter size
 volScalarField taupFilter
 (
     IOobject
     (
         "taup"+filterWidthName,
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     taup
 );
 taupFilter.write();  

 // Dynamic procedure
 if(dynamicAdjustment)
 {
	 if(dynamicAdjustmentVar=="scalarVariance") 
	 {
	    //- Dynamic adjusment; scalar variance
	    #include "dynamicAdjustmentScalarVariance.H"
	 }else if (dynamicAdjustmentVar=="driftVelocity") 
	 {
	    //- Dynamic adjusment; drift velocity
	    #include "dynamicAdjustmentDriftVelocity.H"
         }else
	 {
	    Info << " There is no such variable for dynamic adjusment " << endl;	
	 } 
 }
