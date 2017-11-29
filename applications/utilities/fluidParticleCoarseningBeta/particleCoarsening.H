 // Calculate particle resolved drag force
 Info << "" << endl;
 Info << tab << "Calculating particle resolved drag force..." << endl;
 //resolvedDragF().setForce(barAlpf,tildeUf,nP[me],cellIDs,positions[me],velocities[me],radii[me]);
 resolvedDragF().setForce(barAlpf,tildeUf,nP[me],cellIDs,positions[me],velocities[me],radii[me],particleWeights,neighboringCellIDs);	

 // Calculate particle resolved gas pressure force
 Info << tab << "Calculating particle resolved gas pressure force..." << endl;
 //resolvedDragF().setForce(tildeGradPg,nP[me],cellIDs,positions[me],radii[me]);
 resolvedDragF().setForce(tildeGradPg,nP[me],cellIDs,positions[me],radii[me],particleWeights,neighboringCellIDs);

 label cellI(-1);
 vector position(0,0,0);

 // Lagrangian variables definition
 #include "lagrangianVars.H"			    

 //for(label index = 0; index < particlesInSubDomain.size(); index++)
 for(label index = 0; index < nP[me]; index++)
 {
      for(label iDir=0; iDir<3; iDir++)
      {		

	    // Ratio between filtered drag coefficient and resolved drag coeff.
	    normalizedParticleDragCoeff[index][iDir] = dragF().Fp(index)[iDir];

	    if(filteringWithGradPg) 
	       normalizedParticleDragCoeff[index][iDir] += (           dragF().gradPgFp(index)[iDir]
				                             - resolvedDragF().gradPgFp(index)[iDir] ); 

   	    // Avoid some numerical issues of conditional averaging
	    normalizedParticleDragCoeffNum[index][iDir] = normalizedParticleDragCoeff[index][iDir];
	    normalizedParticleDragCoeffDom[index][iDir] = resolvedDragF().Fp(index)[iDir]; 	

	    if( resolvedDragF().Fp(index)[iDir] != 0 ) 
	       normalizedParticleDragCoeff[index][iDir] /= resolvedDragF().Fp(index)[iDir]; 	

	    // DNS & Resolved Lagrangian relative velocity 		
	    uSlipLag[index][iDir] =	dragF().Ufxp(index)[iDir]; // - dragF().Up(index)[iDir];
	    resolvedUSlipLag[index][iDir] =	resolvedDragF().Ufxp(index)[iDir]; // - dragF().Up(index)[iDir];

	    // Drag correction gij
	    cellI = cellIDs[index];
 	    position[0] = positions[me][index][0];
 	    position[1] = positions[me][index][1];
 	    position[2] = positions[me][index][2];

	    if (cellI>-1 && resolvedUSlipLag[index][iDir] > SMALL && resolvedDragF().voidfractionxp(index) <1 ) 
		gij[index][iDir] =  (   barAlppUfxp.interpolate(position,cellI)[iDir] 
			          - (1.- resolvedDragF().voidfractionxp(index) ) * resolvedDragF().Ufxp(index)[iDir]  )
			       /( (1.- resolvedDragF().voidfractionxp(index)) * ( resolvedUSlipLag[index][iDir] ) );

	    if( gij[index][iDir] >  0 ) gij[index][iDir] =  0.;	
	    if( gij[index][iDir] < -1 ) gij[index][iDir] = -1.;

	    // Drift velocity
	    if (cellI>-1) alppVdi[index][iDir] =  barAlppUfxp.interpolate(position,cellI)[iDir] 
			          - (1.- resolvedDragF().voidfractionxp(index)) * resolvedDragF().Ufxp(index)[iDir];  

	    // Unresolved drag * uSlip 
	    unresolvedDrag[index]  = dragF().Fp(index) - resolvedDragF().Fp(index);

	    // grad(barAlpp)@p
	    if (cellI>-1) gradBarAlppxp[index] = IntGradBarAlppxp.interpolate(position,cellI);

	    // Resolved relaxation time 			     
	    resolvedTaup[index] = resolvedDragF().taup(index);

	    // Resolved U@xp-Us@xp
	    if (cellI>-1) diffResolvedUUsxp[index] = resolvedDragF().Ufxp(index) - tildeUsxp.interpolate(position,cellI);

	    // gradTildeUf@xp
	    //gradTildeUfzxp[index] = IntGradTildeUfzxp.interpolate(position,cellI);

	    // Filtered eulerian force@xp
	    if (cellI>-1) barForceEulerianxp[index] = IntbarForceEulerianxp.interpolate(position,cellI);

	    // Resolved drag
	    resolvedLagDrag[index] = resolvedDragF().Fp(index);

	    // Solid fraction @xp
	    baralppxp[index] = (1.- resolvedDragF().voidfractionxp(index));
      }

      /*
      if(index<5) 
      {
	     Info << " normalizedParticleDragCoeff[index][iDir] = " << normalizedParticleDragCoeff[index] << endl;
	     Info << " Particle force = " << dragF().Fp(index) << endl;
      }
      */
 }

 // Write unresolved drag * uSlip 
 //normalizedParticleDragCoeff.write();
 //normalizedParticleDragCoeffNum.write();
 //normalizedParticleDragCoeffDom.write();
 //unresolvedDrag.write();
 //uSlipLag.write();
 //resolvedUSlipLag.write();
 //gradBarAlppxp.write();
 //resolvedTaup.write();
 //diffResolvedUUsxp.write();
 //gradTildeUfzxp.write();
 //barForceEulerianxp.write();
 //resolvedLagDrag.write();
 //baralppxp.write();

 // Eulerian filtering- //Defined in lagrangianVars
 //vectorField normalizedEulerianDragCoeff(mesh.cells().size(),vector(0,0,0));
 // Avoid some numerical issues of conditional averaging
 //vectorField normalizedEulerianDragCoeffNum(mesh.cells().size(),vector(0,0,0));
 //vectorField normalizedEulerianDragCoeffDom(mesh.cells().size(),vector(0,0,0));
