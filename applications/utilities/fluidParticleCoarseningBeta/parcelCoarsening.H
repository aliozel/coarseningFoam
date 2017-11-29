// Create parcel
forAll(nPInParcel,parI)
{

 //- Create name for parcel 			  
 sprintf(charName, "nP%d",nPInParcel[parI]);
 word parcelName(charName);

 // Create parcelCloud
 parcelCloud parcels(parcelCoarseningDict,mesh,nP[me],positions[me],velocities[me],radii[me],nPInParcel[parI]);					  
 if( readParcels )
 {
   Info << "" << endl;
   Info << tab << "Parcels are reading..." << endl;
   parcels.readParcel();
 }else
 {
   Info << "" << endl;
   Info << tab << "\nParcels are creating..." << endl;
   parcels.createParcel();	
 }	
 // Create parcel resolved drag model class
 dragFParcel = dragModelPost::New(postProcessingDict,mesh);

 Info << tab << "Calculating parcel resolved drag force..." << endl;					  
 dragFParcel().setForceParcel(barAlpf,tildeUf,parcels);

 // Unresolved part of particle&parcel  drag force
 vectorField  normalizedParcelDragCoeff(parcels.numberOfParticles(),vector(0,0,0));
 vectorField  normalizedParcelDragCoeffNum(parcels.numberOfParticles(),vector(0,0,0));
 vectorField  normalizedParcelDragCoeffDom(parcels.numberOfParticles(),vector(0,0,0));

 //
 vector sumResolvedDrag(0,0,0);
 for(label index = 0; index < parcels.numberOfParticles(); index++)
 {
      // parcels.particlesInParcel(index) gives global particleID list
      for(label parI = 0; parI < parcels.particlesInParcel(index).size(); parI++)			       
      {	    
	    label particleIDinParcel = parcels.particlesInParcel(index)[parI];	

	    /*
	    if( index < 5 )
	    {
	      Info << " Parcel = " << index << " particle in parcel " << particleIDinParcel 
		   << " particle position = " << dragF().Up(particleIDinParcel) 		       
		   << " parcel position = " << dragFParcel().Up(index) 
		   << " particle cellIDs = " << dragF().cellID(particleIDinParcel) 		       
		   << " parcel cellIDs = " << dragFParcel().cellID(index) << endl;					    

	      //     << " Particle force = " << dragF().Fp(particleIDinParcel) << endl;
	      //Info << tab << " resolvedDragF().Fp(particleIDinParcel) = " << resolvedDragF().Fp(particleIDinParcel) << endl;
	      //Info << tab << "                dragFParcel().Fp(index) = " << dragFParcel().Fp(index) << endl;
	    }
	    */

	    for(label iDir=0; iDir<3; iDir++)
	    {
	       normalizedParcelDragCoeff[index][iDir] += dragF().Fp(particleIDinParcel)[iDir]; 

	       if(filteringWithGradPg) 
		 normalizedParcelDragCoeff[index][iDir]+= (           dragF().gradPgFp(particleIDinParcel)[iDir]
					                    - resolvedDragF().gradPgFp(particleIDinParcel)[iDir]  ) ;	


	    }

      }
      for(label iDir=0; iDir<3; iDir++)
      {
	  normalizedParcelDragCoeffNum[index][iDir] = normalizedParcelDragCoeff[index][iDir];
	  normalizedParcelDragCoeffDom[index][iDir] = ( parcels.particlesInParcel(index).size() * dragFParcel().Fp(index)[iDir] ) ;

	  if( dragFParcel().Fp(index)[iDir] != 0 )
	    normalizedParcelDragCoeff[index][iDir] /= ( parcels.particlesInParcel(index).size() * dragFParcel().Fp(index)[iDir] ) ;			       
      }

 }

 Info << " userLoopLabel = " << userLoopLabel << endl;

 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[0]+iDirName);
       scalarField scalarVar = normalizedParcelDragCoeff.component(iDir);	
       word caseName(filterWidthName+parcelName+charFilteringGradPg);						  			  
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragFParcel,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
}

for(label iDir=0; iDir<3; iDir++)
{
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[0]+iDirName);
       scalarField scalarVar = normalizedParcelDragCoeffNum.component(iDir);	
       word caseName(filterWidthName+"Num"+parcelName+charFilteringGradPg);						  			  
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragFParcel,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
}

for(label iDir=0; iDir<3; iDir++)
{
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[0]+iDirName);
       scalarField scalarVar = normalizedParcelDragCoeffDom.component(iDir);	
       word caseName(filterWidthName+"Dom"+parcelName+charFilteringGradPg);						  			  
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragFParcel,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
}
}
