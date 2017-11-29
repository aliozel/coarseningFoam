 //- Conditinally averaging the ratio between filtered drag coefficient and resolved drag coeff. 	
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[0]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeff.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       //condAve.calcUserDefine(scalarVar,resolvedDragF,varName,userLoopLabel);	
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel,particleWeights,neighboringCellIDs,maxCellPerPart);
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;	
 }
 
 //- Conditinally averaging the ratio between filtered drag coefficient and resolved drag coeff. Num & Dom	 
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[2]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeffNum.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       //condAve.calcUserDefine(scalarVar,resolvedDragF,varName,userLoopLabel);	
       condAve.calcUserDefine(scalarVar,resolvedDragF,varName,userLoopLabel,particleWeights,neighboringCellIDs,maxCellPerPart);
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

       word varName(variableList[3]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeffDom.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       //condAve.calcUserDefine(scalarVar,resolvedDragF,varName,userLoopLabel);			  	
       condAve.calcUserDefine(scalarVar,resolvedDragF,varName,userLoopLabel,particleWeights,neighboringCellIDs,maxCellPerPart);
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;	
 } 

/*  //- Conditionally averaging lagrangian slip velocity
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[1]+iDirName);
       scalarField scalarVar = uSlipLag.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 } 

 //- Conditionally averaging gij
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[2]+iDirName);
       scalarField scalarVar = gij.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 } 

 //- Conditionally averaging drift velocity
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[3]+iDirName);
       scalarField scalarVar = alppVdi.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 } 
*/
/*
 //- Conditionally averaging eulerian slip velocity
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[5]+iDirName);
       scalarField scalarVar = uSlipEulerian.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 } 

 //- Conditionally averaging eulerian gij
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[6]+iDirName);
       scalarField scalarVar = gijEulerian.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 } 

 //- Conditionally averaging eulerian drift velocity
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(variableList[7]+iDirName);
       scalarField scalarVar = alppVdiEulerian.component(iDir);	
       word caseName(filterWidthName);
       // Conditionally averaging 			  			  
       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       condAve.writeUserDefine(varName,caseName,userLoopLabel);
       userLoopLabel++;
 }  */

/* 		    //- Conditionally averaging user define variables
 for(int ii=3; ii <= variableList.size(); ii++)
 {
	 for(label iDir=0; iDir<3; iDir++)
	 {
	       word iDirName;
	       if(iDir==0){iDirName="X";}
	       else if(iDir==1){iDirName="Y";}
	       else if(iDir==2){iDirName="Z";}

	       word varName(variableList[ii]+iDirName);
	       scalarField scalarVar = normalizedParticleDragCoeff.component(iDir);	
	       word caseName(filterWidthName);
	       label userLoopLabel = iDir + fWidth*3 + fWidth*nParcel*3 + nAveragingVariable*6;
	       //- Debugging			  
	       //Info << " iDir " << iDir << " fWidth " << fWidth << " userLoopLabel " << userLoopLabel << endl;

  	       condAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
	       // Write at the last time step
	       condAve.writeUserDefine(varName,caseName,userLoopLabel);
	 }
 }	 */   
