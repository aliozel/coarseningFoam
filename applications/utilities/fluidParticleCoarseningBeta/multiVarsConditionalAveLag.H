  //- Conditinally averaging the ratio between filtered drag coefficient and resolved drag coeff. 	
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(multVariableList[0]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeff.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //multipleVarsCondAve.calcUserDefine(scalarVar,dragF,varName,multUserLoopLabel);			  	
       multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragF,varName,multUserLoopLabel);	
       // Write every timestep
       multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
       multUserLoopLabel++;	
 }

 //- Conditinally averaging the ratio between filtered drag coefficient and resolved drag coeff. Num & Dom	 
 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(multVariableList[2]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeffNum.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //multipleVarsCondAve.calcUserDefine(scalarVar,dragF,varName,multUserLoopLabel);			  	
       multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragF,varName,multUserLoopLabel);	
       // Write every timestep
       multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
       multUserLoopLabel++;	
 } 

 for(label iDir=0; iDir<3; iDir++)
 {
       word iDirName;
       if(iDir==0){iDirName="X";}
       else if(iDir==1){iDirName="Y";}
       else if(iDir==2){iDirName="Z";}

       word varName(multVariableList[3]+iDirName);
       scalarField scalarVar = normalizedParticleDragCoeffDom.component(iDir);	
       word caseName(filterWidthName+charFilteringGradPg);
       // Conditionally averaging 
       //multipleVarsCondAve.calcUserDefine(scalarVar,dragF,varName,multUserLoopLabel);			  	
       multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragF,varName,multUserLoopLabel);			  	
       // Write every timestep
       multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
       multUserLoopLabel++;	
 } 
