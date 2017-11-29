
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
       multipleVarsCondAve.calcUserDefine(scalarVar,dragF,varName,userLoopLabel);			  	
       // Write every timestep
       multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
       multUserLoopLabel++;	
 }

 
