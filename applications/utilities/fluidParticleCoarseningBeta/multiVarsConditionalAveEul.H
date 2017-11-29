//- Conditionally averaging eulerian drag coefficient
/*
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(multVariableList[4]+iDirName);
      word varName(multVariableList[1]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeff.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //multipleVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,multUserLoopLabel);			  	
      multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,multUserLoopLabel);			  	
      // Write every timestep
      multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
      multUserLoopLabel++;	
}

//- Conditionally averaging eulerian drag coefficient, Num & Dom
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(multVariableList[4]+iDirName);
      word varName(multVariableList[4]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffNum.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      //multipleVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,multUserLoopLabel);			  	
      multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,multUserLoopLabel);				  	
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

      //word varName(multVariableList[4]+iDirName);
      word varName(multVariableList[5]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffDom.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //multipleVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,multUserLoopLabel);			  	
      multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,multUserLoopLabel);	
      // Write every timestep
      multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
      multUserLoopLabel++;	
}

//- Multi-conditional averaging of vdrift
//for(label iDir=0; iDir<3; iDir++)
label iDir = 2;
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(multVariableList[4]+iDirName);
      word varName(multVariableList[6]+iDirName);
      scalarField scalarVar = barAlppVdrift.internalField().component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //multipleVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,multUserLoopLabel);			  	
      multipleVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,multUserLoopLabel);	
      // Write every timestep
      multipleVarsCondAve.writeUserDefine(varName,caseName,multUserLoopLabel);
      multUserLoopLabel++;	
}
*/
//- Conditionally averaging eulerian drag coefficient
forAll(multVariableList,jj)
{
   word variableName(multVariableList[jj]);
   word varName(filterWidthName);
   multipleVarsCondAve.calc(variableName,varName,multAveragingVariableList,jj+multVariableList.size()*fWidth);
}


