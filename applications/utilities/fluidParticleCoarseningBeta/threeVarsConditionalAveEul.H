//- Conditionally averaging eulerian drag coefficient
/*
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(threeVariableList[4]+iDirName);
      word varName(threeVariableList[1]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeff.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //threeVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,threeUserLoopLabel);			  	
      threeVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,threeUserLoopLabel);			  	
      // Write every timestep
      threeVarsCondAve.writeUserDefine(varName,caseName,threeUserLoopLabel);
      threeUserLoopLabel++;	

      Info << "safdg sag a" << endl;
}

//- Conditionally averaging eulerian drag coefficient, Num & Dom
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(threeVariableList[4]+iDirName);
      word varName(threeVariableList[4]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffNum.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      //threeVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,threeUserLoopLabel);			  	
      threeVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,threeUserLoopLabel);				  	
      // Write every timestep
      threeVarsCondAve.writeUserDefine(varName,caseName,threeUserLoopLabel);
      threeUserLoopLabel++;	
}

for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(threeVariableList[4]+iDirName);
      word varName(threeVariableList[5]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffDom.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //threeVarsCondAve.calcUserDefine(scalarVar,dragFEulerian,varName,threeUserLoopLabel);			  	
      threeVarsCondAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,threeUserLoopLabel);	
      // Write every timestep
      threeVarsCondAve.writeUserDefine(varName,caseName,threeUserLoopLabel);
      threeUserLoopLabel++;	
}
*/

forAll(threeVariableList,jj)
{
   word variableName(threeVariableList[jj]);
   word varName(filterWidthName);
   threeVarsCondAve.calc(variableName,varName,threeAveragingVariableList,jj+threeVariableList.size()*fWidth);
}
