//- Conditionally averaging eulerian drag coefficient
/*
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(variableList[4]+iDirName);
      word varName(variableList[1]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeff.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //condAve.calcUserDefine(scalarVar,dragFEulerian,varName,userLoopLabel);			  	
      condAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,userLoopLabel);			  	
      // Write every timestep
      condAve.writeUserDefine(varName,caseName,userLoopLabel);
      userLoopLabel++;	
}

//- Conditionally averaging eulerian drag coefficient, Num & Dom
for(label iDir=0; iDir<3; iDir++)
{
      word iDirName;
      if(iDir==0){iDirName="X";}
      else if(iDir==1){iDirName="Y";}
      else if(iDir==2){iDirName="Z";}

      //word varName(variableList[4]+iDirName);
      word varName(variableList[4]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffNum.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      //condAve.calcUserDefine(scalarVar,dragFEulerian,varName,userLoopLabel);			  	
      condAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,userLoopLabel);				  	
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

      //word varName(variableList[4]+iDirName);
      word varName(variableList[5]+iDirName);
      scalarField scalarVar = normalizedEulerianDragCoeffDom.component(iDir);	
      word caseName(filterWidthName+charFilteringGradPg);
      // Conditionally averaging 
      //condAve.calcUserDefine(scalarVar,dragFEulerian,varName,userLoopLabel);			  	
      condAve.calcUserDefine(scalarVar,resolvedDragFEulerian,varName,userLoopLabel);	
      // Write every timestep
      condAve.writeUserDefine(varName,caseName,userLoopLabel);
      userLoopLabel++;	
}
*/
//- Loop over variables
forAll(variableList,jj)
{
   word variableName(variableList[jj]);
   word varName(filterWidthName);
   condAve.calc(variableName,varName,averagingVariableList,jj+variableList.size()*fWidth);
}

