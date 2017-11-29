//- Loop over variables
forAll(variableList,jj)
{
   word variableName(variableList[jj]);	
   word varName(filterWidthName);
   condAve.calc(variableName,varName,averagingVariableList,jj+variableList.size()*fWidth); 
}

forAll(multVariableList,jj)
{
   word variableName(multVariableList[jj]);	
   word varName(filterWidthName);
   multipleVarsCondAve.calc(variableName,varName,multAveragingVariableList,jj+multVariableList.size()*fWidth); 
}

forAll(threeVariableList,jj)
{
   word variableName(threeVariableList[jj]);	
   word varName(filterWidthName);
   threeVarsCondAve.calc(variableName,varName,threeAveragingVariableList,jj+threeVariableList.size()*fWidth);    
}    
