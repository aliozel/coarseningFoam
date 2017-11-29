//- Dummy word
//word varName(filterWidthName);

//- Conditional averaging
//condAveCoarsening.calc();
//condAveCoarsening.write(varName);

//- Multi-variable conditional averaging
//multipleVarsCondAve.calc();
//multipleVarsCondAve.write(varName);

//- three-variable conditional averaging
//threeVarsCondAve.calc();
//threeVarsCondAve.write(varName);

//- Loop over variables
forAll(variableListCoarsening,jj)
{
   word variableName(variableListCoarsening[jj]);	
   word varName(filterWidthName);
   condAveCoarsening.calc(variableName,varName,averagingVariableListCoarsening,jj+variableListCoarsening.size()*fWidth); 
}

forAll(multVariableListCoarsening,jj)
{
   word variableName(multVariableListCoarsening[jj]);	
   word varName(filterWidthName);
   multipleVarsCondAveCoarsening.calc(variableName,varName,multAveragingVariableListCoarsening,jj+multVariableListCoarsening.size()*fWidth); 
}

forAll(threeVariableListCoarsening,jj)
{
   word variableName(threeVariableListCoarsening[jj]);	
   word varName(filterWidthName);
   threeVarsCondAveCoarsening.calc(variableName,varName,threeAveragingVariableListCoarsening,jj+threeVariableListCoarsening.size()*fWidth);    
}    
