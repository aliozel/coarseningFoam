//- Dummy word
word varName(filterWidthName);

//- Conditional averaging
condAve.calc();
condAve.write(varName);

//- Multi-variable conditional averaging
multipleVarsCondAve.calc();
multipleVarsCondAve.write(varName);

//- three-variable conditional averaging
threeVarsCondAve.calc();
threeVarsCondAve.write(varName);
