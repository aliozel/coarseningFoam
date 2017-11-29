//- Dummy word
word varName("1X1X1");

//- Conditional averaging
condAve.calc();
condAve.write(varName);

//- Multi-variable conditional averaging
multipleVarsCondAve.calc();
multipleVarsCondAve.write(varName);

//- three-variable conditional averaging
threeVarsCondAve.calc();
threeVarsCondAve.write(varName);
