// Post-processing dictionary 
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingPropertiesPostProcessingScalarVariance",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Define couplingProperties dictionary
const dictionary dict(couplingProperties);

// Define postprocessing dictionary
const dictionary postProcessingDict(dict.subDict("postProcessingProps"));					

// Filename for dump files
word filename(postProcessingDict.lookup("couplingFilename"));
fileName relativePath(postProcessingDict.lookup("relativePath"));

// Fluid coarsening
bool fluidCoarsening(false);
bool readFilteredVars(true);
labelList filterWidth(1,1);
dictionary fluidCoarseningDict;
label nFilter(0);
if(postProcessingDict.found("fluidCoarsening")) 
{
	fluidCoarsening = true;
	Info << "\nFluid coarsening will be performed" << endl;
	fluidCoarseningDict = postProcessingDict.subDict("fluidCoarsening");
	filterWidth = labelList(fluidCoarseningDict.lookup("filterWidth"));
	nFilter = filterWidth.size();
	//if(fluidCoarseningDict.found("readFilteredVars")) readFilteredVars = true;
}


//- Number of variables & number of conditional averaging variables
bool conditionalAveraging(false);
wordList variableList;
label nVariable(0);
wordList averagingVariableList;
label nAveragingVariable(0);
dictionary conditionalAveragingDict;
if(postProcessingDict.found("conditionalAve")) 
{
	conditionalAveraging = true;
	Info << "Conditinal averaging will be performed" << endl;
	conditionalAveragingDict = postProcessingDict.subDict("conditionalAve");
	variableList = wordList(conditionalAveragingDict.lookup("variableList"));
	averagingVariableList = wordList(conditionalAveragingDict.lookup("averagingVariableList"));
	nVariable = variableList.size();
	nAveragingVariable = averagingVariableList.size();
}

//- Double-bin conditional averaging
//- Number of variables is 2
bool multConditionalAveraging(false);
wordList multVariableList;
label multNVariable(0);
wordList multAveragingVariableList;
label multNAveragingVariable(0);
dictionary multConditionalAveragingDict;
if(postProcessingDict.found("multiVarsConditionalAve")) 
{
	multConditionalAveraging = true;
	Info << "Two-variables conditional averaging will be performed" << endl;
	multConditionalAveragingDict = postProcessingDict.subDict("multiVarsConditionalAve");
	multVariableList = wordList(multConditionalAveragingDict.lookup("variableList"));
	multAveragingVariableList = wordList(multConditionalAveragingDict.lookup("averagingVariableList"));
	multNVariable = multVariableList.size();
	multNAveragingVariable = 2; //averagingVariableList.size();		
	if(multAveragingVariableList.size()>2)
	{
		Info << "\nMulti-variable conditional averaging works for only two avareging variables " << abort(FatalError);
	}
}

//- Number of total case ( iDir=3 ) + only fluid coarsening
int nTotalCase = 3 * nFilter * nVariable * nAveragingVariable;
int multNTotalCase = 3 * nFilter * multNVariable; 
Info << "Total number of analyses = " << nTotalCase+multNTotalCase << endl;




