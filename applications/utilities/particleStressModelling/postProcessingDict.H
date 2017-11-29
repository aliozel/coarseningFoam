// Post-processing dictionary 

// Define couplingProperties dictionary
const dictionary dict(particleCloud.couplingProperties());

// Define oneWayVTK dictionary
const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));

// Define postprocessing dictionary
//const dictionary postProcessingDict(dict.subDict("postProcessingProps"));					
const dictionary postProcessingDict(dict.subDict("particleStressModellingProps"));

// Maximum number of particles read from oneWayVTK dictionary
int maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));

// Read particle density
const scalar rhop(readScalar(postProcessingDict.lookup("rhoParticle")));

// Read particle diameter
const scalar dp(readScalar(postProcessingDict.lookup("dp")));

// Debugging example particles
bool verbose(false);
labelList exList(1,1);
if(postProcessingDict.found("verbose"))
{
        verbose = true;
        exList = labelList(postProcessingDict.lookup("exList"));
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

//- Number of total case 
int nTotalCase = nVariable * nAveragingVariable;
int multNTotalCase = multNVariable; 
Info << "Total number of analyses = " << nTotalCase+multNTotalCase << endl;
