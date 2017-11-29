// Post-processing dictionary 
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingPropertiesStressPostProcessing",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Define couplingProperties dictionary
const dictionary dict(couplingProperties);

// Define postprocessing dictionary
const dictionary postProcessingDict(dict.subDict("parcelStressModellingProps"));

// Filename for dump files
word filename(postProcessingDict.lookup("couplingFilename"));
fileName relativePath(postProcessingDict.lookup("relativePath"));

// DEM timestep
const scalar DEMts(readScalar(postProcessingDict.lookup("DEMts")));

// Read particle density
const scalar rhoParticle(readScalar(postProcessingDict.lookup("rhoParticle")));

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

// Define fluid coarsening
bool fluidCoarsening(false);
bool reCalcFilterVar(false);
bool readFilteredVars(false);
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
	if(fluidCoarseningDict.found("reCalcFilterVar")) reCalcFilterVar = true;
	if(fluidCoarseningDict.found("readFilteredVars")) readFilteredVars = true;
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

//- Three-bin conditional averaging
//- Number of variables is 3
bool threeConditionalAveraging(false);
wordList threeVariableList;
label threeNVariable(0);
wordList threeAveragingVariableList;
label threeNAveragingVariable(0);
dictionary threeConditionalAveragingDict;
if(postProcessingDict.found("threeVarsConditionalAve")) 
{
	threeConditionalAveraging = true;
	Info << "Three variables conditional averaging will be performed" << endl;
	threeConditionalAveragingDict = postProcessingDict.subDict("threeVarsConditionalAve");
	threeVariableList = wordList(threeConditionalAveragingDict.lookup("variableList"));
	threeAveragingVariableList = wordList(threeConditionalAveragingDict.lookup("averagingVariableList"));
	threeNVariable = threeVariableList.size();
	threeNAveragingVariable = 3; //averagingVariableList.size();		
	if(threeAveragingVariableList.size()>3)
	{
		Info << "\nThree-variable conditional averaging works for only three avareging variables " << abort(FatalError);
	}
}

//- Number of total case 
int nTotalCase = nVariable * nAveragingVariable;
int multNTotalCase = multNVariable; 
int threeNTotalCase = threeNVariable; 
Info << "Total number of analyses = " << nTotalCase+multNTotalCase +threeNTotalCase << endl;

// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// OpenFOAM searching algorithm
Switch treeSearch(postProcessingDict.lookup("treeSearch"));
//meshSearch searchEngine(U.mesh(),polyMesh::FACEPLANES);
meshSearch searchEngine(mesh,polyMesh::FACEPLANES);

// Interpolation of Us OR Us vs Us*
bool interpolation(false);
if(postProcessingDict.found("interpolation")) 
{
	interpolation = true;
	Info << "Us will be interpolated" << endl;
}
bool useUsstar(false);
if(postProcessingDict.found("useUsstar")) 
{
	useUsstar = true;
	Info << "meanUp will be used instead of Us" << endl;	
}	

// Max. number of cells per particle
int maxCellPerPart(30);
if (postProcessingDict.found("maxCellPerPart"))
{
   maxCellPerPart = readScalar(postProcessingDict.lookup("maxCellPerPart"));
   Info << tab << "Max. number of cells per particle " << maxCellPerPart
             << " for mapping & interpolation" << endl;
}

// Min. alpha for mapping
scalar alphaMin(0.3); 

// Not weighting
bool notWeighthing(false);
if(postProcessingDict.found("notWeighting")) 
{
	notWeighthing = true;
	Info << "\nWeigthing is not used" << endl;
	maxCellPerPart = 1;
}
