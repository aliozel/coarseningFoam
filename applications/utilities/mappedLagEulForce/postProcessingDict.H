// Post-processing dictionary 
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingPropertiesPostProcessing",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Define couplingProperties dictionary
const dictionary dict(couplingProperties);

// Define postprocessing dictionary
const dictionary postProcessingDict(dict.subDict("EulerianMappingProps"));					

// Filename for dump files
word filename(postProcessingDict.lookup("couplingFilename"));
fileName relativePath(postProcessingDict.lookup("relativePath"));

// DEM timestep
const scalar DEMts(readScalar(postProcessingDict.lookup("DEMts")));

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

// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// OpenFOAM searching algorithm
Switch treeSearch(postProcessingDict.lookup("treeSearch"));
meshSearch searchEngine(U.mesh(),polyMesh::FACEPLANES);

// Max. number of cells per particle
int maxCellPerPart(30);
// Min. alpha for mapping
scalar alphaMin(0.3); 
