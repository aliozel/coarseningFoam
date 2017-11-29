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
const dictionary postProcessingDict(dict.subDict("postProcessingProps"));					

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
// Pressure term will be accounted or not for filtering
bool filteringWithGradPg(true);
//bool filteringWithGradPg(false);
word charFilteringWithGradPg("");
word charFilteringWithoutGradPg("WithoutGradPg");
word charFilteringGradPg;
if(postProcessingDict.found("filteringWithoutGradPg"))
{
        filteringWithGradPg = false;
	Info << "\nCoarsening will be performed without accounting for gradPg" << endl;
	charFilteringGradPg = charFilteringWithoutGradPg; 
}else
{
	charFilteringGradPg = charFilteringWithGradPg;
	Info << "\nCoarsening will be performed with accounting for gradPg" << endl;
}

// Fluid coarsening
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

// Particle coarsening
bool particleCoarsening(false);
dictionary particleCoarseningDict;
if(postProcessingDict.found("particleCoarsening"))
{
	particleCoarsening = true;
	Info << "Particle coarsening will be performed" << endl;
	particleCoarseningDict = postProcessingDict.subDict("particleCoarsening");
}

// Parcel coarsening
bool parcelCoarsening(false);
labelList nPInParcel(1,1);
dictionary parcelCoarseningDict;
label nParcel(0);
bool readParcels(false);
if(postProcessingDict.found("parcelCoarsening"))
{
	parcelCoarsening = true;
	Info << "Particle coarsening will be performed" << endl;
	parcelCoarseningDict = postProcessingDict.subDict("parcelCoarsening");
	nPInParcel = labelList(parcelCoarseningDict.lookup("numberOfParticlesInParcel"));
	nParcel = nPInParcel.size();
	if(particleCoarseningDict.found("readParcels")) readParcels = true;
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
		Info << "\nMulti-variable conditional averaging works for only two averaging variables " << abort(FatalError);
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
		Info << "\nThree-variable conditional averaging works for only three averaging variables " << abort(FatalError);
	}
}

//- Number of total case 
//int nTotalCase = nFilter * ( nVariable * nAveragingVariable );
//int multNTotalCase = nFilter * multNVariable; 
//int threeNTotalCase = nFilter * threeNVariable; 
//Info << "Total number of analyses = " << nTotalCase+multNTotalCase+threeNTotalCase << endl;

//- Number of total case ( iDir=3 ) + only fluid coarsening
/*
int nTotalCase = 3 * ( nFilter * nParcel + nFilter ) * nVariable * nAveragingVariable;
int multNTotalCase = 3 * ( nFilter * nParcel + nFilter ) * multNVariable; 
int threeNTotalCase = 3 * ( nFilter * nParcel + nFilter ) * threeNVariable; 
*/
int nTotalCase = ( nFilter * nParcel + nFilter ) * nVariable * nAveragingVariable;
int multNTotalCase = ( nFilter * nParcel + nFilter ) * multNVariable; 
int threeNTotalCase = ( nFilter * nParcel + nFilter ) * threeNVariable;
Info << "Total number of analyses = " << nTotalCase + multNTotalCase + threeNTotalCase << endl;

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

// Weighting
Switch weighthing(false);
postProcessingDict.lookup("weighthing") >> weighthing;
dictionary weighthingDict;
Switch IDWeighthing(false);
Switch GaussWeighthing(false);
bool noWeighthing(false);

if(weighthing)
{
	weighthingDict =  postProcessingDict.subDict("weighthingProps");
	weighthingDict.lookup("IDWeighthing")    >> IDWeighthing;
	weighthingDict.lookup("GaussWeighthing") >> GaussWeighthing;	

	if(IDWeighthing)
	{
		Info << "\nIDWeighthing is used" << endl;
	}
	if(GaussWeighthing)
	{
		Info << "\nGaussian weigthing is used" << endl;
	}	
}
if(!weighthing)
{
	Info << "\nWeigthing is not used" << endl;
	noWeighthing = true;
	maxCellPerPart = 1;
}

// Use local variables to calculate Eulerian force
bool usingLocalEul(false);
if(postProcessingDict.found("usingLocalEul")) 
{
	usingLocalEul = true;
	Info << "\nUsing local variables for filtered Eulerian force" << endl;
}


// Dynamic adjustment of model parameter
bool dynamicAdjustment(false);
word dynamicAdjustmentVar("");
dictionary dynamicAdjustmentDict;
wordList dynamicAdjustmentProposedFunctionList;
wordList dynamicAdjustmentResolvedFunctionList;
int testFilterWidth(0);
if(postProcessingDict.found("dynamicAdjustment")) 
{
	dynamicAdjustment = true;
	Info << "\nDynamic adjustment will be applied for the model constant" << endl;
	dynamicAdjustmentDict = postProcessingDict.subDict("dynamicAdjustmentProps");
	dynamicAdjustmentDict.lookup("variable") >> dynamicAdjustmentVar;
	dynamicAdjustmentProposedFunctionList = wordList(dynamicAdjustmentDict.lookup("proposedFunction"));
	dynamicAdjustmentResolvedFunctionList = wordList(dynamicAdjustmentDict.lookup("resolvedFunction"));
	dynamicAdjustmentDict.lookup("testFilterWidth") >> testFilterWidth;
}

