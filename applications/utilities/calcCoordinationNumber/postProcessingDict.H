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
const dictionary postProcessingDict(dict.subDict("coordinationNumberProps"));					

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

// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// OpenFOAM searching algorithm
Switch treeSearch(postProcessingDict.lookup("treeSearch"));
meshSearch searchEngine(mesh,polyMesh::FACEPLANES);

// Search radius
scalar searchRadius(1.2);
if(postProcessingDict.found("searchRadius")) searchRadius = readScalar(postProcessingDict.lookup("searchRadius"));

// Neighboring algorithm class
// Number of particles in a fixed-radius
int k = 20; //20; 

// Dimensions, exact OR approximate  
int dim =3; double eps = 0;

// Squared radius distance
double sqRad;
