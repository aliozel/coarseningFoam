// Read dictionary + dictionaryProps 
const dictionary dict(particleCloud.couplingProperties());
//int maxNumberOfParticles(0);

// oneWayDump dictionary
const dictionary oneWayDumpPropsDict(dict.subDict("oneWayDumpProps"));
//maxNumberOfParticles=readScalar(oneWayDumpPropsDict.lookup("maxNumberOfParticles"));

// Filename for dump files
word filename(oneWayDumpPropsDict.lookup("couplingFilename"));
fileName relativePath(oneWayDumpPropsDict.lookup("relativePath"));

// Read particleStress Sub-dictionary
const dictionary particleStressPropsDict(dict.subDict("particleKineticStressProps"));

// Debuging
bool verbose(false);
// Particle ID for debuging
labelList exList;
if(particleStressPropsDict.found("verbose")) 
{
	verbose = true;	
	//exIndex = readInt(particleStressPropsDict.lookup("exIndex"));	
	exList = labelList(particleStressPropsDict.lookup("exList"));
}


// Number particle class
int nParticleClass = 1;
bool bidisperse(false);
// Bidisperse case
if(particleStressPropsDict.found("nParticleClass"))
{
	nParticleClass = readInt(particleStressPropsDict.lookup("nParticleClass"));
	bidisperse = true;
	Info << " " << endl;
	Pout << "Bi-disperse case, number of particle classes = " << nParticleClass << endl;		
}

// Domain average
bool domainAve(false);
if(particleStressPropsDict.found("domainAve")) domainAve = true;

// Parallel verbose
bool verboseParRun(false);
if(particleStressPropsDict.found("verboseParRun")) verboseParRun = true;

// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// Create local output folder	
fileName outputRelativePath(particleStressPropsDict.lookup("outputRelativePath"));
if( !isDir(mesh.time().path()/outputRelativePath) )
{
	mkDir(mesh.time().path()/outputRelativePath );
}

// Create global output folder	
fileName outputGlobalRelativePath
(
   runTime.processorCase()
 ? ".."/particleStressPropsDict.lookup("outputRelativePath") // 
 : particleStressPropsDict.lookup("outputRelativePath")      // 
);

if( !isDir(mesh.time().path()/outputGlobalRelativePath) )
{
	mkDir(mesh.time().path()/outputGlobalRelativePath );
}
