// Read dictionary + dictionaryProps 
const dictionary dict(particleCloud.couplingProperties());
int maxNumberOfParticles(0);

// oneWayVTK dictionary
const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));
maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));

// Read particleStress Sub-dictionary
const dictionary particleStressPropsDict(dict.subDict("particleStressProps"));
const scalar rhop(readScalar(particleStressPropsDict.lookup("rhoParticle")));					

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

// Collision part debugging
bool verboseColl(false);
if(particleStressPropsDict.found("verboseColl")) 
{
	verboseColl = true;			
}

// Ghost part debugging
bool verboseGhost(false);
if(particleStressPropsDict.found("verboseGhost")) 
{
	verboseGhost = true;			
}	

// Parallel print out part debugging
bool verboseParRun(false);
if(particleStressPropsDict.found("verboseParRun")) 
{
	verboseGhost = true;			
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

// Collision dictionary
scalar collisionDp(1.05);
bool calcCollision(false);
if(particleStressPropsDict.found("calcCollision"))
{
	calcCollision=true;
	if(particleStressPropsDict.found("collisionDp"))
	{
		collisionDp=readScalar(particleStressPropsDict.lookup("collisionDp"));
	}
}

int neighboringListParameters(5);
if(particleStressPropsDict.found("neighboringListParameter"))
{
	neighboringListParameters=readInt(particleStressPropsDict.lookup("neighboringListParameter"));
}	

// Domain average
bool domainAve(false);
if(particleStressPropsDict.found("domainAve")) domainAve = true;


// Parallel computation
int nProcs(1);
if(Pstream::parRun()) nProcs = Pstream::nProcs();

// Create output folder	
fileName outputRelativePath(particleStressPropsDict.lookup("outputRelativePath"));
if( !isDir(mesh.time().path()/outputRelativePath) )
{
	mkDir(mesh.time().path()/outputRelativePath );
}
