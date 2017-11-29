// get path to particle VTK files
string HH = string(filename);
const char* charFilename;
charFilename = HH.c_str();
char index[100];
sprintf(index, charFilename, count);

/*
fileName H
(
   particleCloud.mesh().time().processorCase()
 ? particleCloud.mesh().time().path()/".."/relativePath/index
 : particleCloud.mesh().time().path()/relativePath/index
);
*/

fileName H(relativePath/index);
Info << tab << "Opening file: " << H << endl;

// set file pointer
HH = string(H);
const char * particleFilePath = HH.c_str();
ifstream* inputPtr;
inputPtr = new ifstream(particleFilePath);

// Read data
string just_read = " ";
int just_read_int;
scalar just_read_scalar;
label nP;

for(label ii=0; ii<2; ii++)
{
    *inputPtr >> just_read;   // ITEM: TIMEstep
}    
*inputPtr >> just_read_int;   // number    
for(label ii=0; ii<4; ii++)
{
    *inputPtr >> just_read;   // ITEM: NUMBER OF ATOMS
}    
*inputPtr >> nP;   // number

// give nr of particles to cloud
particleCloud.setNumberOfParticles(nP);

Info << tab << "Number of particles = " << particleCloud.numberOfParticles() << endl;

double ** fieldDummy;
particleCloud.dataExchangeM().allocateArray(fieldDummy,0.,4);

// Allocate vars
particleCloud.dataExchangeM().allocateArray(positions,0.,3);
particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
particleCloud.dataExchangeM().allocateArray(omegas,0.,3);
particleCloud.dataExchangeM().allocateArray(types,0.,1);	
particleCloud.dataExchangeM().allocateArray(radii,0.,1);
particleCloud.dataExchangeM().allocateArray(cellID,-1,particleCloud.voidFractionM().maxCellsPerParticle());

// re-allocate arrays of cloud
//particleCloud.reAllocArrays();
particleCloud.reAllocArraysPost();

for(label ii=0; ii<6; ii++)
{
    *inputPtr >> just_read;   // ITEM: BOX BOUNDS pp pp pp
}
for(label ii=0; ii<2; ii++)
{
    *inputPtr >> just_read_scalar;   // number
}
for(label ii=0; ii<2; ii++)
{
    *inputPtr >> just_read_scalar;   // number
} 
for(label ii=0; ii<2; ii++)
{
    *inputPtr >> just_read_scalar;   // number
}           
for(label ii=0; ii<17; ii++)
{
    *inputPtr >> just_read;   // ITEM: ATOMS id type x y z vx vy vz fx fy fz omegax omegay omegaz radius 
}

for(int index = 0; index < nP; ++index)
{
    *inputPtr >> fieldDummy[index][0]
	      >> types[index][0]
	      >> positions[index][0]
	      >> positions[index][1]
	      >> positions[index][2]
	      >> velocities[index][0]
	      >> velocities[index][1]		  
	      >> velocities[index][2]
	      >> fieldDummy[index][1]
	      >> fieldDummy[index][2]
	      >> fieldDummy[index][3]
	      >> omegas[index][0]
	      >> omegas[index][1]
	      >> omegas[index][2]
	      >> radii[index][0];
}

// Set radii & cellIDs
particleCloud.dataExchangeM().destroy(fieldDummy,4);
