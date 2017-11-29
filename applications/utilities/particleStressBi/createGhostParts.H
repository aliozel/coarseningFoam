// Allocate ghost particle arrays
//ghostParticleCloud.dataExchangeM().allocateArray(ghostPartPositions,0.,3,particlesInSubDomain.size());	
//ghostParticleCloud.dataExchangeM().allocateArray(ghostPartVelocities,0.,3,particlesInSubDomain.size());
//ghostParticleCloud.dataExchangeM().allocateArray(ghostPartOmegas,0.,3,particlesInSubDomain.size());
//ghostParticleCloud.dataExchangeM().allocateArray(ghostPartRadii,0.,1,particlesInSubDomain.size());
//ghostParticleCloud.dataExchangeM().allocateArray(ghostPartTypes,0.,1,particlesInSubDomain.size());

// Create ghost particle list on boundaries
labelList localGhostPartID(maxNumberOfParticles/nProcs+1,-1);
labelListList ghostPartID(nProcs,localGhostPartID); 
labelListList ghostPartTypes(nProcs,localGhostPartID);

// Create ghost particle scalar & vector list 
List<vector> localGhostPartVec(maxNumberOfParticles/nProcs+1,vector(0,0,0));
List<List <vector> > ghostPartPositions(nProcs,localGhostPartVec);
List<List <vector> > ghostPartVelocities(nProcs,localGhostPartVec);
List<List <vector> > ghostPartOmegas(nProcs,localGhostPartVec);

scalarList localGhostPartSca(maxNumberOfParticles/nProcs+1,scalar(-1));
List<scalarList> ghostPartRadii(nProcs,localGhostPartSca);

// Domain lengths
scalar delx = maxX-minX;
scalar dely = maxY-minY;
scalar delz = maxZ-minZ;

//-Initiate number of particles
int nGhostPart = 0;
		
// Search for particles touching boundaries
//for(int index = 0; index < particleCloud.numberOfParticles(); index++)
for(int ii = 0; ii < particlesInSubDomain.size(); ii++)
{							
	label index = particlesInSubDomain[ii];	

	scalar delxMinSq = (particleCloud.position(index).x() - minX)*(particleCloud.position(index).x() - minX);
	scalar delxMaxSq = (particleCloud.position(index).x() - maxX)*(particleCloud.position(index).x() - maxX);
	scalar delyMinSq = (particleCloud.position(index).y() - minY)*(particleCloud.position(index).y() - minY);
	scalar delyMaxSq = (particleCloud.position(index).y() - maxY)*(particleCloud.position(index).y() - maxY);
	scalar delzMinSq = (particleCloud.position(index).z() - minZ)*(particleCloud.position(index).z() - minZ);
	scalar delzMaxSq = (particleCloud.position(index).z() - maxZ)*(particleCloud.position(index).z() - maxZ);
	scalar radsq = particleCloud.radius(index)*particleCloud.radius(index);

	// Find particle on boundaries
	//if(       ( - minX + particleCloud.position(index).x() ) <= particleCloud.radius(index) )
	if(delxMinSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() + delx;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() ;
				
		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);
		
		if (verboseGhost) Pout << " Particle " << index << " on boundary xmax " << endl; 

		nGhostPart++;		  
	}
	//else  if( (   maxX - particleCloud.position(index).x() ) <= particleCloud.radius(index) )
	if(delxMaxSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() - delx;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() ;

		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);

		if (verboseGhost) Pout << " Particle " << index << " on boundary xmax " << endl; 
		nGhostPart++;

	}
	//else  if( ( - minY + particleCloud.position(index).y() ) <= particleCloud.radius(index) )
	if(delyMinSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() + dely;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() ;						

		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);

		if (verboseGhost) Pout << " Particle " << index << " on boundary ymin " << endl; 
		nGhostPart++;

	}	
	//else if( (   maxY - particleCloud.position(index).y() ) <= particleCloud.radius(index) )
	if(delyMaxSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() - dely;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() ;

		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);

		if (verboseGhost) Pout << " Particle " << index << " on boundary ymax " << endl;
		nGhostPart++;

	}	 
	//else if( ( - minZ + particleCloud.position(index).z() ) < particleCloud.radius(index) )
	if(delzMinSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() + delz;					

		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);

		if (verboseGhost) Pout << " Particle " << index << " on boundary zmin " << endl;
		nGhostPart++;

	}	 
	//else if( (   maxZ - particleCloud.position(index).z() ) <= particleCloud.radius(index) )
	if(delzMaxSq <= radsq)
	{
		ghostPartID[Pstream::myProcNo()][nGhostPart] = index;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][0] = particleCloud.position(index).x() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][1] = particleCloud.position(index).y() ;
		ghostPartPositions[Pstream::myProcNo()][nGhostPart][2] = particleCloud.position(index).z() - delz;

		for(int idir=0; idir<3; idir++) 
		{
		  ghostPartVelocities[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.velocity(index)[idir];
		  ghostPartOmegas[Pstream::myProcNo()][nGhostPart][idir] = particleCloud.omega(index)[idir];
		}
		
		ghostPartRadii[Pstream::myProcNo()][nGhostPart] = particleCloud.radius(index);
		ghostPartTypes[Pstream::myProcNo()][nGhostPart] = particleCloud.type(index);

		if (verboseGhost) Pout << " Particle " << index << " on boundary zmax " << endl; 						
		nGhostPart++;
	}		            							      
}

// Parallel computation
// Resize number of ghost particle array
ghostPartID[Pstream::myProcNo()].resize(nGhostPart);

// Distrubute to all processors
Pstream::gatherList(ghostPartID);
Pstream::scatterList(ghostPartID);

// Create global ghost particle IDs
labelList globalGhostPartID = ListListOps::combine<labelList>
(
	ghostPartID,
	accessOp<labelList>()
);

// Number of ghost particles
Pout << tab << "Number of ghost particles = " << globalGhostPartID.size() << endl;

// Ghost particle types
ghostPartPositions[Pstream::myProcNo()].resize(nGhostPart);
ghostPartVelocities[Pstream::myProcNo()].resize(nGhostPart);
ghostPartOmegas[Pstream::myProcNo()].resize(nGhostPart);
ghostPartRadii[Pstream::myProcNo()].resize(nGhostPart);
ghostPartTypes[Pstream::myProcNo()].resize(nGhostPart);

// Distrubute to all processors
Pstream::gatherList(ghostPartPositions);
Pstream::gatherList(ghostPartVelocities);
Pstream::gatherList(ghostPartOmegas);
Pstream::gatherList(ghostPartRadii);
Pstream::gatherList(ghostPartTypes);
Pstream::scatterList(ghostPartPositions);
Pstream::scatterList(ghostPartVelocities);
Pstream::scatterList(ghostPartOmegas);
Pstream::scatterList(ghostPartRadii);
Pstream::scatterList(ghostPartTypes);

// Create global ghost particle variables
List<vector> globalGhostPartPositions = ListListOps::combine< List<vector> >
(
	ghostPartPositions,
	accessOp< List<vector> >()
);
List<vector> globalGhostPartVelocities = ListListOps::combine< List<vector> >
(
	ghostPartVelocities,
	accessOp< List<vector> >()
);
List<vector> globalGhostPartOmegas = ListListOps::combine< List<vector> >
(
	ghostPartOmegas,
	accessOp< List<vector> >()
);
scalarList globalGhostPartRadii = ListListOps::combine<scalarList>
(
	ghostPartRadii,
	accessOp<scalarList>()
);
labelList globalGhostPartTypes = ListListOps::combine<labelList>
(
	ghostPartTypes,
	accessOp<labelList>()
);

/*
for(int index = 0; index < globalGhostPartID.size(); index++)
	Info << " index = " << index << " types = " << globalGhostPartTypes[index] << endl;
*/

//
if(verboseGhost)
{
	Pout << " Ghost ID " << tab 
	     << "Position[x]" << tab << "Real pos.[x]" << tab
	     << "Position[y]" << tab << "Real pos.[y]" << tab
	     << "Position[z]" << tab << "Real pos.[z]" << endl;;
	for(int index = 0; index < globalGhostPartID.size(); index++)
	{
		Pout << index << tab;
		Pout << globalGhostPartID[index] << tab;
		for(int idir=0; idir<3; idir++)
		{
			Pout << ghostPartPositions[index][idir] << tab;
			Pout << particleCloud.position(globalGhostPartID[index])[idir] <<tab;
		}
		Pout << endl;
	}	
}
			
