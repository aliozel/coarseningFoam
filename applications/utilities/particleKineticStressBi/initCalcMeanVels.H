// Particle mass + total mass
scalar mass_p(0);
scalar totalmass_p(0);			

// Initiate
int iPartClass = 0;
//npPartClass[iPartClass] = particleCloud.numberOfParticles();
npPartClass[iPartClass] = nPartInSubVolume;  

sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0); 

if( bidisperse ) 
{
	for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	{
		// Initiate particle number of each class
		npPartClass[iPartClass] = 0;

		// Initiate velocities
		for(int idir=0; idir<3; idir++) 
		{
			meanVel[iPartClass][idir] = 0 ;
			if(domainAve) globalMeanVel[iPartClass][idir] = 0 ;
		}

		// Initiate particle kinetic stresses
		sigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
		// Initiate global particle kinetic stresses
		globalSigmaKin[iPartClass] = symmTensor(0,0,0,0,0,0);
	}
}		

// Create particle list
//for(int index = 0; index <  particleCloud.numberOfParticles(); index++)
for(int ii = 0; ii <  particlesInSubDomain.size(); ii++)
{							
	// Index of particle
	label index = particlesInSubDomain[ii];
	// Total velocity of particles
	int iPartClass = 0;	
	//for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
	for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += mass_p * particleCloud.velocity(index)[idir];			
	totalmass_p += mass_p;

        // Total velocity of particles in each cell
        label cellID;
        cellID = particleCloud.cellIDs()[index][0];
        if(cellID>-1) 
	{
		nPCell[cellID] += 1.0;
		for(int idir=0; idir<3; idir++)
                {
                        meanVelCell[cellID][idir] += particleCloud.velocity(index)[idir];
                }

		if( bidisperse )
		{
			iPartClass = particleCloud.type(index);

			for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
			npPartClass[iPartClass]++;	
		}	
	}
}	

// Global domain average
if (domainAve) 
{
	for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	{
		for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] = meanVel[iPartClass][idir];

		// Parallel computation
		reduce(globalMeanVel[iPartClass], sumOp<vector>());
		// Normalize
		//for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] /= particleCloud.numberOfParticles();
                for(int idir=0; idir<3; idir++) globalMeanVel[iPartClass][idir] /= nPartInSubVolume;
        }
}			

// Normalize sub-volume velocities
if(npPartClass[iPartClass]!=0)
{
	 for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
}else
{
	for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] = 0;
}

if( bidisperse )
{
	for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	{
		for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
	}
}

// Normalize velocity at each cell
forAll(mesh.cells(),cellI)
{
    if(nPCell[cellI]!=0)
    {
        for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir]/= nPCell[cellI];
    }else
    {
        for(int idir=0; idir<3; idir++) meanVelCell[cellI][idir] = 0.;
    }
}

if(verboseParRun)
{
	iPartClass = 0;									
	Info << " " << endl;
	Pout << " Particle class = " << iPartClass << endl;
	Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
	Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
	Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

	if (domainAve && Pstream::master())
	{				
		Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
		Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
		Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
	}

	if( bidisperse ) 
	{
		for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        {
			Info << " " << endl;
			Pout << " Particle class = " << iPartClass << endl;
			Pout << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
			Pout << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
			Pout << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

			if (domainAve && Pstream::master())
			{				
				Pout << " Domain <u_p,x> = " << globalMeanVel[iPartClass][0] << endl;
				Pout << " Domain <u_p,y> = " << globalMeanVel[iPartClass][1] << endl;
				Pout << " Domain <u_p,z> = " << globalMeanVel[iPartClass][2] << endl;
			}

		}
	}

}
