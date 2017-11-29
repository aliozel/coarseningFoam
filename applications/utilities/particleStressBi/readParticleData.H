// Create particle list
labelList particlesInSubDomain(particleCloud.numberOfParticles(),-1);

// Number of particles in sub-domain
int nPartInSubVolume(0);
scalar totalPartVol(0);	

if ( runTime.timeName() != "0" )
{

	int count = runTime.value() / particleCloud.dataExchangeM().DEMts();				

	//Read only if master
	if(Pstream::master())
	{
		Info<< " " << endl;
		particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
		Info<< tab << "Reading particle velocities" << endl;
		Info<< " " << endl;

		particleCloud.dataExchangeM().getData("omega","vector-atom",omegas,count);
		Info<< tab << "Reading particle angular velocities" << endl;
		Info<< " " << endl;
				
		particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
		Info<< tab << "Reading particle positions" << endl;
		Info<< " " << endl;

		particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
		Info<< tab << "Reading particle radius" << endl;		
		Info<< " " << endl;

		particleCloud.dataExchangeM().getData("type","scalar-atom",types,count);
		Info<< tab << "Reading particle types " << endl;
		Info<< " " << endl;

		if ( verbose )
		{
			particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
			Info<< " 	Reading forces on particles" << endl;
			Info<< " " << endl;
		}
	}
	
	// Send to child CPUs
	if(Pstream::parRun())
	{
		for(int index = 0; index < maxNumberOfParticles; index++)
		{
			for(int idir = 0; idir < 3; idir++)
			{
				Pstream::scatter(velocities[index][idir]);
				Pstream::scatter(positions[index][idir]);
				Pstream::scatter(omegas[index][idir]);
			}
			Pstream::scatter(radii[index][0]);
			Pstream::scatter(types[index][0]);
		}
	}

        particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
	particleCloud.setPos(positions);
	particleCloud.setVel(velocities);
	particleCloud.setOmega(omegas);
	particleCloud.setType(types);	

	// Particle IDs in sub-domain
	for(int index = 0; index < particleCloud.numberOfParticles(); index++)
	{
		label cellI = particleCloud.cellIDs()[index][0];					
		if(cellI > -1)
		{
			totalPartVol += 4./3. * constant::mathematical::pi
					      * particleCloud.radius(index) 
					      * particleCloud.radius(index)
					      * particleCloud.radius(index); 

                    	particlesInSubDomain[nPartInSubVolume] = index;
                    	nPartInSubVolume++;
		}
	}	

	// Resize partID list in sub-domain
        particlesInSubDomain.resize(nPartInSubVolume);
	// Sub-volume of domain
	scalar domainVol(0);	
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Pout << tab << "Sub-domain volume[m^3] = " << domainVol << endl;
	alppInSubVolume = totalPartVol / domainVol;
	Pout << tab << "Number of particles in sub-domain = " << nPartInSubVolume <<endl;
	Pout << tab << "Solid volume fraction in sub-domain = " << alppInSubVolume <<endl;

	if(verbose)
	{
		for(label ii=0; ii < exList.size(); ii++)
		{
			Info << "" << endl;
			int exIndex = exList[ii]; 
			label cellI = particleCloud.cellIDs()[exIndex][0];				
			if(cellI > -1)
			{					
				Pout << tab << "index  = " << exIndex << endl;
				Pout << tab << "rp     = " << particleCloud.radius(exIndex) << endl;
				Pout << tab << "Type   = " << particleCloud.type(exIndex) << endl;
				Pout << tab << "Vp     = " << particleCloud.velocity(exIndex) << endl;
				Pout << tab << "Omegap = " << particleCloud.omega(exIndex) << endl;				
				Pout << tab << "Xp     = " << particleCloud.position(exIndex) << endl;
				Pout << tab << "Fp     = " << forces[exIndex][0] << " " << forces[exIndex][1] << " " << forces[exIndex][2] << endl;
				Pout << tab << "CellID = " << cellI << endl;
			}	
		}
	}

}
			
