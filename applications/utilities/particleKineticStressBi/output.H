// Write local output
writeOutput(	    nParticleClass,
			   runTime,
			      mesh,
		outputRelativePath,
		     outputKinFile,
			bidisperse,
			  sigmaKin,
			 domainVol,
		   alppInSubVolume	);

// Normalize kinetic/collisional stress in each cell
forAll(mesh.cells(),cellI)
{
    for(int ii=0; ii<6; ii++ )
    {
         sigmaKinCell[cellI][ii] /= mesh.V()[cellI];
    }
}
// Write kinetic/collisional and mean particle velocity in each cell
meanVelCell.write();
sigmaKinCell.write();
nPCell.write();

if( domainAve )
{
	if(Pstream::parRun())
	{
	    for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	    {
		    for(label ii=0; ii < 6; ii++)
		    {
			    // Parallel computation
			    reduce(globalSigmaKin[iPartClass][ii], sumOp<scalar>());
		    }
	    }
	
	}
	
	if(Pstream::parRun()&&Pstream::master())
	{
	    // Write global output
	    writeOutput(	      nParticleClass,
			        	     runTime,
			        		mesh,
			    outputGlobalRelativePath,
				 outputGlobalKinFile,
					  bidisperse,
				      globalSigmaKin,
				     globalDomainVol,
			       globalAlppInSubVolume	);
	}
}

if(verbose)
{	

	int iPartClass = 0;
	Info << " " << endl;
	Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
	Pout << " sigmaKin_xx = " << sigmaKin[iPartClass][0]/domainVol << endl;
	Pout << " sigmaKin_yy = " << sigmaKin[iPartClass][1]/domainVol << endl;
	Pout << " sigmaKin_zz = " << sigmaKin[iPartClass][2]/domainVol << endl;			
	Pout << " sigmaKin_xy = " << sigmaKin[iPartClass][3]/domainVol << endl;
	Pout << " sigmaKin_xz = " << sigmaKin[iPartClass][4]/domainVol << endl;
	Pout << " sigmaKin_yz = " << sigmaKin[iPartClass][5]/domainVol << endl;

	if( domainAve && Pstream::master() )
	{						
		Info << " Domain sigmaKin_xx = " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
		Info << " Domain sigmaKin_yy = " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
		Info << " Domain sigmaKin_zz = " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
		Info << " Domain sigmaKin_xy = " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
		Info << " Domain sigmaKin_xz = " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
		Info << " Domain sigmaKin_yz = " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;													
	}

	if( bidisperse ) 
	{
		for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
		{
			Info << " " << endl;
			Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
			Pout << " sigmaKin_xx = " << sigmaKin[iPartClass][0]/domainVol << endl;
			Pout << " sigmaKin_yy = " << sigmaKin[iPartClass][1]/domainVol << endl;
			Pout << " sigmaKin_zz = " << sigmaKin[iPartClass][2]/domainVol << endl;			
			Pout << " sigmaKin_xy = " << sigmaKin[iPartClass][3]/domainVol << endl;
			Pout << " sigmaKin_xz = " << sigmaKin[iPartClass][4]/domainVol << endl;
			Pout << " sigmaKin_yz = " << sigmaKin[iPartClass][5]/domainVol << endl;

			if( domainAve && Pstream::master() )
			{						
				Info << " Domain sigmaKin_xx = " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
				Info << " Domain sigmaKin_yy = " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
				Info << " Domain sigmaKin_zz = " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
				Info << " Domain sigmaKin_xy = " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
				Info << " Domain sigmaKin_xz = " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
				Info << " Domain sigmaKin_yz = " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;												
			}											
		}
	}
}
