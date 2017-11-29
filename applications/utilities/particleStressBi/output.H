// Write output
writeOutput(	    nParticleClass,
			   runTime,
			      mesh,
		outputRelativePath,
		     outputKinFile,
		    outputCollFile,
			bidisperse,
		     calcCollision,
			  sigmaKin,
			 sigmaColl,
			 domainVol,
		   alppInSubVolume	);

// Normalize kinetic/collisional stress in each cell
forAll(mesh.cells(),cellI)
{
    for(int ii=0; ii<6; ii++ )
    {
         sigmaKinCell[cellI][ii] /= mesh.V()[cellI];
        sigmaCollCell[cellI][ii] /= mesh.V()[cellI];
    }
}
// Write kinetic/collisional and mean particle velocity in each cell
  meanVelCell.write();
 sigmaKinCell.write();
sigmaCollCell.write();

if( domainAve )
{
	for(int iPartClass = 0; iPartClass <= nParticleClass ; iPartClass++)
	{
		for(label ii=0; ii < 6; ii++)
		{
			globalSigmaColl[iPartClass][ii] = sigmaColl[iPartClass][ii];
			// Parallel computation
			reduce(globalSigmaKin[iPartClass][ii], sumOp<scalar>());
			reduce(globalSigmaColl[iPartClass][ii], sumOp<scalar>());					
		}
	}
}

if(verbose)
{
	if(calcCollision)
	{
		for(label ii=0; ii < exList.size(); ii++)
		{
			Info << "" << endl;
			int exIndex = exList[ii]; 
			Info << "" << endl;
			if ( particleCloud.cellIDs()[exIndex][0] > -1 )
				Info << tab << "Fcoll = " << fcoll[exIndex][0] << " " << fcoll[exIndex][1] << " " << fcoll[exIndex][2] << endl;
		}
	}
	if(cohesion)
	{
		for(label ii=0; ii < exList.size(); ii++)
		{
			Info << "" << endl;
			int exIndex = exList[ii]; 
			if ( particleCloud.cellIDs()[exIndex][0] > -1 )
				Info << " Fcoh = " <<  fcoh[exIndex][0] << " " <<  fcoh[exIndex][1] << " " <<  fcoh[exIndex][2] << endl;
			Info << "" << endl;
		}	
	}	

	int iPartClass = 0;
	Info << " " << endl;
	Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
	Pout << " sigmaKin_xx= " << sigmaKin[iPartClass][0]/domainVol << endl;
	Pout << " sigmaKin_yy= " << sigmaKin[iPartClass][1]/domainVol << endl;
	Pout << " sigmaKin_zz= " << sigmaKin[iPartClass][2]/domainVol << endl;			
	Pout << " sigmaKin_xy= " << sigmaKin[iPartClass][3]/domainVol << endl;
	Pout << " sigmaKin_xz= " << sigmaKin[iPartClass][4]/domainVol << endl;
	Pout << " sigmaKin_yz= " << sigmaKin[iPartClass][5]/domainVol << endl;

	if( domainAve && Pstream::master() )
	{						
		Info << " Domain sigmaKin_xx= " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
		Info << " Domain sigmaKin_yy= " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
		Info << " Domain sigmaKin_zz= " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
		Info << " Domain sigmaKin_xy= " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
		Info << " Domain sigmaKin_xz= " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
		Info << " Domain sigmaKin_yz= " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;													
	}

	if(calcCollision)
	{
		Info << " " << endl;
		Pout << " Particle class = " << iPartClass << " collisional stresses" << endl;
		Pout << " sigmaColl_xx= " << sigmaColl[iPartClass][0]/domainVol << endl;
		Pout << " sigmaColl_yy= " << sigmaColl[iPartClass][1]/domainVol << endl;
		Pout << " sigmaColl_zz= " << sigmaColl[iPartClass][2]/domainVol << endl;						
		Pout << " sigmaColl_xy= " << sigmaColl[iPartClass][3]/domainVol << endl;
		Pout << " sigmaColl_xz= " << sigmaColl[iPartClass][4]/domainVol << endl;			
		Pout << " sigmaColl_yz= " << sigmaColl[iPartClass][5]/domainVol << endl;

		if( domainAve && Pstream::master() )
		{						
			Info << " Domain sigmaColl_xx= " << globalSigmaColl[iPartClass][0]/globalDomainVol << endl;
			Info << " Domain sigmaColl_yy= " << globalSigmaColl[iPartClass][1]/globalDomainVol << endl;
			Info << " Domain sigmaColl_zz= " << globalSigmaColl[iPartClass][2]/globalDomainVol << endl;						
			Info << " Domain sigmaColl_xy= " << globalSigmaColl[iPartClass][3]/globalDomainVol << endl;
			Info << " Domain sigmaColl_xz= " << globalSigmaColl[iPartClass][4]/globalDomainVol << endl;			
			Info << " Domain sigmaColl_yz= " << globalSigmaColl[iPartClass][5]/globalDomainVol << endl;														
		}

	}



	if( bidisperse ) 
	{
		for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
		{
			Info << " " << endl;
			Pout << " Particle class = " << iPartClass << " kinetic stresses" << endl;
			Pout << " sigmaKin_xx= " << sigmaKin[iPartClass][0]/domainVol << endl;
			Pout << " sigmaKin_yy= " << sigmaKin[iPartClass][1]/domainVol << endl;
			Pout << " sigmaKin_zz= " << sigmaKin[iPartClass][2]/domainVol << endl;			
			Pout << " sigmaKin_xy= " << sigmaKin[iPartClass][3]/domainVol << endl;
			Pout << " sigmaKin_xz= " << sigmaKin[iPartClass][4]/domainVol << endl;
			Pout << " sigmaKin_yz= " << sigmaKin[iPartClass][5]/domainVol << endl;

			if( domainAve && Pstream::master() )
			{						
				Info << " Domain sigmaKin_xx= " << globalSigmaKin[iPartClass][0]/globalDomainVol << endl;
				Info << " Domain sigmaKin_yy= " << globalSigmaKin[iPartClass][1]/globalDomainVol << endl;
				Info << " Domain sigmaKin_zz= " << globalSigmaKin[iPartClass][2]/globalDomainVol << endl;			
				Info << " Domain sigmaKin_xy= " << globalSigmaKin[iPartClass][3]/globalDomainVol << endl;
				Info << " Domain sigmaKin_xz= " << globalSigmaKin[iPartClass][4]/globalDomainVol << endl;
				Info << " Domain sigmaKin_yz= " << globalSigmaKin[iPartClass][5]/globalDomainVol << endl;												
			}						

			if(calcCollision)
			{
				Info << " " << endl;
				Pout << " Particle class = " << iPartClass << " collisional stresses" << endl;
				Pout << " sigmaColl_xx= " << sigmaColl[iPartClass][0]/domainVol << endl;
				Pout << " sigmaColl_yy= " << sigmaColl[iPartClass][1]/domainVol << endl;
				Pout << " sigmaColl_zz= " << sigmaColl[iPartClass][2]/domainVol << endl;						
				Pout << " sigmaColl_xy= " << sigmaColl[iPartClass][3]/domainVol << endl;
				Pout << " sigmaColl_xz= " << sigmaColl[iPartClass][4]/domainVol << endl;			
				Pout << " sigmaColl_yz= " << sigmaColl[iPartClass][5]/domainVol << endl;

				if( domainAve && Pstream::master() )
				{						
					Info << " Domain sigmaColl_xx= " << globalSigmaColl[iPartClass][0]/globalDomainVol << endl;
					Info << " Domain sigmaColl_yy= " << globalSigmaColl[iPartClass][1]/globalDomainVol << endl;
					Info << " Domain sigmaColl_zz= " << globalSigmaColl[iPartClass][2]/globalDomainVol << endl;						
					Info << " Domain sigmaColl_xy= " << globalSigmaColl[iPartClass][3]/globalDomainVol << endl;
					Info << " Domain sigmaColl_xz= " << globalSigmaColl[iPartClass][4]/globalDomainVol << endl;			
					Info << " Domain sigmaColl_yz= " << globalSigmaColl[iPartClass][5]/globalDomainVol << endl;														
				}

			}						
		}
	}
}
