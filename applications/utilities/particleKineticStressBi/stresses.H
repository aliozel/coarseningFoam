// Index of particle
label index_j = particlesInSubDomain[ii];

// Type of particle
label typeJ = particleCloud.velocity(index_j)[0];
	
// Particle mass
mass_p = 4./3.*rhoParticle*constant::mathematical::pi*particleCloud.radius(index_j)*particleCloud.radius(index_j)*particleCloud.radius(index_j);	

// Particle kinetic stress
iPartClass = 0;

sigmaKin[iPartClass][0] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] )
                                    * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] );

sigmaKin[iPartClass][1] += mass_p   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] )
                                    * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

sigmaKin[iPartClass][2] += mass_p   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] )
                                    * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

sigmaKin[iPartClass][3] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] )
                                    * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

sigmaKin[iPartClass][4] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] )
                                    * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

sigmaKin[iPartClass][5] += mass_p   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] )
                                    * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

label cellID;
cellID = particleCloud.cellIDs()[index_j][0];
if(cellID>-1)
{
	// Particle kinetic stress in each cell
	sigmaKinCell[cellID][0] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVelCell[cellID][0] )
                                	    * ( particleCloud.velocity(index_j).x() - meanVelCell[cellID][0] );

	sigmaKinCell[cellID][1] += mass_p   * ( particleCloud.velocity(index_j).y() - meanVelCell[cellID][1] )
                                	    * ( particleCloud.velocity(index_j).y() - meanVelCell[cellID][1] );

	sigmaKinCell[cellID][2] += mass_p   * ( particleCloud.velocity(index_j).z() - meanVelCell[cellID][2] )
                                	    * ( particleCloud.velocity(index_j).z() - meanVelCell[cellID][2] );

	sigmaKinCell[cellID][3] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVelCell[cellID][0] )
                                	    * ( particleCloud.velocity(index_j).y() - meanVelCell[cellID][1] );

	sigmaKinCell[cellID][4] += mass_p   * ( particleCloud.velocity(index_j).x() - meanVelCell[cellID][0] )
                                	    * ( particleCloud.velocity(index_j).z() - meanVelCell[cellID][2] );

	sigmaKinCell[cellID][5] += mass_p   * ( particleCloud.velocity(index_j).y() - meanVelCell[cellID][1] )
                                	    * ( particleCloud.velocity(index_j).z() - meanVelCell[cellID][2] );

}

if(domainAve)
{
	 globalSigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
						 * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ); 

	 globalSigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
						 * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ); 

	 globalSigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ) 
						 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ); 

	 globalSigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
						 * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] );

	 globalSigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
						 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );

	 globalSigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
						 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );				
}				  

if( bidisperse ) 
{				
	 iPartClass = typeJ + 1;

	 sigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ); 

	 sigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
					   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ); 

	 sigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ) 
					   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ); 

	 sigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

	 sigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

	 sigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
					   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );	

	 if(domainAve)
	 {
		 globalSigmaKin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					        	 * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ); 

		 globalSigmaKin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
							 * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ); 

		 globalSigmaKin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ) 
							 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] ); 

		 globalSigmaKin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					        	 * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] );

		 globalSigmaKin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - globalMeanVel[iPartClass][0] ) 
					        	 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );

		 globalSigmaKin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - globalMeanVel[iPartClass][1] ) 
							 * ( particleCloud.velocity(index_j).z() - globalMeanVel[iPartClass][2] );						

	}

}	
