 // Loop over ghost particles 	
for (int index_i = 0; index_i < globalGhostPartID.size(); index_i++) //
{

	// Neighboring particle variables	
	posI = globalGhostPartPositions[index_i];
	velI = globalGhostPartVelocities[index_i];
        omeI = globalGhostPartOmegas[index_i];
	radI = globalGhostPartRadii[index_i];
	typeI = globalGhostPartTypes[index_i] - 1;

	//Info << " Ghost loop, index_j " << index_j << " index_i " << index_i << endl;

	// Calculate collision forces
	calcForce( index_j,velJ,posJ,omeJ,radJ,typeJ,
		   index_i,velI,posI,omeI,radI,typeI,
		   collisionModelI,
		   youngsModulus,
		   poissonsRatio,
		   coefficientRestitution,
		   coefficientFriction,
   	           k_n, 
		   k_t, 
		   gamma_n, 
		   gamma_t, 
                   e_n,
                   mu_f,		
  		   fcoll, 
		   ftan,				
		   dt, 
		   tangential_history, 
		   liq, 
		   liquid_transfer, 
		   liquidVol, 
		   surf_tension, 
		   fluid_visc, 
		   fcap, 
		   fvisc,
		   first_touch,
		   cohesion,
		   minimumDistanceVdW,
		   cohEnergyDens,
		   fcoh,
		   rhop,
		   bidisperse,
		   sigmaCollJI,
		   verboseColl,
		   numberOfPairs			);
		   	
	iPartClass = 0;
        sigmaColl[iPartClass][0] += sigmaCollJI.xx();
        sigmaColl[iPartClass][1] += sigmaCollJI.yy();
        sigmaColl[iPartClass][2] += sigmaCollJI.zz();
        sigmaColl[iPartClass][3] += sigmaCollJI.xy();
        sigmaColl[iPartClass][4] += sigmaCollJI.xz();
        sigmaColl[iPartClass][5] += sigmaCollJI.yz();

	if( bidisperse )
	{							  						
                if( typeJ == typeI )
                {

			iPartClass = typeJ + 1 ;
			sigmaColl[iPartClass][0] += sigmaCollJI.xx();
                        sigmaColl[iPartClass][1] += sigmaCollJI.yy();
                        sigmaColl[iPartClass][2] += sigmaCollJI.zz();
                        sigmaColl[iPartClass][3] += sigmaCollJI.xy();
                        sigmaColl[iPartClass][4] += sigmaCollJI.xz();
                        sigmaColl[iPartClass][5] += sigmaCollJI.yz();
                }

	}
}		 

