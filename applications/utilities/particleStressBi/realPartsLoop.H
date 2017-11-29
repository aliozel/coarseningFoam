for (int jj = 1; jj < k; jj++) // k=0 particle itself
//for (int index_i = 0; index_i < k; index_i++) // k=0 particle itself ???
{	
        label index_i = particlesInSubDomain[nnIdx[jj]]; //nnIdx[jj];
	
	//Info << " index_i " << index_i << endl;

	// Neighboring particle variables
	velI  = particleCloud.velocity(index_i);
	posI  = particleCloud.position(index_i);
	omeI  = particleCloud.omega(index_i);
	radI  = particleCloud.radius(index_i);
	typeI = particleCloud.type(index_i) - 1; // Just nclass starts from "0"									      

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
