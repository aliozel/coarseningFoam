#include "calcCollisionalForce.H"

namespace Foam
{

	// Calculate liquid volumes in liquid bridges
	void calcLiqIn(scalarField& liq, int& i, int& j, vectorField& liquidVol)
	{
        	liquidVol[i][j] = liq[i] / 6. + liq[j] / 6.;
        	liq[i] -= liq[i] / 6.;
        	liq[j] -= liq[j] / 6.;
	}

	// Re-distribute liquid to particles
	void calcLiqBack(scalarField& liq, int& i, int& j, vectorField& liquidVol)
	{
        	liq[i] += liquidVol[i][j]/2.;
        	liq[j] += liquidVol[i][j]/2.;
        	liquidVol[i][j] = 0;
	}

	void calcForce(	int& i, vector& velJ, vector& posJ, vector& omeJ, scalar& radJ, int& typeJ, 
			int& j, vector& velI, vector& posI, vector& omeI, scalar& radI, int& typeI, 	
			word& collisionModelI,
			scalarField& youngsModulus,
			scalarField& poissonsRatio,
			vectorField& coefficientRestitution,
			vectorField& coefficientFriction,
			scalar& k_n, 
			scalar& k_t, 
			scalar& gamma_n, 
			scalar& gamma_t,
			scalar& e_n,		  
			scalar& mu_f,
			vectorField& fcoll, 
			vectorField& ftan,
			const scalar& dt, 
			bool& tangential_history, 
			scalarField& liq, 
			bool& liquid_transfer, 
			vectorField& liquidVol, 
			scalar& surf_tension, 
			scalar& fluid_visc, 
			vectorField& fcap, 
			vectorField& fvisc,
			bool& first_touch,
			bool& cohesion,
			scalar& minimumDistanceVdW,
			vectorField& cohEnergyDens,
			vectorField& fcoh,
			const scalar& rhop,	   	    					        
			const bool& bidisperse,
			symmTensor& sigma_coll_JI,
			const bool& verbose,
			scalar& numberOfPairs		)
	{	
		
		const double pi = constant::mathematical::pi; 
		
		double ddist(0);

		scalar n_vec[3];
		scalar t_vec[3];

		scalar dist[3];
		scalar delta_n(0);

		scalar v_r[3];
		scalar v_n[3];
		scalar v_t[3];

		scalar w_r[3];
		scalar v_trel[3];

		scalar F_n[3];
		scalar F_t[3];	
		
		scalar delta_t[3];

		double shrmag(0); double fs(0); double fmu(0);

		// Liquid bridge variables
		double tilde_h(0);
		double tilde_liquidVol(0); 
		double const_A(0);
		double const_B(0); 
		double const_C(0); 
		double Ca(0);
		double ksi(0);

		// Cohesion variables
		double ri(0); double rj(0); 
		double smin(0); double s(0); double smax(0);

		// Collisional stress
		sigma_coll_JI = symmTensor(0,0,0,0,0,0);
		
		// Initialise
		for (int dir=0;dir<3;dir++)
		{				 			 
			F_n[dir] = 0;
			F_t[dir] = 0;
			n_vec[dir] = 0;
			t_vec[dir] = 0;

			dist[dir] = 0;

			v_r[dir] = 0;
			v_n[dir] = 0;
			v_t[dir] = 0;

			w_r[dir] = 0;
			v_trel[dir] = 0;
			
			delta_t[dir] = 0;
		}
		
		for (int dir=0;dir<3;dir++) 
			dist[dir] = posJ[dir]-posI[dir]; 

		double distsq = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
		ddist = sqrt( distsq );

		for (int dir=0;dir<3;dir++) 
			n_vec[dir] = dist[dir] / ( ddist + 1.e-25 );
								
        	delta_n = ( radI + radJ ) - ddist;
        
		for (int dir=0;dir<3;dir++) 
			v_r[dir] = velJ[dir]-velI[dir];

		for (int dir=0;dir<3;dir++) 
			v_n[dir] = ( v_r[0] * n_vec[0] + v_r[1] * n_vec[1] + v_r[2] * n_vec[2] ) * n_vec[dir]; 	
			
		double radsum = radI + radJ;
		double radsumsq = radsum*radsum;

		// If there is collision					
		if( distsq<=radsumsq )
		{	
			numberOfPairs++;

			double mi = 4./3.*pi*rhop*radI*radI*radI;
        		double mj = 4./3.*pi*rhop*radJ*radJ*radJ;

       			double meff = mi * mj / ( mi + mj );

			if(collisionModelI=="hertz")
        		{
				e_n = coefficientRestitution[typeJ][typeI];
				mu_f = coefficientFriction[typeJ][typeI];

				double R_star = radI *radJ / ( radI + radJ ); 
				double Y_star = 1. / ( (1-pow(poissonsRatio[typeI],2))/youngsModulus[typeI] + (1-pow(poissonsRatio[typeJ],2))/youngsModulus[typeJ] );
				double G_star = 1. / ( 2*(2+poissonsRatio[typeI])*(1-poissonsRatio[typeI])/youngsModulus[typeI] + 2*(2+poissonsRatio[typeJ])*(1-poissonsRatio[typeJ])/youngsModulus[typeJ] );
				double Beta = Foam::log(e_n)/Foam::sqrt(Foam::log(e_n)*Foam::log(e_n)+pi*pi);

				double sqrtval = Foam::sqrt(R_star*mag(delta_n));

                		double S_n = 2.*Y_star*sqrtval;
                		double S_t = 8.*G_star*sqrtval;

                		k_n = 4./3.*Y_star*sqrtval;
                		k_t = S_t; //8*G_star*sqrt(R_star*mag(delta_n));

                		gamma_n = -2.*Foam::sqrt(5./6.)*Beta*Foam::sqrt(S_n*meff);
                		gamma_t = -2.*Foam::sqrt(5./6.)*Beta*Foam::sqrt(S_t*meff);


        		 }
        		 else if (collisionModelI=="hooke")
        		 {
                		  gamma_n *= meff;
                		  gamma_t *= meff;
        		 }

			if(verbose)
			{
				Pout << " " << endl;
                       		Pout << " k_n = " << k_n;
                        	Pout << " gamma_n = " << gamma_n;
                                Pout << " k_t = " << k_t;
                        	Pout << " gamma_t = " << gamma_t << endl;
				Pout << " During collision relative distance = " << delta_n << " between particle " << j << "(type "<<typeJ+1<<")" << " and " << i << "(type "<<typeI+1<<")" << endl;

			}
			          
            		for (int dir=0;dir<3;dir++)
			{	 
				v_t[dir]  = v_r[dir] - v_n[dir];
				F_n[dir]  = k_n * delta_n * n_vec[dir] - gamma_n * v_n[dir];
            		}

			// Liquid tranfer: Cohesion & viscous forces here 
			if (liquid_transfer)
			{
				/*
				if(!first_touch[i][j])
				{
					first_touch[i][j] = true;
					first_touch[j][i] = true;
					calcLiqIn(liq, i, j, liquidVol);
				}
				*/
				if(!first_touch)
				{
					first_touch = true;
					calcLiqIn(liq, i, j, liquidVol);
				}				

				// Capillary force model
				if ( liquidVol[i][j] > 0 )
				{	
					if(verbose) Pout << " Volume of liquid bridge = " << liquidVol[i][j] << endl;
					
					tilde_h = delta_n / radJ;
					tilde_liquidVol = liquidVol[i][j] / ( radJ * radJ * radJ );
					const_A = - 1.1 * pow( tilde_liquidVol, -0.53 );
					const_B = - 0.019 * log( tilde_liquidVol ) + 0.48;
					const_C = 0.0042 * log ( tilde_liquidVol ) + 0.078;
					for (int dir=0;dir<3;dir++)
						fcap[j][dir] = ( exp ( const_A * tilde_h + const_B ) + const_C ) * n_vec[dir] * pi * radJ * surf_tension ;

					// Viscous force model
					Ca = fluid_visc * sqrt(   ( v_r[0] * n_vec[0] ) * ( v_r[0] * n_vec[0] )	
										    + ( v_r[1] * n_vec[1] ) * ( v_r[1] * n_vec[1] )	
										    + ( v_r[2] * n_vec[2] ) * ( v_r[2] * n_vec[2] ) ) / surf_tension;

					ksi = 1 - 1./sqrt( 1 + 2 * tilde_liquidVol / ( pi * tilde_h * tilde_h) );

					for (int dir=0;dir<3;dir++)
						fvisc[j][dir] = ( 3./2. * Ca / tilde_h * ksi * ksi ) * n_vec[dir] * pi * radJ * surf_tension ;

					// Add to normal force
					for (int dir=0;dir<3;dir++)
						F_n[dir] -= fcap[j][dir]; // + fvisc[j][dir];	
				}		
			}

			if(cohesion)
			{
				ri = radI;
				rj = radJ;
				smin = minimumDistanceVdW;    

				for (int dir=0;dir<3;dir++)
				{
					fcoh[j][dir] = cohEnergyDens[typeJ][typeI]/3.0 
							* (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
							* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
							* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir] ;			
					F_n[dir] -= fcoh[j][dir];
				}
			}	 

			/*
            for (int dir=0;dir<3;dir++)
				w_r[dir] = ( radJ - 0.5 * ddist ) * omeJ[dir] / radJ;

            */
            
            		const double crj = radJ - 0.5 * mag(delta_n);
            		const double cri = radI - 0.5 * mag(delta_n);
             
            		for (int dir=0;dir<3;dir++)
                		w_r[dir] = ( crj * omeJ[dir] + cri * omeI[dir] ) / ( ddist + 1.e-25 ) ;
             
			// Tangential relative position
            		v_trel[0] = v_t[0] - ( dist[2]*w_r[1] - dist[1]*w_r[2] );
            		v_trel[1] = v_t[1] - ( dist[0]*w_r[2] - dist[2]*w_r[0] );
            		v_trel[2] = v_t[2] - ( dist[1]*w_r[0] - dist[0]*w_r[1] );
            
			// shear history effects
			if (tangential_history)
			{
				for (int dir=0;dir<3;dir++)
						//delta_t[i][j][dir] += v_trel[dir] * dt;
						//delta_t[dir] += v_trel[dir] * dt;
                        		delta_t[dir] = v_trel[dir] * dt;
                
                // Rotate shear displacement
                		double rsht = delta_t[0] * n_vec[0] + delta_t[1] * n_vec[1] + delta_t[2] * n_vec[2];
                		for (int dir=0;dir<3;dir++)
                        		delta_t[dir] -= rsht * n_vec[dir];

                
				// rotate shear displacements
					//for (int dir=0;dir<3;dir++)
					//delta_t[i][j][dir] -= ( delta_t[i][j][0] * dist[0] + delta_t[i][j][1] * dist[1] + delta_t[i][j][2] * dist[2] ) 
					//				* dist[dir] / radJ / radJ ;
					//delta_t[dir] -= ( delta_t[0] * dist[0] + delta_t[1] * dist[1] + delta_t[2] * dist[2] )
					//				* dist[dir] / radJ / radJ ;

                //}

                //shrmag = sqrt( delta_t[i][j][0] * delta_t[i][j][0] + delta_t[i][j][1] * delta_t[i][j][1] + delta_t[i][j][2] * delta_t[i][j][2]);
                		shrmag = sqrt( delta_t[0] * delta_t[0] + delta_t[1] * delta_t[1] + delta_t[2] * delta_t[2]);

                // tangential forces = shear + tangential position damping
                		for (int dir=0;dir<3;dir++)
                    //F_t[dir] = - k_t * delta_t[i][j][dir]  	
                
                    //Edited by GU
                        		F_t[dir] = - ( k_t * delta_t[dir] );
                        //F_t[dir] = - k_t * delta_t[dir] - gamma_t * v_t[dir];
                    
                // rescale frictional displacements and forces if needed
                		fs = k_t * shrmag; //sqrt( F_t[0] * F_t[0] + F_t[1] * F_t[1] + F_t[2] * F_t[2] );
                		fmu = mu_f * sqrt( F_n[0] * F_n[0] + F_n[1] * F_n[1] + F_n[2] * F_n[2] );

                // energy loss from sliding or damping
                		if (fs > fmu) 
                		{
                    			if (shrmag != 0.0)
                    			{	
                        			for (int dir=0;dir<3;dir++)
                        			{
                                			F_t[dir] *= fmu/fs;
                               				//delta_t[i][j][dir] = - F_t[dir] / k_t;
                             				delta_t[dir] = - F_t[dir] / k_t;
                        			}
                    			}
                    			else F_t[0] = F_t[1] = F_t[2] = 0.0;
               			 }
                		else
                		{
                //Edited by GU	
                    			for (int dir=0;dir<3;dir++)
                        			F_t[dir] -= gamma_t * v_trel[dir];
                //Edited by GU
                		}
            		}else
        		{
                	    fmu = mu_f * sqrt( F_n[0] * F_n[0] + F_n[1] * F_n[1] + F_n[2] * F_n[2] );
                	    const double vtrelmag = sqrt( v_trel[0]*v_trel[0] + v_trel[1]*v_trel[1] + v_trel[2]*v_trel[2] );
                	    const double Ft_damping = gamma_t * vtrelmag;
                	    if( vtrelmag != 0.0 )
                	    {
                		for (int dir=0;dir<3;dir++)
                        	    F_t[dir] = - min( fmu, Ft_damping ) * v_trel[dir] / vtrelmag;
                	    }else
                	    {
                		for (int dir=0;dir<3;dir++)
                        	    F_t[dir] = 0.0;
                	    }
        		}

			for (int dir=0;dir<3;dir++)
			{	
				fcoll[j][dir] += ( F_n[dir] + F_t[dir] );
				 ftan[j][dir] += F_t[dir];

				if(verbose)
				{
					Pout << " F_n[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << F_n[dir];
					Pout << " F_t[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << F_t[dir];
					Pout << " v_j[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << velJ[dir] << " " 
					     << " v_i[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << velI[dir];
					Pout << " x_j[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << posJ[dir] << " " 
					     << " x_i[" << dir << "]= " << setw(IOstream::defaultPrecision() + 6) << posI[dir] << endl;
				}
				//fcoll[i][dir] += ( F_n[dir] + F_t[dir] );
				// ftan[i][dir] += F_t[dir];						 
			}
			
		}
		else if(delta_n < 0)
		{ 
			if(liquid_transfer) //|| ( particles.liquidOn(i) > 0 ) || ( particles.liquidOn(j) > 0 ) )
			{
				if(liquidVol[i][j] > 0 )
				{
					if(verbose) Pout << " During breaking up the bridge relative distance= " << delta_n << " between particle " << i << "\t and " << j << endl;
					for (int dir=0;dir<3;dir++)
					{
						if(verbose) Pout << " position of particle= " << i << " " << posI[dir] << "\t and " << j << " " << posJ[dir] << endl;
					}
					if(verbose) Pout << " During breaking up volume of liquid bridge = " << liquidVol[i][j] << endl;

					double dist_rup = pow(liquidVol[i][j], 1./3.);
					if( delta_n < dist_rup )
					{
						// Capillary force model
						tilde_h = delta_n / radJ;
						tilde_liquidVol = liquidVol[i][j] / ( radJ * radJ * radJ );
						const_A = - 1.1 * pow( tilde_liquidVol, -0.53 );
						const_B = - 0.019 * log( tilde_liquidVol ) + 0.48;
						const_C = 0.0042 * log ( tilde_liquidVol ) + 0.078;
						for (int dir=0;dir<3;dir++)
							fcap[j][dir] = ( exp ( const_A * tilde_h + const_B ) + const_C ) * n_vec[dir] * pi * radJ * surf_tension ;

						// Viscous force model
						Ca = fluid_visc * sqrt(   ( v_r[0] * n_vec[0] ) * ( v_r[0] * n_vec[0] )	
											    + ( v_r[1] * n_vec[1] ) * ( v_r[1] * n_vec[1] )	
											    + ( v_r[2] * n_vec[2] ) * ( v_r[2] * n_vec[2] ) ) / surf_tension;

						ksi = 1 - 1./sqrt( 1 + 2 * tilde_liquidVol / ( pi * tilde_h * tilde_h) );

						for (int dir=0;dir<3;dir++)
						{	
							fvisc[j][dir] = ( 3./2. * Ca / tilde_h * ksi * ksi ) * n_vec[dir] * pi * radJ * surf_tension ;
							fvisc[i][dir] = fvisc[j][dir] ; // For post-processing
						}
						// Add to normal force

						for (int dir=0;dir<3;dir++)
						{
							fcoll[j][dir] -= ( - fcap[j][dir] - fvisc[j][dir] );						
							//fcoll[i][dir] += ( - fcap[j][dir] - fvisc[j][dir] );
						}

					}
					else
					{
						if(verbose) Pout << " Test delta_n > dist_rup " << endl;
						calcLiqBack(liq, i, j, liquidVol);
						//first_touch[i][j] = false;
						//first_touch[j][i] = false;
						first_touch = true;	
					}
				}						
			}

			else if(cohesion)
			{
				ri = radI;
				rj = radJ;
				smin = minimumDistanceVdW;  
				s = - delta_n; // separating distance between the surfaces of the two interacting particles
				smax = (ri + rj) / 4.0; // To speed up the simulation, a maxmimum cutoff separation equal to d/4 is introduced. Beyond smax, the van der Waals cohesive force is not considered.	


				if (s < smax)
				{		
					if (s < smin)
					{	
        					for (int dir=0;dir<3;dir++)
						{
							fcoh[j][dir] =  cohEnergyDens[typeJ][typeI]/3.0 
									* (2.0*ri*rj*(ri+rj+smin))/(smin*smin*(2*ri + 2*rj + smin)*(2*ri + 2*rj + smin)) 
									* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1)   
									* ((smin * (2*ri + 2*rj + smin)) / ((ri + rj + smin)*(ri + rj + smin) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir];
						      fcoll[j][dir] -= fcoh[j][dir];		
						}			
        				} 
        				else
					{ 	
        					for (int dir=0;dir<3;dir++)
						{
							fcoh[j][dir] =  cohEnergyDens[typeJ][typeI]/3.0 
									* (2.0*ri*rj*(ri+rj+s))/(s*s*(2*ri + 2*rj + s)*(2*ri + 2*rj + s)) 
									* ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1)
									* ((s * (2*ri + 2*rj + s)) / ((ri + rj + s)*(ri + rj + s) - (ri-rj)*(ri-rj)) - 1) * n_vec[dir];
						      fcoll[j][dir] -= fcoh[j][dir];		
						}			

					}

				}		

			}
		
		}
		
		// Particle collisional stresses
		/*
		sigma_coll.xx() += ( posI[0] - posJ[0] ) * ( F_n[0] + F_t[0] ); // fcoll[j][0]; 													
		sigma_coll.xy() += ( posI[0] - posJ[0] ) * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll.xz() += ( posI[0] - posJ[0] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2]; 													
		sigma_coll.yy() += ( posI[1] - posJ[1] ) * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll.yz() += ( posI[1] - posJ[1] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		sigma_coll.zz() += ( posI[2] - posJ[2] ) * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		*/
				
		sigma_coll_JI.xx() = 0.5 * dist[0] * ( F_n[0] + F_t[0] ); // fcoll[j][0]; 																										
		sigma_coll_JI.yy() = 0.5 * dist[1] * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll_JI.zz() = 0.5 * dist[2] * ( F_n[2] + F_t[2] ); // fcoll[j][2];				
		sigma_coll_JI.xy() = 0.5 * dist[0] * ( F_n[1] + F_t[1] ); // fcoll[j][1]; 													
		sigma_coll_JI.xz() = 0.5 * dist[0] * ( F_n[2] + F_t[2] ); // fcoll[j][2]; 
		sigma_coll_JI.yz() = 0.5 * dist[1] * ( F_n[2] + F_t[2] ); // fcoll[j][2];
		
		if( delta_n >= 0 )
		{
			if(verbose)
			{
				Pout << " " << endl;
				Pout << " sigma_coll_JI.xx() = " << sigma_coll_JI.xx() << endl; 																										
				Pout << " sigma_coll_JI.yy() = " << sigma_coll_JI.yy() << endl;													
				Pout << " sigma_coll_JI.zz() = " << sigma_coll_JI.zz() << endl;				
				Pout << " sigma_coll_JI.xy() = " << sigma_coll_JI.xy() << endl; 													
				Pout << " sigma_coll_JI.xz() = " << sigma_coll_JI.xz() << endl;
				Pout << " sigma_coll_JI.yz() = " << sigma_coll_JI.yz() << endl;	
			}	
		}
		
	
	}
	
}
