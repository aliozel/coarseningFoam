#include "coarseningPar.H"
#include "dragForce.H"
 
namespace Foam
{ 

 	void coarseningPar
	(
		cfdemCloud& parMainsm,
		cfdemCloud& parsm,
		const volScalarField& rho_,		
		const volScalarField& alpf_,
		const volScalarField& baralpf_,
		const volVectorField& Uf_,
		const volVectorField& tildeUf_,
		const volScalarField& p_,
		const volScalarField& barPg_,
		vectorField& pargii_,
		const labelListList& parMainList_,
		const labelListList& parList_,
		const bool& verbose_,
        	const word& par_,
		const bool& velSlipParcel_,
		const volVectorField& Up_,
		const bool& velSlipParcelInt_  		
	)
 	{
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = parMainsm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = parMainsm.turbulence().nu();
		#endif
		
		// Distance vector
		vector dist(0,0,0);
				
		// Parcel/Particle averaged force 
		vectorField parForce(parMainList_.size(),vector(0,0,0));
				
		// Parcel/Particle averaged drag force
		vectorField parDragForce(parMainList_.size(),vector(0,0,0));
		
		// Parcel/Particle averaged filtered relative velocity 
		vectorField parFilteredRelativeVelocity(parMainList_.size(),vector(0,0,0));

		// Parcel/Particle averaged relative velocity
		vectorField parRelativeVelocity(parMainList_.size(),vector(0,0,0));	
			
		volVectorField gradp_ = fvc::grad(p_);
		volVectorField gradbarPg_ = fvc::grad(barPg_);
			   
		interpolationCellPoint<vector> gradPgInterpolator_(gradp_);
		interpolationCellPoint<vector> gradbarPgInterpolator_(gradbarPg_);
		
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		interpolationCellPoint<vector> tildeUInterpolator_(tildeUf_);
		
		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		interpolationCellPoint<scalar> barvoidfractionInterpolator_(baralpf_);
				
		// Calculate slip velocity by interpolation from the slip velocity at the cell center (Parcel approach)
		volVectorField velSlip = Uf_ - Up_;
		volVectorField tildeVelSlip = tildeUf_ - Up_;
		
		interpolationCellPoint<vector> velSlipInterpolator_(velSlip);
		interpolationCellPoint<vector> tildeVelSlipInterpolator_(tildeVelSlip);
		
		// Interpolate particle velocity at cell center
		interpolationCellPoint<vector> UpMainIntInterpolator_(Up_);			

        	vector gradPgMain(0,0,0);
        	vector gradbarPgMain(0,0,0);
        	vector positionMain(0,0,0);
        	vector UpMain(0,0,0);
        	vector UfluidMain(0,0,0);
        	vector tildeUfluidMain(0,0,0);
        	vector UrMain(0,0,0);
        	vector tildeUrMain(0,0,0);
        
		vector gradPg(0,0,0);	
	   	vector gradbarPg(0,0,0);		
		vector position(0,0,0);
		vector Up(0,0,0);
		vector Ufluid(0,0,0);
		vector tildeUfluid(0,0,0);
		vector Ur(0,0,0);
		vector tildeUr(0,0,0);

       		scalar VolpMain(0);
        	scalar dsMain(0);
        
		scalar Volp(0);
		scalar ds(0);
		
		vector WenYuDragMain(0,0,0);
		scalar voidfractionMain(0);
        
        	vector WenYuDrag(0,0,0);
        	scalar voidfraction(0);

        	scalar rhofMain(0);
        	scalar nufMain(0);
        
		scalar rhof(0);
		scalar nuf(0);
		
        	label cellIMain(-1);
		label cellI(-1);
		
		vector WenYuDragMainTildeMicro(0,0,0);
		vector WenYuDragTildeMicro(0,0,0);
		
		scalar barvoidfraction(0);
		scalar barvoidfractionMain(0);
		
		//Pout << " " << par_ << " coarsening working " << endl;
		
		for(int parI=0; parI < parMainList_.size(); parI++ )
		{
            
        	    int parIID = parMainList_[parI][0];

        	    // Main particle/parcel variables
        	    cellIMain = parMainsm.cellIDs()[parIID][0];
		    			
        	    //if( cellIMain > -1 )
               	    //{			    
			    
			    positionMain = parMainsm.position(parIID);

        		       gradPgMain =    gradPgInterpolator_.interpolate(positionMain,cellIMain);
        		    gradbarPgMain = gradbarPgInterpolator_.interpolate(positionMain,cellIMain);

        		         UfluidMain =      UInterpolator_.interpolate(positionMain,cellIMain);
        		    tildeUfluidMain = tildeUInterpolator_.interpolate(positionMain,cellIMain); //tildeUf_[cellIMain]

        		    UpMain = parMainsm.velocity(parIID);
			    
			    if( velSlipParcelInt_ && par_ == "parcel" ) UpMain = UpMainIntInterpolator_.interpolate(positionMain,cellIMain);
			    
        		         UrMain =      UfluidMain - UpMain;
        		    tildeUrMain = tildeUfluidMain - UpMain;

			    // Calculate slip velocity by interpolation from the slip velocity at the cell center (Parcel approach)
			    if( velSlipParcel_ && par_ == "parcel" ) tildeUrMain = tildeVelSlipInterpolator_.interpolate(positionMain,cellIMain);
			    //if( velSlipParcel_ ) tildeUrMain = tildeVelSlipInterpolator_.interpolate(positionMain,cellIMain);		
				
        		    dsMain = 2*parMainsm.radius(parIID);			    			    
        		    VolpMain = dsMain*dsMain*dsMain*M_PI/6;

        		    // Calculate WenYu Drag
        		       voidfractionMain =    voidfractionInterpolator_.interpolate(positionMain,cellIMain);
        		    barvoidfractionMain = barvoidfractionInterpolator_.interpolate(positionMain,cellIMain);
			    
			    nufMain = nufField[cellIMain];
        		    rhofMain = rho_[cellIMain];
			    
        		    WenYuDragForce(     UrMain,dsMain,rhofMain,nufMain,   voidfractionMain,WenYuDragMain);
			    WenYuDragForce(tildeUrMain,dsMain,rhofMain,nufMain,barvoidfractionMain,WenYuDragMainTildeMicro);
			    
			    if(verbose_ && parI==0)
			    //if(verbose_)
			    {
			    	    Pout << " " << endl;
				    Pout << " " << endl;
				    Pout << " " << endl;	
				    if(par_ == "particle") Pout << " Particle 				= " << parIID            			<< endl; 
				    if(par_ == "parcel"  ) Pout << " Parcel 				= " << parIID            			<< endl; 
			            Pout << " cellIMain				= " << cellIMain         			<< endl; 
				    Pout << " positionMain@Xp			= " << positionMain      			<< endl;
				    Pout << " gradPg_intMain@Xp			= " << gradPgMain    			<< endl;
				    Pout << " gradbarPg_intMain			= " << gradbarPgMain 			<< endl;
				    Pout << " UfluidMain@Xp				= " << UfluidMain        			<< endl; 
        		            Pout << " tildeUfluidMain@Xp			= " << tildeUfluidMain   			<< endl;
        		            Pout << " UpMain					= " << UpMain 					<< endl;
        		            Pout << " voidfractionMain@Xp			= " << voidfractionMain           		<< endl;
				    Pout << " voidfraction@Cell			= " << alpf_[cellIMain]           		<< endl;	
        		            Pout << " WenYuDragMain 				= " << WenYuDragMain              		<< endl; 	
			    }

        		    if ( par_ == "particle" )
        		    {
				for(int j=0;j<3;j++)
                		{

                		    parForce[parI][j] = (   VolpMain * gradbarPgMain[j]
                                        		  - VolpMain *    gradPgMain[j]
                                        		  + WenYuDragMain[j]                )  ;

                		    parFilteredRelativeVelocity[parI][j] = tildeUrMain[j];

                		    //parDragForce[parI][j] = WenYuDragMain[j];
				    parDragForce[parI][j] = WenYuDragMainTildeMicro[j];
				    	
                		    //parRelativeVelocity[parI][j] = UrMain[j];
				    
                		}
			        if(verbose_ && parI==0)
				//if(verbose_)
			        {
				    Pout << " parForceMain				= " << parForce[parI]             		<< endl;
			            Pout << " parFilteredRelativeVelocityMain	= " << parFilteredRelativeVelocity[parI]  	<< endl;   
                		    Pout << " parDragForceMain			= " << parDragForce[parI]                 	<< endl; 
				    Pout << " parRelativeVelocityMain		= " << parRelativeVelocity[parI]          	<< endl;
			    	}			
        		    }
        		    else if ( par_ == "parcel" )
        		    {
				// Sub-particles list
                		labelList parList_L = parList_[parI];
				
				//Pout << " Parcel coarsening test " << endl;

                		forAll(parList_L,particleI)
                		{				
                		    label particleII = parList_L[particleI]; 

                		    cellI = parsm.cellIDs()[particleII][0];

                		    if( cellI > -1 )
                		    {
                        		position = parsm.position(particleII); 

                        		   gradPg =    gradPgInterpolator_.interpolate(position,cellI); 
                        		gradbarPg = gradbarPgInterpolator_.interpolate(position,cellI); 

                        		     Ufluid =      UInterpolator_.interpolate(position,cellI); 
                        		tildeUfluid = tildeUInterpolator_.interpolate(position,cellI);	

                        		Up = parsm.velocity(particleII);
                  
		        	 	     Ur = Ufluid-Up;
                        		tildeUr = tildeUfluid-Up;

                        		ds = 2*parsm.radius(particleII);
                        		Volp = ds*ds*ds*M_PI/6;

                        		// Calculate WenYu Drag 
                        		   voidfraction =    voidfractionInterpolator_.interpolate(position,cellI);
                        		barvoidfraction = barvoidfractionInterpolator_.interpolate(position,cellI);
					
					nuf = nufField[cellI];
                        		rhof = rho_[cellI];
						
                        		WenYuDragForce(     Ur,ds,rhof,nuf,   voidfraction,WenYuDrag);										
					WenYuDragForce(tildeUr,ds,rhof,nuf,barvoidfraction,WenYuDragTildeMicro);	
					
					/*
					Pout << " parcel ID = " << parIID << " particle ID = " << particleII << " parI = " << parI << " of " 
					     << parList_.size() << " dsMain = " << dsMain << " ds = " << ds 
					     << " cellIMain = " << cellIMain << " cellI = " << cellI << endl;
					*/
				 
                        		for(int j=0;j<3;j++)
                        		{

					    parForce[parI][j] += (   VolpMain * gradbarPgMain[j]
                                                		   - Volp     *    gradPg[j]
                                                		   + WenYuDrag[j]	         	)  ;	
								   								   				
                        		    parFilteredRelativeVelocity[parI][j] += tildeUrMain[j];

                        		    //parDragForce[parI][j] += WenYuDrag[j];
					    parDragForce[parI][j] += WenYuDragTildeMicro[j];

                        		    // parRelativeVelocity[parI][j] += Ur[j] ; 						

                        		}
					
					if(verbose_&&parI==0)
					//if(verbose_)
					{
					    Pout << " " << endl;		
					    Pout << " particle in parcel               	= " << particleII            			<< endl; 
			        	    Pout << " cellI                       		= " << cellI         		        	<< endl; 
					    Pout << " position                    		= " << position      				<< endl;
					    Pout << " gradPg_int@Xp				= " << gradPg    				<< endl;
					    Pout << " gradbarPg_int@Xp               	= " << gradbarPg 				<< endl;
					    Pout << " Ufluid@Xp                   		= " << Ufluid        				<< endl; 
        		        	    Pout << " tildeUfluid@Xp              		= " << tildeUfluid   				<< endl;
        		        	    Pout << " Up                          		= " << Up	 				<< endl;
        		        	    Pout << " voidfraction@Xp             		= " << voidfraction           			<< endl;
					    Pout << " voidfraction@Cell           		= " << alpf_[cellI]           			<< endl;	
        		        	    Pout << " WenYuDrag                   		= " << WenYuDrag              			<< endl; 			
					    Pout << " parForce                    		= " << parForce[parI]         			<< endl;
					    Pout << " parFilteredRelativeVelocity 		= " << parFilteredRelativeVelocity[parI]  	<< endl;   
                			    Pout << " parDragForce                		= " << parDragForce[parI]                 	<< endl; 
					    Pout << " parRelativeVelocity         		= " << parRelativeVelocity[parI]          	<< endl;
					}					

                		  }

                		}

        		    }
			   
		    //}		    

		}				

		// Parcel/Particle averaged filtered drag coefficient
		vectorField parbarBetai_(parMainList_.size(),vector(0,0,0));
		
		// Drag correction with particle coarsening
		pargii_.resize(parMainList_.size());
		
		//forAll(parList_,parII)
		for(int ii=0; ii < parMainList_.size(); ii++ )
		{ 	
			for(int j=0;j<3;j++) 
			{				
			 //parbarBetai_[ii][j] =        parForce[ii][j] / ( parFilteredRelativeVelocity[ii][j] + SMALL ) ;
			 //     pargii_[ii][j] =    parbarBetai_[ii][j] *           parRelativeVelocity[ii][j] 
			 //	                       / ( parDragForce[ii][j] + SMALL ) ;				

			 parbarBetai_[ii][j] = parForce[ii][j] / ( parFilteredRelativeVelocity[ii][j] + SMALL ) ;
			      pargii_[ii][j] = parForce[ii][j] / ( parDragForce[ii][j] + SMALL ) ;				 

			}
			
			if(verbose_ && ii==0)
			//if(verbose_)
			{
				//int index = 0;			
				Pout << " " << endl;
				if(par_ == "particle") Pout << " Particle ";
				if(par_ == "parcel"  ) Pout << " Parcel "  ;
				Pout << "barBetai = " << parbarBetai_[ii]  << endl;
				if(par_ == "particle") Pout << " Particle ";
				if(par_ == "parcel"  ) Pout << " Parcel "  ;
				Pout << "gii      = " << pargii_[ii]       << endl;		
			}

		}		
			
	}
	
}	
