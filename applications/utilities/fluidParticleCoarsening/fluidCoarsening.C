#include "fluidCoarsening.H"
#include "dragForce.H"

namespace Foam
{

	void filteringEulerianScalar
	(
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volScalarField& phi_,
		volScalarField& fphi_ 
	)
	{

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction = 0;
			fphi_[cellI] = 0; 			
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 

			forAll(fcellI,filtercellI)
			{
				barvoidfraction +=  voidfraction_[fcellI[filtercellI]] 
					               * mesh.V()[fcellI[filtercellI]]; 				
				fphi_[cellI]    +=  voidfraction_[fcellI[filtercellI]] 
					               *     phi_[fcellI[filtercellI]] 
						       * mesh.V()[fcellI[filtercellI]]; 
			}
			if ( barvoidfraction > 0 )
			{
				fphi_[cellI] /= barvoidfraction;
			}
			else 
			{
				fphi_[cellI] = 0;
			}
		} 

	}	

	void filteringEulerianVector
	(
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volVectorField& phi_,
		volVectorField& fphi_ 
	)
	{

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction = 0;
			fphi_[cellI] = vector(0,0,0); 			
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 

			forAll(fcellI,filtercellI)
			{
				barvoidfraction +=  voidfraction_[fcellI[filtercellI]] 
					               * mesh.V()[fcellI[filtercellI]]; 				
				fphi_[cellI]    +=  voidfraction_[fcellI[filtercellI]] 
					               *     phi_[fcellI[filtercellI]] 
						       * mesh.V()[fcellI[filtercellI]]; 
			}
			if ( barvoidfraction > 0 )
			{
				fphi_[cellI] /= barvoidfraction;
			}
			else 
			{
				fphi_[cellI] = vector(0,0,0);
			}
		} 

	}
	    
	void filteringEulerianVariables
	(
		const argList& args, 
		const Time& runTime, 
		const fvMesh& mesh,
		const labelListListList& stencillistfilter,
		const int& filterwidth,
		const volScalarField& voidfraction_,
		const volVectorField& U_,
		const volScalarField& p_,
		volScalarField& baralpf_,
		volVectorField& tildeUf_,
		volScalarField& barPg_,
		const volVectorField& Up_,
		volVectorField& tildeUp_,
		const volVectorField& Us_,
		volVectorField& tildeUs_,
		const bool EulerEulerFiltering_,
		const bool EulerianVslipBin_  
	)
	{

		char charfilterwidth[100]; 
		int Filter = 2*(filterwidth+1)+1;
		
		// Filtering volume fraction of gas phase		
		Info<< " 	Filtering voidfraction" << endl;
		sprintf(charfilterwidth, "barvoidfraction_%dX%dX%d",Filter,Filter,Filter);	

		volScalarField barvoidfraction
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
		    voidfraction_	
		);
				
		
		forAll(mesh.cells(),cellI)
		{
		        scalar total_volume = 0;
			barvoidfraction[cellI] = 0; 
			const labelList fcellI = stencillistfilter[filterwidth][cellI]; 
											
			forAll(fcellI,filtercellI)
			{
				total_volume           +=       mesh.V()[fcellI[filtercellI]];
				barvoidfraction[cellI] +=  voidfraction_[fcellI[filtercellI]] 
					                          * mesh.V()[fcellI[filtercellI]]; 
			}
		        if( total_volume > 0 )
			{
				barvoidfraction[cellI] /= total_volume; 
			}
			else 
			{
				barvoidfraction[cellI] = 0;
			}
		} 
		
		Info<< " 	Writing filtered voidfraction" << endl;
	    	barvoidfraction.write();		
		baralpf_ = barvoidfraction;
		
		// Filtering gas velocity		
		Info<< " 	Filtering U" << endl;
		sprintf(charfilterwidth, "tildeU_%dX%dX%d",Filter,Filter,Filter);	
		
		volVectorField tildeU
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
		    U_
		);
		
		filteringEulerianVector
		(
			mesh,
			stencillistfilter,
			filterwidth,
			voidfraction_,
			U_,
			tildeU 
		);
			
		Info<< " 	Writing filtered U" << endl;
	    	tildeU.write();		
		tildeUf_ = tildeU;
		
		// Filtering volume fraction of gas phase		
		Info<< " 	Filtering gas pressure" << endl;
		sprintf(charfilterwidth, "barPg_%dX%dX%d",Filter,Filter,Filter);	
		
		volScalarField barp
		(
		    IOobject
		    (
		        charfilterwidth,
		        runTime.timeName(),
		        mesh,
		        IOobject::NO_READ
		    ),
		    p_	
		);
		
		filteringEulerianScalar
		(
			mesh,
			stencillistfilter,
			filterwidth,
			voidfraction_,
			p_,
			barp 
		);
				
		Info<< " 	Writing filtered gas pressure" << endl;
	    	barp.write();		
		barPg_ = barp;	
				
		if( EulerEulerFiltering_ || EulerianVslipBin_ )
		{			
			// Filtering particle velocity		
			Info<< " 	Filtering Up" << endl;
			sprintf(charfilterwidth, "tildeUp_%dX%dX%d",Filter,Filter,Filter);	

			volVectorField tildeUp
			(
			    IOobject
			    (
		        	charfilterwidth,
		        	runTime.timeName(),
		        	mesh,
		        	IOobject::NO_READ
			    ),
			    Up_
			);
			
			filteringEulerianVector
			(
				mesh,
				stencillistfilter,
				filterwidth,
				1-voidfraction_,
				Up_,
				tildeUp 
			);

			Info<< " 	Writing filtered Up" << endl;
	    		tildeUp.write();		
			tildeUp_ = tildeUp;		

			// Filtering cfdem eulerian velocity		
			Info<< " 	Filtering Us" << endl;
			sprintf(charfilterwidth, "tildeUs_%dX%dX%d",Filter,Filter,Filter);	

			volVectorField tildeUs
			(
			    IOobject
			    (
		        	charfilterwidth,
		        	runTime.timeName(),
		        	mesh,
		        	IOobject::NO_READ
			    ),
			    Us_
			);
			
			filteringEulerianVector
			(
				mesh,
				stencillistfilter,
				filterwidth,
				1-voidfraction_,
				Us_,
				tildeUs 
			);

			Info<< " 	Writing filtered Us" << endl;
	    		tildeUs.write();		
			tildeUs_ = tildeUs;

		}

    }
    
	void filteredEulerEulerDragCoefficient
	(
		cfdemCloud& sm,
		const fvMesh& mesh,
		const volScalarField& p_,	
		const volScalarField& barPg_,
		const volVectorField& Uf_,
		const volVectorField& Up_,		
		const volVectorField& tildeUf_,
		const volVectorField& tildeUp_,		
		const volScalarField& alpf_,
		const volScalarField& baralpf_,		
		const volScalarField& rho_,
		const volVectorField& MappedForce_,
		vectorField& pargii_,
		const bool& verbose_
	)

	{

		// get viscosity field
		#ifdef comp
		    const volScalarField& nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif

		// Parcel/Particle averaged force 
		vectorField parForce(mesh.cells().size(),vector(0,0,0));
				
		// Parcel/Particle averaged drag force
		vectorField parDragForce(mesh.cells().size(),vector(0,0,0));
		
		// Parcel/Particle averaged filtered relative velocity 
		vectorField parFilteredRelativeVelocity(mesh.cells().size(),vector(0,0,0));

		// Parcel/Particle averaged relative velocity
		vectorField parRelativeVelocity(mesh.cells().size(),vector(0,0,0));
						
		volVectorField gradPg = fvc::grad(p_);
                volVectorField gradbarPg = fvc::grad(barPg_);
		
		vector gradPgII(0,0,0);
		vector gradbarPgII(0,0,0);		

                vector Ufluid(0,0,0);
                vector Up(0,0,0);
                vector Ur(0,0,0);
		vector tildeUr(0,0,0);
		
		scalar voidfraction(0);
		scalar barvoidfraction(0);

		scalar ds(0);
		scalar Volp(0);
		
		scalar nuf(0);
		scalar rhof(0);
		
		vector WenYuDrag(0,0,0);
		vector MappedForceCell(0,0,0);
		
		//scalar parForce(0);
		//scalar Beta_fil(0);		
		
		// Correction in vector array
		forAll(mesh.cells(),ii)
		{

                        Ufluid = Uf_[ii]; 
                        Up = Up_[ii];
                        Ur = Ufluid-Up;

			if( mag(Up) > SMALL )
			{

				gradPgII = gradPg[ii]; 
                        	gradbarPgII = gradbarPg[ii];				
	
				tildeUr = tildeUf_[ii] - tildeUp_[ii]; 

                        	ds = 2*sm.radius(0); // Be careful for polydisperse case
				Volp = ds*ds*ds*M_PI/6;

                        	voidfraction = alpf_[ii];
				barvoidfraction = baralpf_[ii];					
				//nuf = nufField[ii];
                        	//rhof = rho_[ii];

                        	// Calculate WenYu Drag 
				//WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);
				//if( ( 1.-alpf_[ii] ) > SMALL ) // If there is no particle in the cell then no drag & correction ...   
				
				MappedForceCell = MappedForce_[ii];
				
				for(int j=0;j<3;j++)
                        	{

					/*
					parForce[ii][j] =   (1.-voidfraction)*gradPgII[j]
			      			   	  - (1.-barvoidfraction)*gradbarPgII[j]
			      			   	  + (1.-voidfraction)*MappedForceCell[j]/Volp;	//+ (1.-voidfraction)*WenYuDrag[j]/Volp; 
					*/
					/*
					parForce[ii][j] =   (1.-voidfraction)*gradPgII[j]
			      			   	  - (1.-barvoidfraction)*gradbarPgII[j]
			      			   	  + MappedForceCell[j]/Volp;					
					*/
					parForce[ii][j] =                  
			      			   	  - (1.-barvoidfraction)*gradbarPgII[j] 
			      			   	  +                  MappedForceCell[j] / mesh.V()[ii];					
							
                		    	parFilteredRelativeVelocity[ii][j] = tildeUr[j] ;

                		    	parDragForce[ii][j] = MappedForceCell[j]/ mesh.V()[ii];		

                		    	parRelativeVelocity[ii][j] = Ur[j];				
				}
				
				if(verbose_  && ii==0 )
			        //if(verbose_)
				{
				    Pout << "CellID = " << ii << endl;
				    Pout << " parForceMain				= " << parForce[ii]             		<< endl;
			            Pout << " parFilteredRelativeVelocityMain	= " << parFilteredRelativeVelocity[ii]  	<< endl;   
                		    Pout << " parDragForceMain			= " << parDragForce[ii]                 	<< endl; 
				    Pout << " parRelativeVelocityMain		= " << parRelativeVelocity[ii]          	<< endl;
			    	    Pout << " " << endl;
				}				

			}
		}

		// Parcel/Particle averaged filtered drag coefficient
		vectorField parbarBetai_(mesh.cells().size(),vector(0,0,0));
		
		// Drag correction with particle coarsening
		pargii_.resize(mesh.cells().size());
		
		forAll(mesh.cells(),ii)
		{ 	
			for(int j=0;j<3;j++) 
			{				
			 parbarBetai_[ii][j] =        parForce[ii][j] / ( parFilteredRelativeVelocity[ii][j] + SMALL ) ;
			      pargii_[ii][j] =    parbarBetai_[ii][j] *           parRelativeVelocity[ii][j] 
				                       / ( parDragForce[ii][j] + SMALL ) ;				
			}
			
                        Up = Up_[ii];
			//if(verbose_ && ii==0)
			if(verbose_)
			{
				if( mag(Up) > SMALL && ii==0 )
				{
					//int index = 0;			
					Pout << " Eulerian barBetai = " << parbarBetai_[ii]  << endl;
					Pout << " Eulerian gii      = " << pargii_[ii]       << endl;
					Pout << " " << endl;		
				}
			}
		}			
													 
   	}    
}	
