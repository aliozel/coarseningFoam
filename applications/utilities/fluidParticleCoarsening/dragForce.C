#include "dragForce.H"

namespace Foam
{	
	void WenYuDragForce
	(		
		const vector& Ur_,
		const scalar& ds_,
		const scalar& rho_,
		const scalar& nuf_,
		const scalar& voidfraction_,
		vector& WenYuDrag_
	)
	{
		scalar Rep(0);
		scalar magUr(0);
		scalar CD(0);
		scalar betaP(0);
		
		scalar Volp(0);	
		
		magUr = mag(Ur_);		
		Volp = ds_*ds_*ds_*M_PI/6;

        	if(ds_ <= 0) Pout << " ds " << ds_ << endl;
		if(voidfraction_> 1.0 ) Pout << " Int. voidfraction > 1 " << " value= " << voidfraction_ << endl;				
		
		if (magUr > 0 && ds_ > 0)
		{
			// calc particle Re number
			Rep = voidfraction_*ds_*magUr/(nuf_+SMALL);
            
			// calc CD
			if (Rep < 1000)
			{
				CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction_,-1.65);
			}
			else
			{
				CD = 0.44*pow(voidfraction_,-1.65); 
			}  	

			// calc drag coefficient 
			betaP = 3./4.*rho_*CD*magUr/ds_;  	

			// calc particle's drag
			WenYuDrag_ = Volp*betaP*Ur_;

		}else
		{
			WenYuDrag_ = vector(0,0,0);
		}	
						
	}

	void CalculateDragForce
	(
		cfdemCloud& sm,
		const volScalarField& alpf_,
		const volVectorField& Uf_,
		const volScalarField& rho_,
		const bool& verbose_,
		vectorField& DragForce_,
		const labelListList& particleList_
	)
	{
		
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif

		// Local variables	
		label  cellI(-1);
		vector drag(0,0,0);
		vector Ufluid(0,0,0);
		
		vector position(0,0,0);
		scalar voidfraction(1);
		
		vector Up(0,0,0);
		vector Ur(0,0,0);
		scalar ds(0);
		
		scalar nuf(0);
		scalar rhof(0);
		
		vector WenYuDrag(0,0,0);
		
		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
				
		// 
		//_AO_Parallel
		DragForce_.resize(particleList_.size());
					
	    	for(int ii =0; ii < particleList_.size(); ii++)
		{
			int index = particleList_[ii][0];
			cellI = sm.cellIDs()[index][0];
			drag = vector(0,0,0);
			Ufluid = vector(0,0,0);
			WenYuDrag = vector(0,0,0);
			DragForce_[ii] = vector(0,0,0);
			    
			if (cellI > -1) // particle Found
			{
				position = sm.position(index);
			
				if ( alpf_[cellI] > 1. ) Pout << " voidfraction > 1 " << alpf_[cellI] << endl;
			
				voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
				Ufluid = UInterpolator_.interpolate(position,cellI);
			
				if ( voidfraction > 1. ) 
				{
						Pout << " Int. voidfraction > 1 " << " value= " << voidfraction;
						voidfraction = alpf_[cellI];
						Pout << " mod. value = " << voidfraction << endl;
				}		
				Up = sm.velocity(index);
			
				Ur = Ufluid-Up;				
				ds = 2*sm.radius(index);		
				rhof = rho_[cellI];
				nuf = nufField[cellI];
			
				// Drag force
				WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);
						
				if(verbose_ && index <= 1)
				{
					Info << "" << endl;
					Pout << " index = " << index << endl;
					Pout << " position = " << position << endl; 
					Pout << " Up = " << Up << endl;
					Pout << " Ur = " << Ur << endl;
					Pout << " dp = " << ds << endl;
					Pout << " rho = " << rhof << endl;
					Pout << " nuf = " << nuf << endl;
					Pout << " voidfraction = " << voidfraction << endl;
					Pout << " drag = " << WenYuDrag << endl;
					Info << " " << endl;
				}
			}	
			for(int j=0;j<3;j++) DragForce_[ii][j] = WenYuDrag[j];
		}	
	}
	
	void EulerianParticleVelocityForce
	(
		cfdemCloud& sm,			
		const fvMesh& mesh,
		volVectorField& Uf_,
		volVectorField&	Up_,
		volScalarField& rho_,
		volScalarField& alpf_,
		volScalarField& Pg_,
		volVectorField& MappedDragForce_,
		const labelListList& particleList_,
		const bool& weighting_
	)
	{		
		// Neighbouring cells
		CPCCellToCellStencil neighbourCells(mesh);
				
		// get viscosity field
		#ifdef comp
		    const volScalarField nufField = sm.turbulence().mu()/rho_;
		#else
		    const volScalarField& nufField = sm.turbulence().nu();
		#endif

		// Gas pressure gradient
		volVectorField gradPg_ = fvc::grad(Pg_);
		interpolationCellPoint<vector> gradPgInterpolator_(gradPg_);

		// Local variables	
		label  cellID(-1);
		vector drag(0,0,0);
		vector Ufluid(0,0,0);
		
		vector position(0,0,0);
		scalar voidfraction(1);
		
		vector Up(0,0,0);
		vector Ur(0,0,0);
		scalar ds(0);
		
		scalar nuf(0);
		scalar rhof(0);
		
		vector WenYuDrag(0,0,0);
		
		interpolationCellPoint<scalar> voidfractionInterpolator_(alpf_);
		interpolationCellPoint<vector> UInterpolator_(Uf_);	
		
		scalar dist_s(0);
		scalar sumWeights(0);
		
		scalarField               weightScalar(27,scalar(0.0));
		Field <Field <scalar> >   particleWeights(particleList_.size(),weightScalar);
		
		//Info << " particle size " << particleList_.size() << endl;
		
		// Number of particle in a cell
		scalarField np(mesh.cells().size(),scalar(0));
		
		// Particle volume
		scalar Volp(0);
		vector gradPg_int(0,0,0);
		
		for(int ii = 0; ii < particleList_.size(); ii++)
		{
			int index = particleList_[ii][0];
			
			cellID = sm.cellIDs()[index][0];
			position = sm.position(index);			    

                        Ufluid = UInterpolator_.interpolate(position,cellID); 
			Up = sm.velocity(index);
                        Ur = Ufluid-Up;

                        ds = 2*sm.radius(index);

                        // Calculate WenYu Drag 
                        voidfraction = voidfractionInterpolator_.interpolate(position,cellID);
                        nuf = nufField[cellID];
                        rhof = rho_[cellID];	
                        WenYuDragForce(Ur,ds,rhof,nuf,voidfraction,WenYuDrag);	
    
        		Volp = ds*ds*ds*M_PI/6;
			gradPg_int = gradPgInterpolator_.interpolate(position,cellID);
			
			//if (cellID > -1)  // particle centre is in domain
            		//{
				if(weighting_)
				{
					labelList& cellsNeigh = neighbourCells[cellID];
					sumWeights = 0;
					dist_s = 0;

					//Info << " index = " << index << " ii = " << ii << " cellID = " << cellID << endl;

					forAll(cellsNeigh,jj)
					{
						// Find distances between particle and neighbouring cells					
						dist_s = mag(sm.mesh().C()[cellsNeigh[jj]]-position)/pow(sm.mesh().V()[cellsNeigh[jj]],1./3.);

						if(dist_s <= 0.5)
						{		
							particleWeights[ii][jj] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
						}
						else if (dist_s > 0.5 && dist_s <= 1.5)
						{		
							particleWeights[ii][jj] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
						}
						else if (dist_s > 1.5 && dist_s <= 2.5)
						{		
							particleWeights[ii][jj] =  pow(2.5-dist_s,4)/24.;
						}
						else
						{		
							particleWeights[ii][jj] = 0;
						}

						sumWeights += particleWeights[ii][jj];

					}	

					forAll(cellsNeigh,jj)
					{	
						if ( sumWeights != 0 )
						{
							Up_[cellID] 	         +=  Up*particleWeights[ii][jj]/sumWeights;
							MappedDragForce_[cellID] += (WenYuDrag + Volp * gradPg_int) * particleWeights[ii][jj]/sumWeights;
						}
						else
						{
							Up_[cellID] 		 = vector(0,0,0);
							MappedDragForce_[cellID] = vector(0,0,0);	
						}
					}
				}
				else
				{
					Up_[cellID] 	         += Up;
					MappedDragForce_[cellID] += WenYuDrag + Volp * gradPg_int;						
				}

		}
		
	}

}	
