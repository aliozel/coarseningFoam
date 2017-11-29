#include "createParcel.H"
 
namespace Foam
{ 
 	void createParcels
	(
		const fvMesh& mesh,			
		cfdemCloud& sm,
        	cfdemCloud& parcelsm,
		const int& nparticle_,
		const bool verbose_,
		const scalar& minX_,
		const scalar& maxX_,
		const scalar& minY_,
		const scalar& maxY_,
		const scalar& minZ_,
		const scalar& maxZ_,
	        //ANNpointArray& dataPts,
	        //ANNpoint& queryPt,
	        //ANNkd_tree*& kdTree,
		double **& parcelPositions,
                double **& parcelVelocities,
                double **& parcelRadii,
                double **& parcelCellID,
		const labelListList& particleList_,
        	labelListList& parcelList_,
        	labelListList& parcelParticleList_,
		const bool& useAllParticles_
	)
	{		

		if ( nparticle_ > particleList_.size() )
		{
			FatalError << " Number of particles in a parcel > number of particles" << abort(FatalError);
		}
		
		if ( nparticle_ < 1 )
		{
			FatalError << " Number of particles < 0 in a parcel " << abort(FatalError);
		}
		
		// Number of particles in a parcel
		//int k = nparticle_; 		
		int k = 3 * nparticle_; 		
		//int k = particleList_.size() / nparticle_ / 4;
		if (useAllParticles_) k = nparticle_; 
		
		// Dimensions, exact OR approximate  
		int dim =3; double eps = 0;
						
		// Number of points
		int nPts;
		nPts = particleList_.size();
		
		Pout << " Number of particles in a parcel = " << nparticle_ << endl;		
		
		// Squared radius
		const scalar sqRad = pow( (maxX_ - minX_) * (maxX_ - minX_) 
				     	 +(maxY_ - minY_) * (maxY_ - minY_)
					 +(maxZ_ - minZ_) * (maxZ_ - minZ_), 0.5); // / 10.;
		
		Pout << " Squared radius = " << sqRad << endl;
				
		// Data points
		ANNpointArray dataPts;
		// Query points
		ANNpoint queryPt;
		
		ANNidxArray	nnIdx;              // 	near neighbour indices
		ANNdistArray 	dists;			//	near neighbour distances
		ANNkd_tree* 	kdTree;			//	search structure
	
		Pout << " Created kdTree variables " << endl;
	
		// Allocate 
		queryPt = annAllocPt(dim);
		dataPts = annAllocPts(nPts, dim);
		nnIdx = new ANNidx[k];
		dists = new ANNdist[k];

		Pout << " Allocated kdTree variables " << endl;
        
        	labelList particleCreateParcelList(particleList_.size());
				        
        	for(int ii =0; ii < particleList_.size(); ii++)
		{
				label particleGlobalID = particleList_[ii][0];
				dataPts[ii][0] = sm.position(particleGlobalID).x();
				dataPts[ii][1] = sm.position(particleGlobalID).y();
				dataPts[ii][2] = sm.position(particleGlobalID).z();
                		particleCreateParcelList[ii] = particleList_[ii][0];
		}


        	Pout << " Creating kdTree..." << endl;
        	kdTree = new ANNkd_tree(dataPts, nPts, dim);
    
        	// Parcel center of mass
        	vectorField parcelCenter(particleList_.size(),vector(0,0,0));
        	// Calculate weighted velocities
        	vectorField parcelVel(particleList_.size(),vector(0,0,0));
        	// Parcel diameter
        	scalarField parcelDiameter(particleList_.size(),scalar(0.0));
        
        	// Local parcel variable
        	vector Up(0,0,0);
        	vector position(0,0,0);
       
		Pout << " Entering particle loops to create parcel" << endl; 
 
    		label parcelI = 0;
    		forAll(particleList_,ii)
		{				
	    		label particleGlobalID = particleList_[ii][0];
        
        		if ( particleCreateParcelList[ii] > -1 )
        		{
        		    queryPt[0] = sm.position(particleGlobalID).x();
        		    queryPt[1] = sm.position(particleGlobalID).y();
        		    queryPt[2] = sm.position(particleGlobalID).z();

        		    kdTree->annkFRSearch(
                                		    queryPt,			// query point					
                                		    sqRad,				// squared radius
                                		    k,                  // number of the near neighbours to return
                                		    nnIdx,				// nearest neighbor array
                                		    dists,				// dist to near neighbours
                                		    eps			);

        		    int partSum = 0;
			    int i = 0;
        		    //for (int i = 0; i < k; i++)       		    
			    while( i < k && partSum < nparticle_ )
			    {
                		//Pout << " nnIdx[i] = " << nnIdx[i] << " particleCreateParcelList[nnIdx[i]] = " << particleCreateParcelList[nnIdx[i]] << endl;
                		
				if ( particleCreateParcelList[nnIdx[i]] != -1 )
                		{
                		    if (!useAllParticles_) particleCreateParcelList[nnIdx[i]] = -1 ;

                		    // Parcel velocity, position, diameter
                		    Up = sm.velocity(particleList_[nnIdx[i]][0]);
                		    position = sm.position(particleList_[nnIdx[i]][0]);

                		    for(int j=0;j<3;j++)
                		    {
                        		// Parcel velocity
                        		parcelVel[parcelI][j]    += Up[j];
                        		// Parcel center of mass
                        		parcelCenter[parcelI][j] += position[j];

                		    }
                		    parcelDiameter[parcelI] = 2.*sm.radius(particleList_[nnIdx[i]][0]);
                		    parcelParticleList_[parcelI][partSum] = particleList_[nnIdx[i]][0];
                		    partSum++;
				    if( partSum == nparticle_ ) parcelI++;
                		};
				i++;				
        		    };
			    
        		}
        
    		}
       
		//Pout << " parcelParticleList_.size " << parcelParticleList_.size() << endl;		
		parcelParticleList_.resize(parcelI);
		//Pout << " Particles are located in parcels..." << " local number of parcels " << parcelI << endl;
		//Pout << " parcelParticleList_.size " << parcelParticleList_.size() << endl;
        	// Parallel computation
        	//label totalNpParcel = particleList_.size();
        	int totalNpParcel = parcelI;
        	if(Pstream::parRun())
        	{
        	    reduce(totalNpParcel, sumOp<int>());
        	}

        	Pout << " Total number of parcels = " << totalNpParcel << endl;

        	// Set number of parcels
        	parcelsm.setNumberOfParticles(totalNpParcel);
        	//parcelsm.reAllocArraysPost();
		
		/*				
		List < scalar > listScalar(1);	 
		List < List < scalar > > listScalarScalar(particleList_.size(),listScalar);
		List < List < List < scalar > > > listParcelRadiis(Pstream::nProcs(),listScalarScalar);
		
		List < vector > listVector(1);	 
		List < List < vector > > listVectorVector(particleList_.size(),listVector);
		List < List < List < vector > > > listParcelPositions(Pstream::nProcs(),listVectorVector);					
		List < List < List < vector > > > listParcelVelocities(Pstream::nProcs(),listVectorVector);	
		*/
		
		List < scalar > listScalar(1);	 
		List < List < scalar > > listScalarScalar(parcelI,listScalar);
		List < List < List < scalar > > > listParcelRadiis(Pstream::nProcs(),listScalarScalar);
		
		List < vector > listVector(1);	 
		List < List < vector > > listVectorVector(parcelI,listVector);
		List < List < List < vector > > > listParcelPositions(Pstream::nProcs(),listVectorVector);					
		List < List < List < vector > > > listParcelVelocities(Pstream::nProcs(),listVectorVector);		
			
		//for(int index = 0; index <  npLocalParcel; index++)
		//for(int index = 0; index < particleList_.size(); index++)
        	for(int index = 0; index < parcelParticleList_.size(); index++)
		{
			
			listParcelRadiis[Pstream::myProcNo()][index][0] = parcelDiameter[index]/2.;
            
            		listParcelPositions[Pstream::myProcNo()][index] = vector(parcelCenter[index][0]/parcelParticleList_[index].size(),
                                                                     		 parcelCenter[index][1]/parcelParticleList_[index].size(),
                                                                     		 parcelCenter[index][2]/parcelParticleList_[index].size()   );
                                                                     
 			//listParcelVelocities[Pstream::myProcNo()][index] = vector(parcelVel[index][0],parcelVel[index][1],parcelVel[index][2]);

            		listParcelVelocities[Pstream::myProcNo()][index] = vector(parcelVel[index][0]/parcelParticleList_[index].size(),
                                                                      		  parcelVel[index][1]/parcelParticleList_[index].size(),
                                                                      		  parcelVel[index][2]/parcelParticleList_[index].size());		
																	  
		}
				
		// Distrubute to all processors
		Pstream::gatherList(listParcelRadiis);
		Pstream::scatterList(listParcelRadiis);

		Pstream::gatherList(listParcelPositions);
		Pstream::scatterList(listParcelPositions);

		Pstream::gatherList(listParcelVelocities);
		Pstream::scatterList(listParcelVelocities);

		Pout << " Pstream gather/scatter parcel vel/pos/radii done .. " << endl;
											
		// Create global parcels
		List < List < scalar > > globalListParcelRadiis;
	
		globalListParcelRadiis = ListListOps::combine< List < List < scalar > > >
		(
			listParcelRadiis,
			accessOp< List < List < scalar > > >()
		);			
		
		// Create global parcels
		List < List < vector > > globalListParcelPositions;

		globalListParcelPositions = ListListOps::combine< List < List < vector > > >
		(
			listParcelPositions,
			accessOp< List < List < vector > > >()
		);			
					
		// Create global parcels
		List < List < vector > > globalListParcelVelocities;

		globalListParcelVelocities = ListListOps::combine< List < List < vector > > >
		(
			listParcelVelocities,
			accessOp< List < List < vector > > >()
		);			
		
        	// Define Parcel Cloud positions, velocities, radius and cellIDs

		Pout << " Before global list streaming ..." << endl;
	
		if(Pstream::master)
		{
			for(int index = 0; index <  parcelsm.numberOfParticles(); index++)
			{	
				for(int j=0;j<3;j++) 
				{			
					parcelPositions[index][j] = globalListParcelPositions[index][0][j];
					parcelVelocities[index][j] = globalListParcelVelocities[index][0][j] ;		
				}							
				parcelRadii[index][0] = globalListParcelRadiis[index][0] ;
			}
			
		}

		Pout << " After Pstream master " << endl;
        
		// Send to child CPUs
		if(Pstream::parRun())
		{
			for(int index = 0; index < parcelsm.numberOfParticles(); index++)
			{
				for(int idir = 0; idir < 3; idir++)
				{
					Pstream::scatter(parcelPositions[index][idir]);
					Pstream::scatter(parcelVelocities[index][idir]);
				}
				Pstream::scatter(parcelRadii[index][0]);
			}
		}
	

		Pout << " Parcel pos/vel/radii distributed ..." << endl;
						
       	 	parcelsm.locateM().findCell(NULL,parcelPositions,parcelCellID,parcelsm.numberOfParticles());
        	parcelsm.setPos(parcelPositions);
        	parcelsm.setVel(parcelVelocities);

		Pout << " After setting parcel pos/vel/cellID " << endl;
	
        	int  NparcelCellID = 0;
        	for(int parcelII = 0; parcelII < parcelsm.numberOfParticles(); parcelII++)
        	{
        	    label cellI = parcelsm.cellIDs()[parcelII][0];
		    //Pout << " parcel radii in createParcel.C " << parcelRadii[parcelII][0] << endl;
		    //if( cellI > -1 ) Pout << " Parcel = " << parcelII << " cellI = " << cellI << endl;	
        	    if( cellI > -1 )
        	    {
                	parcelList_[NparcelCellID][0] = parcelII;
                	NparcelCellID++;
			//Pout << " Parcel = " << parcelII << " cellI = " << cellI << " NparcelCellID = " << NparcelCellID << endl;
			//Pout << " Parcel = " << parcelII << " cellI = " << parcelsm.cellIDs()[parcelII][0] << " radii = " << parcelsm.radius(parcelII) << endl;
        	    }
        	}
		Pout << " Number of parcels = " << parcelsm.numberOfParticles() << endl;
		Pout << " Local number of parcels = " << parcelI << endl;
        	parcelList_.resize(NparcelCellID);
		Pout << " NparcelCellID = " << NparcelCellID << endl;
        	Pout << " Parcel list resize = " << parcelList_.size() << endl;
		Pout << " parcelList_.size() - parcelI =  " << parcelList_.size() - parcelI << endl; 
		
		// We lose parcel somewhere
		if ( parcelList_.size() - parcelI > 0 ) 
		{
			parcelList_.resize(parcelI);
			Pout << " parcelList_.size() - parcelI =  " << parcelList_.size() - parcelI << endl; 
		}
        
		//if ( verbose_ && nparticle_ <=2)
		if ( verbose_ && nparticle_ <=8)
		{							
			int index = 0;
			labelList parcelPart = parcelParticleList_[index];
			Pout << " Parcel particle list " << parcelPart << endl; 
			Pout << " Parcel center     " <<  parcelsm.position(index)  		<< endl;
			Pout << " Parcel velocity   " <<  parcelsm.velocity(index)  		<< endl;
			Pout << " Parcel diameter   " << 2.*parcelsm.radius(index)              << endl;
			Pout << " Parcel cell ID    " <<  parcelsm.cellIDs()[index][0]          << endl;
					
			forAll(parcelPart, ii)
			{
				Pout << " Particle " << parcelPart[ii] << endl;
				Pout << " Particle center    " <<        sm.position(parcelPart[ii])   << endl;
				Pout << " Particle velocity  " <<        sm.velocity(parcelPart[ii])   << endl;
				Pout << " Particle diameter  " <<     2.*sm.radius(parcelPart[ii])     << endl;
			        Pout << " Particle cell ID   " <<        sm.cellIDs()[parcelPart[ii]][0]        << endl;
			}			
		}
		
		// Deallocate kdTree
		//delete kdTree;
        
        // Deallocate memories
        //parcelsm.dataExchangeM().destroy(parcelPositions,3);
        //parcelsm.dataExchangeM().destroy(parcelVelocities,3);
        //parcelsm.dataExchangeM().destroy(parcelRadii,1);
        //parcelsm.dataExchangeM().destroy(parcelCellID,1);
        
	}
 

}	
