CPCCellToCellStencil wideStencil(mesh);
for(label index = 0; index < nP[me]; index++)
{
    // Debugging, particle ID
    label testCellID_ = exList[0]; //0;

    // Particle radius
    scalar radius = radii[me][index];    
    // Particle volume
    scalar volume = 4./3.*radius*radius*radius*pi;
    // Particle position
    vector positionCenter = positions[me][index];

    // Define weighting function
    scalarField particleWeightsGauss(maxCellPerPart,scalar(0.0));
    // Total weigthings
    scalar sumWeights(0);

    // Weighting of cell center
    label particleCenterCellID = cellIDs[index];
    
    if (particleCenterCellID > -1 )
    {

	    if(verbose && index == testCellID_) 
	    {
		Info << tab << "Gaussian mapping - Particle " << testCellID_ << endl;
	    }

	    // Distance to neighbor cells
	    scalar dist_s(0);

	    // Create neighbor cells list
	    const labelList& cellsNeigh = wideStencil[particleCenterCellID]; //mesh.cellCells()[particleCenterCellID]; //neighbourCells[particleCenterCellID];

	    // Loop over neighborg cells
	    forAll(cellsNeigh,jj)
	    {

		    // Fill particle-neighbor cells list
		    neighboringCellIDs[index][jj] = cellsNeigh[jj];

		    // Find distances between particle and neighbouring cells					
		    dist_s = mag(mesh.C()[cellsNeigh[jj]]-positionCenter)/Foam::pow(mesh.V()[cellsNeigh[jj]],1./3.);

		    if(dist_s <= 0.5)
		    {		
			    particleWeightsGauss[jj] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
		    }
		    else if (dist_s > 0.5 && dist_s <= 1.5)
		    {		
			    particleWeightsGauss[jj] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
		    }
		    else if (dist_s > 1.5 && dist_s <= 2.5)
		    {		
			    particleWeightsGauss[jj] =  pow(2.5-dist_s,4)/24.;
		    }
		    else
		    {		
			    particleWeightsGauss[jj] = 0;
		    }

		    sumWeights += particleWeightsGauss[jj];

		    // debug
                    if(verbose && index == testCellID_) 
		    {
			Info << tab << "neighbors[" << jj << "]" << " ID = " << cellsNeigh[jj]
			     << " gaussian weight = " << particleWeightsGauss[jj] << endl;
		    }	
	    }
	    // debug
            if(verbose && index == testCellID_) Info << tab << "sum of gaussian weights = " << sumWeights << endl;
	    
	    forAll(cellsNeigh,jj)
	    {	

		    if ( sumWeights != 0 )
		    {
                	    // volume occupied in every covered cell
                	    scalar occupiedVolume = volume*particleWeightsGauss[jj]/sumWeights;

                	    // correct volumefraction of centre
                	    mappedVoidfraction[cellsNeigh[jj]] -= occupiedVolume/mesh.V()[cellsNeigh[jj]];

			    if(verbose && index == testCellID_)
                    	    {   
				Info << tab << "void fraction contribution into cell " << cellsNeigh[jj] 
				     << " " << occupiedVolume/mesh.V()[cellsNeigh[jj]] << endl;
                    	    }	

                	    particleWeights[index][jj] = particleWeightsGauss[jj]/sumWeights;
                	    particleVolumes[index][jj] = occupiedVolume;
		    }
		    else
		    {
			    particleWeights[index][jj] = 0.;
                	    particleVolumes[index][jj] = 0.;	
		    }
	    }

    }

}// end loop all particles

Info << tab << "End of Gaussian mapping " << endl;
