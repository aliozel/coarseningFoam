{
    scalar pi = M_PI;
    vector position(0,0,0);
    label cellID(-1);
    scalar radius(-1.);
    scalar cellVol(0);
    scalar tooMuch(0.);
    
    // Lists for mapping function
    cellsPerParticle.setSize(nP[me]);
    neighboringCellIDs.setSize(nP[me],subCellIDsList);
    
    particleWeights.setSize(nP[me],subCellList);
    particleVolumes.setSize(nP[me],subCellList);
    
    if(IDWeighthing)
    {	
            #include "IDWeighthing.H"
    }
    else if(GaussWeighthing)
    { 	
	    #include "GaussWeighthing.H"
    }
    else if(noWeighthing)
    { 	
	    for(label index = 0; index < nP[me]; index++)
	    {

        	    cellsPerParticle[index] = 1.;
        	    cellID = cellIDs[index];

        	    if (cellID > -1)  
        	    {
                	cellVol = U.mesh().V()[cellID];
        	    	radius = radii[me][index];
        	    	scalar volume =  4./3.*pi*radius*radius*radius;;

                	// Cell centered weight
                	particleWeights[index][0] = 1.;
			
			// store particleVolume for each particle
                	particleVolumes[index][0] = volume;

                	mappedVoidfraction[cellID] -= volume/cellVol;
        	    }
	    }
    }
    
    Info << tab << "Writing mapped volume fraction... " << endl;
    mappedVoidfraction.write();
    // Update mappedVoidfraction
    //forAll(mesh.cells(),index) mappedVoidfraction[index] = 1.;  

}
