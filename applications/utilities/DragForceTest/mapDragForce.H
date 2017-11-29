{
    label cellI(-1);
    vector valueVec(0,0,0);
    vector valueVec2(0,0,0);
    scalar weightP(0.);
    
    bool debug(false);
    //bool debug(true);

    if(!usingLocalEul)
    {
	 for(int index=0; index < nP[me]; index++)
	 {
	     // Centering cell
	     cellI = cellIDs[index];
             if (cellI > -1)
             {
        	 valueVec = dragF().Fp(index);
        	 valueVec2 = dragF().gradPgFp(index);
		 weightP = particleWeights[index][0];	    

        	 if(filteringWithGradPg)
		 {
	    	     dragFpEulerian[cellI] += ( valueVec + valueVec2 ) * weightP;
		     if(cellI==0) Info << tab << "dragFpEulerian = " 
	    				      << ( valueVec + valueVec2 ) * weightP / mesh.V()[cellI]
	                        	      << endl;	
        	 }else
		 {
	    	     dragFpEulerian[cellI] += valueVec*weightP;
		     if(cellI==0) 
		     {
			     Info << tab << "dragFpEulerian = " 
	    		  		 << valueVec*weightP / mesh.V()[cellI]
					 << " Ufxp = " << dragF().Ufxp(index)
					 << " Up = " << dragF().Up(index)
					 << " voidfractionxp = " << dragF().voidfractionxp(index)
					 << " weightP-1 = " << weightP
	                        	 << endl;
		     }		    	
		 }

	     }

	     if(debug) Info << " index = " << index << " cellI = " << cellI << " weightP = " << weightP;

	     // Neighboring cells    
	     for(int subCell=1 ; subCell < maxCellPerPart; subCell++)
	     {
        	 cellI = neighboringCellIDs[index][subCell];
        	 if (cellI > -1)
        	 {
        	     valueVec = dragF().Fp(index);
        	     valueVec2 = dragF().gradPgFp(index); 
		     weightP = particleWeights[index][subCell];

		     if(filteringWithGradPg)
	    	     {
	    		     dragFpEulerian[cellI] += ( valueVec + valueVec2 ) * weightP;
            	     }else
	    	     {
	    		     dragFpEulerian[cellI] += valueVec*weightP;
			     if(cellI==0) 
			     {
				     Info << tab << "dragFpEulerian = " 
	    		  			 << valueVec*weightP / mesh.V()[cellI]
						 << " Ufxp = " << dragF().Ufxp(index)
						 << " Up = " << dragF().Up(index)
						 << " voidfractionxp = " << dragF().voidfractionxp(index)
						 << " weightP-2 = " << weightP
	                        		 << endl;
			     }	
	    	     }

 	             if(debug) Info << " weight2 = " << weightP;
        	 }
	     }
	     if(debug) Info << endl;
	     

	}
	forAll(mesh.cells(),index)
	{
	     // Divide force by volume of mesh
	     dragFpEulerian[index] /= mesh.V()[index];
        }		
    }else
    {
	forAll(mesh.cells(),index)
	{
	   dragFpEulerian[index] = dragFEulerian().FpEulerian(index);
        }
    }		
    
    Info << tab << "Force was mapped... " << endl;
    Info << tab << "cellI = 0 " << endl;
    Info << tab << "dragFpEulerian = " << dragFpEulerian[0] << endl;
    Info << tab << "alpp = " << 1.-mappedVoidfraction[0] << endl;
    Info << tab << "U = " << U[0] << endl;
    Info << tab << "mappedUs = " << mappedUs[0] << endl;
    Info << tab << "" << endl;
}    
