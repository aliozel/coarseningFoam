//-Calculate stresses
volSymmTensorField sigmaKin
(
    IOobject
    (
	"sigmaKin",
	mesh.time().timeName(),
	mesh,
	IOobject::NO_READ, //READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
);

vector velPart(0,0,0);

scalar pi = constant::mathematical::pi;

scalar weightP(0);
scalar UsxInt(0);
scalar UsyInt(0);
scalar UszInt(0);

for(int index = 0; index < nP[me]; index++)
{
   label cellI = cellIDs[index];
   
   weightP = 0.;
   UsxInt = 0.;
   UsyInt = 0.;
   UszInt = 0.;
   
   if( cellI > -1 )
   {
      // Particle vel  post
      velPart = velocities[me][index];

      // Particle mass      
      scalar rad = radii[me][index];
      scalar mass_p = 4./3.*rhoParticle*pi*rad*rad*rad;     
            
      if(useUsstar)
      {
        UsxInt = meanVelCell[0][cellI];
        UsyInt = meanVelCell[1][cellI];
        UszInt = meanVelCell[2][cellI];
      }else
      {
        // Centering cell
        weightP = particleWeights[index][0];
	UsxInt += weightP*mappedUs[cellI][0];
        UsyInt += weightP*mappedUs[cellI][1];
        UszInt += weightP*mappedUs[cellI][2]; 

	if( index==0 ) Info << tab << "index = " << index << " weightP = " << weightP;

	// Neighboring cells    
	for(int subCell=1 ; subCell < maxCellPerPart; subCell++)
	{
            label cellJ = neighboringCellIDs[index][subCell];
            if (cellJ > -1)
            {
        	weightP = particleWeights[index][subCell];
		if( index == 0 ) Info << tab << " weightP-2 = " << weightP;
		UsxInt += weightP*mappedUs[cellJ][0];
                UsyInt += weightP*mappedUs[cellJ][1];
                UszInt += weightP*mappedUs[cellJ][2]; 	      		
	    }
	}
		           
      } 


      if( index == 0 ) Info << endl;      

      // Particle kinetic stress in each cell

      sigmaKin[cellI].xx() += mass_p * ( velPart[0] - UsxInt )
                                     * ( velPart[0] - UsxInt );

      sigmaKin[cellI].yy() += mass_p * ( velPart[1] - UsyInt )
                                     * ( velPart[1] - UsyInt );

      sigmaKin[cellI].zz() += mass_p * ( velPart[2] - UszInt )
                                     * ( velPart[2] - UszInt );

      sigmaKin[cellI].xy() += mass_p * ( velPart[0] - UsxInt )
                                     * ( velPart[1] - UsyInt );

      sigmaKin[cellI].xz() += mass_p * ( velPart[0] - UsxInt )
                                     * ( velPart[2] - UszInt );

      sigmaKin[cellI].yz() += mass_p * ( velPart[1] - UsyInt )
                                     * ( velPart[2] - UszInt );					  
					
   }   
} 

// Normalize kinetic stresses at each cell
forAll(mesh.cells(),cellI)
{
    for(int idir=0; idir<6; idir++) 
    	sigmaKin[cellI][idir]/= mesh.V()[cellI];    
}
