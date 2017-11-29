// Filtering solid volume fraction
//barAlpp = convKernel.coarseningSca("barAlpp",filterWidth[fWidth],1-voidfraction);
barAlpp = convKernel.coarseningSca("barAlpp",filterWidth[fWidth],1-mappedVoidfraction);

// Filtering solid velocity
//tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-voidfraction,Us);
tildeUs = convKernel.coarseningVecFavre("tildeUs",filterWidth[fWidth],1.-mappedVoidfraction,mappedUs);

// Reynolds stresses
//volScalarField tildeAlppUsxUsx = convKernel.coarseningScaFavre("tildeAlppUsxUsx",filterWidth[fWidth],(1-voidfraction),Us.component(0)*Us.component(0));
//volScalarField tildeAlppUsyUsy = convKernel.coarseningScaFavre("tildeAlppUsyUsy",filterWidth[fWidth],(1-voidfraction),Us.component(1)*Us.component(1));
//volScalarField tildeAlppUszUsz = convKernel.coarseningScaFavre("tildeAlppUszUsz",filterWidth[fWidth],(1-voidfraction),Us.component(2)*Us.component(2));
volScalarField tildeAlppUsxUsx = convKernel.coarseningScaFavre("tildeAlppUsxUsx",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(0));
volScalarField tildeAlppUsyUsy = convKernel.coarseningScaFavre("tildeAlppUsyUsy",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(1)*mappedUs.component(1));
volScalarField tildeAlppUszUsz = convKernel.coarseningScaFavre("tildeAlppUszUsz",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(2)*mappedUs.component(2));

volScalarField tildeAlppUsxUsy = convKernel.coarseningScaFavre("tildeAlppUsxUsy",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(1));
volScalarField tildeAlppUsxUsz = convKernel.coarseningScaFavre("tildeAlppUsxUsz",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(0)*mappedUs.component(2));
volScalarField tildeAlppUsyUsz = convKernel.coarseningScaFavre("tildeAlppUsyUsz",filterWidth[fWidth],(1-mappedVoidfraction),mappedUs.component(1)*mappedUs.component(2));

forAll(mesh.cells(),cellI)
{
	diffTildeeAlppUsUs[cellI] = symmTensor
					(		
			  		  tildeAlppUsxUsx[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][0] ),
			  		  tildeAlppUsyUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][1] ),
			  		  tildeAlppUszUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][2]*tildeUs[cellI][2] ),
			  		  tildeAlppUsxUsy[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][1] ),
			  		  tildeAlppUsxUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][0]*tildeUs[cellI][2] ),
			  		  tildeAlppUsyUsz[cellI] - barAlpp[cellI]*( tildeUs[cellI][1]*tildeUs[cellI][2] )
					);
}

// Calculate Umean in a filter control volume
volScalarField nPmean = convKernel.coarseningSum("nPmeanFilter",filterWidth[fWidth],nPCell);
volScalarField UxmeanFilter = convKernel.coarseningSum("UxmeanFilter",filterWidth[fWidth],nPCell*meanVelCell.component(0));
volScalarField UymeanFilter = convKernel.coarseningSum("UymeanFilter",filterWidth[fWidth],nPCell*meanVelCell.component(1));
volScalarField UzmeanFilter = convKernel.coarseningSum("UzmeanFilter",filterWidth[fWidth],nPCell*meanVelCell.component(2));

//- Normalize
forAll(mesh.cells(),cellI)
{
    if(nPmean[cellI]!=0)
    {
        UxmeanFilter[cellI] /= nPmean[cellI];
	UymeanFilter[cellI] /= nPmean[cellI];
	UzmeanFilter[cellI] /= nPmean[cellI];
    }else
    {
        UxmeanFilter[cellI] = 0.;
	UymeanFilter[cellI] = 0.;
	UzmeanFilter[cellI] = 0.;
    } 	
}

//-Calculate stresses

volSymmTensorField tildeSigmaKin
(
    IOobject
    (
	"tildeSigmaKin"+filterWidthName,
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
        UsxInt = UxmeanFilter[cellI];
        UsyInt = UymeanFilter[cellI];
        UszInt = UzmeanFilter[cellI];
      }else
      {
        // Centering cell
        weightP = particleWeights[index][0];
	UsxInt += weightP*tildeUs[cellI][0];
        UsyInt += weightP*tildeUs[cellI][1];
        UszInt += weightP*tildeUs[cellI][2]; 

	if( index==0 ) Info << tab << "index = " << index << " weightP = " << weightP;

	// Neighboring cells    
	for(int subCell=1 ; subCell < maxCellPerPart; subCell++)
	{
            label cellJ = neighboringCellIDs[index][subCell];
            if (cellJ > -1)
            {
        	weightP = particleWeights[index][subCell];
		if( index == 0 ) Info << tab << " weightP-2 = " << weightP;
		UsxInt += weightP*tildeUs[cellJ][0];
                UsyInt += weightP*tildeUs[cellJ][1];
                UszInt += weightP*tildeUs[cellJ][2]; 	      		
	    }
	}
		           
      } 


      if( index == 0 ) Info << endl;      

      // Particle kinetic stress in each cell

      tildeSigmaKin[cellI].xx() += mass_p * ( velPart[0] - UsxInt )
                                	  * ( velPart[0] - UsxInt );

      tildeSigmaKin[cellI].yy() += mass_p * ( velPart[1] - UsyInt )
                                	  * ( velPart[1] - UsyInt );

      tildeSigmaKin[cellI].zz() += mass_p * ( velPart[2] - UszInt )
                                	  * ( velPart[2] - UszInt );

      tildeSigmaKin[cellI].xy() += mass_p * ( velPart[0] - UsxInt )
                                	  * ( velPart[1] - UsyInt );

      tildeSigmaKin[cellI].xz() += mass_p * ( velPart[0] - UsxInt )
                                	  * ( velPart[2] - UszInt );

      tildeSigmaKin[cellI].yz() += mass_p * ( velPart[1] - UsyInt )
                                	  * ( velPart[2] - UszInt );					  
					
   }   
} 

// Normalize kinetic stresses at each cell
forAll(mesh.cells(),cellI)
{
    for(int idir=0; idir<6; idir++) 
    	tildeSigmaKin[cellI][idir]/= mesh.V()[cellI];    
}
