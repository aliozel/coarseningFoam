/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Calculate correlation coefficient

-----------------------------------------------------------------------------*/

#include "correlationCoeff.H"

namespace Foam {
    defineTypeNameAndDebug(correlationCoeff,0);
}

Foam::correlationCoeff::correlationCoeff
(
    const objectRegistry& obr
)
:
    	obr_(obr),
	readFromFile_(false)
{
	// Init
	init();
}

Foam::correlationCoeff::~correlationCoeff()
{}

void Foam::correlationCoeff::init()
{
	// Do nothing
}

Foam::scalar Foam::correlationCoeff::calc
(
	const volScalarField& phiA,
	const volScalarField& phiB
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	scalar r(0);
	scalar meanA(0);
	scalar meanB(0);		
	
	// Calculate mean
	forAll(mesh.cells(),cellI)
	{
		meanA += phiA[cellI];
		meanB += phiB[cellI];
	} 
	
	// Normalize
	meanA /= mesh.cells().size();
	meanB /= mesh.cells().size();

	// Calculate the coefficient
	scalar rE(0);
	scalar rStdA(0);
	scalar rStdB(0);
	
	forAll(mesh.cells(),cellI)
	{
		rE    += ( phiA[cellI] - meanA ) * ( phiB[cellI] - meanB);
	
		rStdA += ( phiA[cellI] - meanA ) * ( phiA[cellI] - meanA );
		rStdB += ( phiB[cellI] - meanB ) * ( phiB[cellI] - meanB );		
	}	
	
	if(rStdA != 0 && rStdB != 0 )
	{
		r = rE / ( sqrt(rStdA)*sqrt(rStdB) );
	}else
	{
		r = 0;
	}
	
	return r;		
}

Foam::scalar Foam::correlationCoeff::calc
(
	const scalarField& phiA,
	const scalarField& phiB
)
{
	scalar r(0);
	scalar meanA(0);
	scalar meanB(0);		
	
	// Calculate mean
	forAll(phiA,cellI)
	{
		meanA += phiA[cellI];
		meanB += phiB[cellI];
	} 
	
	// Normalize
	meanA /= phiA.size();
	meanB /= phiB.size();

	// Calculate the coefficient
	scalar rE(0);
	scalar rStdA(0);
	scalar rStdB(0);
	
	forAll(phiA,cellI)
	{
		rE    += ( phiA[cellI] - meanA ) * ( phiB[cellI] - meanB);
	
		rStdA += ( phiA[cellI] - meanA ) * ( phiA[cellI] - meanA );
		rStdB += ( phiB[cellI] - meanB ) * ( phiB[cellI] - meanB );		
	}	
	
	if(rStdA != 0 && rStdB != 0 )
	{
		r = rE / ( sqrt(rStdA)*sqrt(rStdB) );
	}else
	{
		r = 0;
	}
	
	return r;		
}
	

void Foam::correlationCoeff::read()
{
    if(readFromFile_) 
    {
	// Do nothing
    }

}

void Foam::correlationCoeff::write()
{
    if(readFromFile_) 
    {
	// Do nothing
    }

}


// ************************************************************************* //
