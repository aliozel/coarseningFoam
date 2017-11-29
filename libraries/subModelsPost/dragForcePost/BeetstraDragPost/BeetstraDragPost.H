/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    Beetstra drag model

Class
    BeetstraDragPost

SourceFiles
    BeetstraDragPost.C

\*---------------------------------------------------------------------------*/

#ifndef BeetstraDragPost_H
#define BeetstraDragPost_H

#include "dragModelPost.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class BeetstraDragPost
:
    public dragModelPost
{
    
    dictionary propsDict_;

    scalar rhoParticle_;
    
    word densityFieldName_; 
    
    const volScalarField& rho_;

    bool verbose_;

    bool interpolation_; 
    
    scalar maxCellPerPart_;
    
    bool taupFuncPhi_;
    
    bool taupConst_;
    
public:

    //- Runtime type information
    TypeName("BeetstraDragPost");


    // Constructors

     //- Construct from components
     BeetstraDragPost
     (
         const dictionary& dict,
         const objectRegistry& obr
     );

    // Destructor

    ~BeetstraDragPost();

    // Member Functions
    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_
    ) const; 

    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	vectorField& partPositions_,
	vectorField& partVelocities_,
	scalarList& partRadii_
    ) const;

    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	List<vector>& partVelocities_,
	List<scalar>& partRadii_
    ) const;
        
    void setForceParcel
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	parcelCloud& parcelCloud_
    ) const;     

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_,
	const volVectorField& Us_	
    ) const;

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	scalarField& partRadii_,
	const volVectorField& Us_	
    ) const;

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	List<scalar>& partRadii_,
	const volVectorField& Us_	
    ) const;

    void setForce
    (
    	const volVectorField& gradPg_,
	cfdemCloud& particleCloud_
    ) const; 

    void setForce
    (
	const volVectorField& gradPg_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	scalarList& partRadii_
    ) const;
            
    void setForceParcel
    (
    	const volVectorField& gradPg_,
	parcelCloud& parcelCloud_
    ) const;
    
    void setForceEulerian
    (
	const volScalarField& voidfraction_, 
	const volVectorField& gradPg_	    
    ) const;   
    
    //- New functions with weigthing 
    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	List<vector>& partVelocities_,
	List<scalar>& partRadii_,
	List<List<scalar> >& particleWeights_,
	labelListList& neighboringCellIDs_
    ) const;     

    void setForce
    (
	const volVectorField& gradPg_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	scalarList& partRadii_,
	List<List<scalar> >& particleWeights_,
	labelListList& neighboringCellIDs_	
    ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
