/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    Carrier phase pressure gradient

Class
     Carrier phase pressure gradient

SourceFiles
    gradPgPost.C

-----------------------------------------------------------------------------*/

#ifndef gradPgPost_H
#define gradPgPost_H

#include "dragModelPost.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gradPgPost Declaration
\*---------------------------------------------------------------------------*/

class gradPgPost
:
    public dragModelPost
{
    
    dictionary propsDict_;

    bool verbose_;

    bool interpolation_; 
    
public:

    //- Runtime type information
    TypeName("gradPgPost");


    // Constructors

    //- Construct from components
    gradPgPost
    (
	const dictionary& dict,
	const objectRegistry& obr
    );

    // Destructor

    ~gradPgPost();

    // Member Functions
    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_
    ) const{}; 

    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	vectorField& partPositions_,
	vectorField& partVelocities_,
	scalarList& partRadii_
    ) const{};

    void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	List<vector>& partVelocities_,
	List<scalar>& partRadii_
    ) const{};
        
    void setForceParcel
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	parcelCloud& parcelCloud_
    ) const{};     

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_,
	const volVectorField& Us_	
    ) const{};

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	scalarField& partRadii_,
	const volVectorField& Us_	
    ) const{};

    void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	List<scalar>& partRadii_,
	const volVectorField& Us_	
    ) const{};

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
	scalarList& partRadii_
    ) const;
            
    void setForceParcel
    (
    	const volVectorField& gradPg_,
	parcelCloud& parcelCloud_
    ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
