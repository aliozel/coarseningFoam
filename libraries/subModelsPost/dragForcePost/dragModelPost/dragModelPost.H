/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    Locate model

Class
    dragModelPost

SourceFiles
    dragModelPost.C

-----------------------------------------------------------------------------*/

#ifndef dragModelPost_H
#define dragModelPost_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "parcelCloud.H"

#if defined(version21) || defined(version16ext)
    #include "turbulenceModel.H"
#elif defined(version15)
    #include "RASModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dragModelPost Declaration
\*---------------------------------------------------------------------------*/

class dragModelPost
{

protected:

    // Protected data
    const dictionary& dict_;

    const objectRegistry& obr_;

    mutable vectorField dragForce_;
    
    mutable vectorField gradPgForce_;    

    mutable scalarField taup_;

    mutable vectorField Ufxp_;
    
    mutable labelList cellIDs_;

    mutable vectorField velocities_;
        
    mutable vectorField positions_;
    
    mutable scalar nP_;
    
    mutable scalarField voidfractionxp_; 
    
    mutable vectorField dragForceEulerian_; 
    
    mutable vectorField gradPgForceEulerian_;
    
    mutable scalarField taupEulerian_;      

    //- Turbulence model
    const word turbulenceModelType_;

    #if defined(version21) || defined(version16ext)
        #ifdef comp
            const compressible::turbulenceModel& turbulence_;
        #else
            const incompressible::turbulenceModel& turbulence_;
        #endif
    #elif defined(version15)
        const incompressible::RASModel& turbulence_;
    #endif 
    
    mutable scalar rhoParticle_;
    	
public:

    //- Runtime type information
    TypeName("dragModelPost");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dragModelPost,
        dictionary,
        (
            const dictionary& dict,
            const objectRegistry& obr
        ),
        (dict,obr)
    );

    // Constructors

    //- Construct from components
    dragModelPost
    (
        const dictionary& dict,
        const objectRegistry& obr
    );
    
    // Destructor

    virtual ~dragModelPost();

    // Selector

    static autoPtr<dragModelPost> New
    (
        const dictionary& dict,
        const objectRegistry& obr
    );

    // Member Functions
    virtual void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_	    
    ) const = 0;

    virtual void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	vectorField& partPositions_,
	vectorField& partVelocities_,
	scalarList& partRadii_	
    ) const = 0;

    virtual void setForce
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	List<vector>& partVelocities_,
	List<scalar>& partRadii_	
    ) const = 0;
        
    virtual void setForceParcel
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	parcelCloud& parcelCloud_	    
    ) const = 0;
 
    virtual void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_,
	const volVectorField& Us_	    
    ) const = 0;

    virtual void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	scalarField& partRadii_,
	const volVectorField& Us_	    
    ) const = 0;

    virtual void setForceEulerian
    (
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	List<scalar>& partRadii_,
	const volVectorField& Us_	    
    ) const = 0;

    //- Grad Pg
    virtual void setForce
    (
    	const volVectorField& gradPg_,
	cfdemCloud& particleCloud_
    ) const = 0;

    virtual void setForce
    (
	const volVectorField& gradPg_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	scalarList& partRadii_
    ) const = 0;
            
    virtual void setForceParcel
    (
    	const volVectorField& gradPg_,
	parcelCloud& parcelCloud_
    ) const = 0;

    virtual void setForceEulerian
    (
	const volScalarField& voidfraction_, 
	const volVectorField& gradPg_	    
    ) const = 0;

    //- New functions with weigthing 
    virtual void setForce
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
    ) const = 0;   

    virtual void setForce
    (
	const volVectorField& gradPg_,
	const int & partNP_,
	labelList& partCellIDs_,
	List<vector>& partPositions_,
	scalarList& partRadii_,
	List<List<scalar> >& particleWeights_,
	labelListList& neighboringCellIDs_	
    ) const = 0;
                
    inline vector& Fp(const int& index) const{ return dragForce_[index]; };

    inline vector& gradPgFp(const int& index) const{ return gradPgForce_[index]; };
    
    inline scalar& taup(const int& index) const{ return taup_[index]; };

    inline vector& Ufxp(const int& index) const{ return Ufxp_[index]; };
    
    inline label&  cellID(const int& index) const{ return cellIDs_[index]; };

    inline vector& xp(const int& index) const{ return positions_[index]; };
        
    inline vector& Up(const int& index) const{ return velocities_[index]; };
    
    inline scalar& nP() const{ return nP_; };
    
    inline scalar& voidfractionxp(const int& index) const{ return voidfractionxp_[index]; };

    inline vector& FpEulerian(const int& index) const{ return dragForceEulerian_[index]; };

    inline vector& gradPgFpEulerian(const int& index) const{ return gradPgForceEulerian_[index]; };

    inline scalar& taupEulerian(const int& index) const{ return taupEulerian_[index]; };
    
    #if defined(version21) || defined(version16ext)
         #ifdef comp
             inline const compressible::turbulenceModel& turbulence() const{ return turbulence_; };
         #else
             inline const incompressible::turbulenceModel& turbulence() const{ return turbulence_; };
         #endif
    #elif defined(version15)
        inline const incompressible::RASModel& turbulence() const{ return turbulence_; };
    #endif 

    inline scalar& rhoParticle() const{ return rhoParticle_; };
    
    inline tmp<volScalarField> nuf() const{ return turbulence().nu();};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
