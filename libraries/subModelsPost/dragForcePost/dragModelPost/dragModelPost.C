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
    localModelPost.C

-----------------------------------------------------------------------------*/

#include "error.H"
#include "dragModelPost.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dragModelPost, 0);

defineRunTimeSelectionTable(dragModelPost, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dragModelPost::dragModelPost
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    dict_(dict),
    obr_(obr),
    dragForce_(1,vector(0,0,0)),
    gradPgForce_(1,vector(0,0,0)),
    taup_(1,scalar(0)),
    Ufxp_(1,vector(0,0,0)), 
    cellIDs_(1,-1),
    velocities_(1,vector(0,0,0)),    
    positions_(1,vector(0,0,0)),
    nP_(0),
    voidfractionxp_(1,scalar(0)),
    dragForceEulerian_(1,vector(0,0,0)),
    taupEulerian_(1,scalar(0)),
    turbulenceModelType_(dict.lookup("turbulenceModelType")),
    turbulence_
    (
        #if defined(version21) || defined(version16ext)
            #ifdef comp
                refCast<const fvMesh>(obr).lookupObject<compressible::turbulenceModel>
            #else
                refCast<const fvMesh>(obr).lookupObject<incompressible::turbulenceModel>
            #endif
        #elif defined(version15)
            refCast<const fvMesh>(obr).lookupObject<incompressible::RASModel>
        #endif
        (
            turbulenceModelType_
        )
    ),
    rhoParticle_(1.)    	      
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dragModelPost::~dragModelPost()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
