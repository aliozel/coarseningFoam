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
    standardSearchPost

SourceFiles
    standardSearchPost.C

\*---------------------------------------------------------------------------*/

#ifndef standardSearchPost_H
#define standardSearchPost_H

#include "locateModelPost.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class standardSearchPost
:
    public locateModelPost
{

public:

    //- Runtime type information
    TypeName("standardPost");


    // Constructors

        //- Construct from components
        standardSearchPost
        (
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Destructor

        ~standardSearchPost();


    // Member Functions
        label findCell
        (
	    double** const& mask,
	    //double**& positions,
	    vectorField& positions,
	    //double**& cellIDs,
	    labelList&	cellIDs,
	    //int size
	    label size,
	    scalarField& radius
        ) const;

        label findSingleCell
        (
            vector& position,
            label& oldCellID
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
