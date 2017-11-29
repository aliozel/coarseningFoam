/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    Search

Class
    engineSearchPost

SourceFiles
    engineSearchPost.C

-----------------------------------------------------------------------------*/
#ifndef engineSearchPost_H
#define engineSearchPost_H

#include "locateModelPost.H"

#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class engineSearchPost
:
    public locateModelPost
{

private:

    dictionary propsDict_;

    Switch faceDecomp_;

    Switch treeSearch_;

    meshSearch searchEngine_;

public:

    //- Runtime type information
    TypeName("enginePost");


    // Constructors

        //- Construct from components
        engineSearchPost
        (
            const dictionary& dict,
            const objectRegistry& obr
        );

    // Destructor

        ~engineSearchPost();


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
