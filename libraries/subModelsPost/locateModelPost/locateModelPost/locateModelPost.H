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
    locateModelPost

SourceFiles
    localModelPost.C

-----------------------------------------------------------------------------*/

#ifndef locateModelPost_H
#define locateModelPost_H

#include "fvCFD.H"
#include "parcelCloud.H"
//#include "agglomerateCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class locateModelPost Declaration
\*---------------------------------------------------------------------------*/

class locateModelPost
{

protected:

    // Protected data
        const dictionary& dict_;

        const objectRegistry& obr_;

public:

    //- Runtime type information
    TypeName("locateModelPost");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            locateModelPost,
            dictionary,
            (
                const dictionary& dict,
                const objectRegistry& obr
            ),
            (dict,obr)
        );


    // Constructors

        //- Construct from components
        locateModelPost
        (
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Destructor

        virtual ~locateModelPost();


    // Selector

        static autoPtr<locateModelPost> New
        (
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions
        virtual label findCell
        (
	    double** const& mask,
	    //double**& positions,
	    vectorField& positions,
	    //double**& cellIDs,
	    labelList&	cellIDs,
	    //int size
	    label size,
	    scalarField& radius
        ) const = 0;

        virtual label findSingleCell
        (
            vector& position,
            label& oldCellID
        ) const = 0;

        virtual label intersection
        (
                const point& pStart,
                const point& pEnd
        ) const {return -1;};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
