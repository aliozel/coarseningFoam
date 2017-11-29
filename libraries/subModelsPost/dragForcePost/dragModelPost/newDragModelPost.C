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

-----------------------------------------------------------------------------*/

#include "error.H"

#include "dragModelPost.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<dragModelPost> dragModelPost::New
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    word dragModelPostType
    (
        dict.lookup("dragModel")
    );

    Info << tab << "Selecting dragModelPost "
                << dragModelPostType << endl;


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dragModelPostType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dragModelPost::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown dragModelPostType type "
            << dragModelPostType
            << ", constructor not in hash table" << endl << endl
            << "    Valid dragModelPost types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<dragModelPost>(cstrIter()(dict,obr));
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
