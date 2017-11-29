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

-----------------------------------------------------------------------------*/

#include "error.H"

#include "locateModelPost.H"
#include "standardSearchPost.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<locateModelPost> locateModelPost::New
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    word locateModelPostType
    (
        dict.lookup("locateModelPost")
    );

    Info << tab << "Selecting locateModelPost "
                << locateModelPostType << endl;


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(locateModelPostType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "locateModelPost::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown locateModelPostType type "
            << locateModelPostType
            << ", constructor not in hash table" << endl << endl
            << "    Valid locateModelPost types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<locateModelPost>(cstrIter()(dict,obr));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
