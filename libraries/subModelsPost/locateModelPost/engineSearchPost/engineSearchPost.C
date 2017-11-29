/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
                 Locate model

-----------------------------------------------------------------------------*/

#include "error.H"
#include "engineSearchPost.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(engineSearchPost, 0);

addToRunTimeSelectionTable
(
    locateModelPost,
    engineSearchPost,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
engineSearchPost::engineSearchPost
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    locateModelPost(dict,obr),
    propsDict_(dict.subDict(typeName + "Props")),
    //faceDecomp_(propsDict_.lookup("faceDecomp")),
    treeSearch_(propsDict_.lookup("treeSearch")),
    #ifdef version16ext
        //searchEngine_(particleCloud_.mesh(),false) //(particleCloud_.mesh(),faceDecomp_)
        searchEngine_(refCast<const fvMesh>(obr),false)
    #elif defined(version21)
        //searchEngine_(particleCloud_.mesh(),polyMesh::FACEPLANES) // FACEPLANES or FACECENTRETETS; FACEDIAGTETS not stable
        searchEngine_(refCast<const fvMesh>(obr),polyMesh::FACEPLANES) 
    #endif
    //searchEngine_(particleCloud_.mesh(),faceDecomp_) // only 2.0.x
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

engineSearchPost::~engineSearchPost()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label engineSearchPost::findCell
(
    double** const& mask,
    //double**& positions,
    vectorField& positions,
    //double**& cellIDs,
    labelList&	cellIDs,
    //int size
    label size,
    scalarField& radius
) const
{
    vector position;
    for(int index = 0;index < size; ++index)
    {
        //cellIDs[index][0]=-1;
	cellIDs[index] = -1;
        //if(mask[index][0] && particleCloud_.radius(index) > SMALL)
        //if(particleCloud_.radius(index) > SMALL)
        if(radius[index] > SMALL)
	{

            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];
            // find cell
            //cellIDs[index][0] =searchEngine_.findCell(position,cellIDs[index][0],treeSearch_);
            cellIDs[index] =searchEngine_.findCell(position,cellIDs[index],treeSearch_); 
	}
    }
    return 1;
}

label engineSearchPost::findSingleCell
(
    vector& position,
    label& oldCellID
) const
{
    // find cell
    return searchEngine_.findCell(position,oldCellID,treeSearch_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
