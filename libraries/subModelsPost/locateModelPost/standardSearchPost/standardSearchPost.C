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

#include "error.H"
#include "standardSearchPost.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardSearchPost, 0);

addToRunTimeSelectionTable
(
    locateModelPost,
    standardSearchPost,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardSearchPost::standardSearchPost
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    locateModelPost(dict,obr)//,
    //obr_(obr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardSearchPost::~standardSearchPost()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label standardSearchPost::findCell
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
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    vector position;
    for(int index = 0;index < size; ++index)
    {

        //cellIDs[index][0]=-1;
	cellIDs[index] =-1;

        //if(mask[index][0] && particleCloud_.radius(index) > SMALL)
        //if(particleCloud_.radius(index) > SMALL)
        if(radius[index] > SMALL)
	{
            // create pos vector
            for(int i=0;i<3;i++) position[i] = positions[index][i];

            // find cell
            #ifdef version16ext
                //cellIDs[index][0] = particleCloud_.mesh().findCell(position);
            	//cellIDs[index] = particleCloud_.mesh().findCell(position);
	        cellIDs[index] = mesh.findCell(position);
	    #elif defined(version21)
                //cellIDs[index][0] = particleCloud_.mesh().findCell(position, polyMesh::FACEPLANES);
            	//cellIDs[index] = particleCloud_.mesh().findCell(position, polyMesh::FACEPLANES);
            	cellIDs[index] = mesh.findCell(position, polyMesh::FACEPLANES);
	    #endif
        }
    }

    return 1;
}

label standardSearchPost::findSingleCell
(
    vector& position,
    label& oldCellID
) const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    // find cell
    #ifdef version16ext
        //return particleCloud_.mesh().findCell(position);
        return mesh.findCell(position);
    #elif defined(version21)
        //return particleCloud_.mesh().findCell(position, polyMesh::FACEPLANES);
        return mesh.findCell(position, polyMesh::FACEPLANES);
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
