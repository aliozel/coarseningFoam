/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Parcel

Description
    Function object that agglomerateCloud

SourceFiles
    agglomerateCloudI.H

-----------------------------------------------------------------------------*/

#include "locateModelPost.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const locateModelPost& agglomerateCloud::locateM() const
{
    return locateModel_;
}

#if defined(version21) || defined(version16ext)
    #ifdef comp
        inline const compressible::turbulenceModel& agglomerateCloud::turbulence() const
    #else
        inline const incompressible::turbulenceModel& agglomerateCloud::turbulence() const
    #endif
#elif defined(version15)
    inline const incompressible::RASModel& agglomerateCloud::turbulence() const
#endif
{
    return turbulence_;
}

}
// ************************************************************************* //
