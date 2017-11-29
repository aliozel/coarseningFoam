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
    Function object that parcelCloud

SourceFiles
    parcelCloudI.H

-----------------------------------------------------------------------------*/

#include "locateModelPost.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if defined(version21) || defined(version16ext)
    #ifdef comp
        inline const compressible::turbulenceModel& parcelCloud::turbulence() const
    #else
        inline const incompressible::turbulenceModel& parcelCloud::turbulence() const
    #endif
#elif defined(version15)
    inline const incompressible::RASModel& parcelCloud::turbulence() const
#endif
{
    return turbulence_;
}

}
// ************************************************************************* //
