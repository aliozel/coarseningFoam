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
    parcelCloud.C

-----------------------------------------------------------------------------*/

#ifndef parcelCloud_H
#define parcelCloud_H

#include "dictionary.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "Ostream.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "volFields.H"

#if defined(version21) || defined(version16ext)
    #include "turbulenceModel.H"
#elif defined(version15)
    #include "RASModel.H"
#endif

#include "cfdemCloud.H"
// Neighboring algorithm 
#include "ANN.h"

#include "meshSearch.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class fvMesh;
class regIOobject;
class ANN;
class locateModelPost;

class parcelCloud
{
    //- Disallow default bitwise copy construct
    parcelCloud(const parcelCloud&);

    //- Disallow default bitwise assignment
    void operator=(const parcelCloud&);

protected:

    // Private data

    // Parcel dictionary 
    const dictionary& dict_;
    	
    // Mesh
    const objectRegistry& obr_;
    
    // Particles' variables
    const int  partNP_;
    List<vector>& partPositions_;
    List<vector>& partVelocities_;
    List<scalar>& partRadii_;	
	
    // Number of particles in parcel
    const int& nPInParcel_;    
    
    // Parcel particles, velocities, positions, types, cellIDs, turbulence model
    labelListList particlesIn_;
    vectorField velocities_;
    vectorField positions_;
    scalarField radius_;
    labelList types_;
    labelList cellIDs_;
    
    // Use all particle for parcel
    bool useAllParticles_;

    // Init
    void init();

    //- Find the closest particles 
    virtual labelListList findParticlesIn();
    
    //- Mesh volume
    int nPParcel_;
    
    //- Locate model
    Switch treeSearch_;	

    //- Turbulence model
    const word turbulenceModelType_;

    #if defined(version21) || defined(version16ext)
        #ifdef comp
            const compressible::turbulenceModel& turbulence_;
        #else
            const incompressible::turbulenceModel& turbulence_;
        #endif
    #elif defined(version15)
        const incompressible::RASModel& turbulence_;
    #endif 	

    //- Locate model
    inline const locateModelPost& locateM() const;

    //- Calculate parcel radius
    virtual scalarField parcelRadius();        
    
    //- Calculate parcel velocity
    virtual vectorField parcelVel();    
    
    //- Calculate parcel position
    virtual vectorField parcelPos(); 
 
    //- Locate parcels
    virtual labelList parcelLocate(); 
        
    //- Dictionary for parcel
    dictionary parcelCoarseningDict_;  

public:

    //- Runtime type information
    TypeName("parcelCloud");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    parcelCloud
    (
    	const dictionary& dict,
	const objectRegistry& obr,
	const int & partNP,
	List<vector>& partPositions,
	List<vector>& partVelocities,
	List<scalar>& partRadii,	
	const int& nPInParcel
    );

    //- Destructor
    virtual ~parcelCloud();

    // Member Functions
    //- Create parcel cloud
    virtual void createParcel();

    //- Read parcel cloud
    virtual void readParcel();
                    
    //- Read filtered variable from file
    virtual void read();    

    //- Write stencils into the file
    virtual void write();
    
    // Access to parcel data
    int numberOfParticles();
    vector velocity(int&);
    vector position(int&);
    scalar radius(int&);
    label type(int&);
    label cellIDs(int&);
    labelList particlesInParcel(int& index);
    
    #if defined(version21) || defined(version16ext)
         #ifdef comp
             inline const compressible::turbulenceModel& turbulence() const;
         #else
             inline const incompressible::turbulenceModel& turbulence() const;
         #endif
    #elif defined(version15)
        inline const incompressible::RASModel& turbulence() const;
    #endif    

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parcelCloudI.H"

#endif
