/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Agglomerate

Description
    Function object that agglomerateCloud

SourceFiles
    agglomerateCloud.C

-----------------------------------------------------------------------------*/

#ifndef agglomerateCloud_H
#define agglomerateCloud_H

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
#include "parcelCloud.H"
// Neighboring algorithm 
//#include "ANN.h"

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

class agglomerateCloud
{
    //- Disallow default bitwise copy construct
    agglomerateCloud(const agglomerateCloud&);

    //- Disallow default bitwise assignment
    void operator=(const agglomerateCloud&);

protected:

    // Private data

    // Agglomerate dictionary 
    const dictionary& dict_;
    	
    // Mesh
    const objectRegistry& obr_;
    
    // Particle cloud
    cfdemCloud& sm_;
    
    // Particle list in the sub-domain
    //labelList particlesInSubDomain_;

    // Time-step
    const int count_;
    
    //- Boolean to locate agglomerate
    int locateAgglomerate_;    
    
    // Number of particles in agglomerate
    labelList nPInAgglomerate_;    
    
    // Agglomerate particles, velocities, positions, types, cellIDs, turbulence model
    labelListList particlesIn_;
    vectorField velocities_;
    vectorField positions_;
    scalarField radius_;
    labelList types_;
    labelList cellIDs_;
    	    
    // Read from file
    bool readFromFile_;

    //- Number of agglomerate
    int nPAgglomerate_;
    
    // Connectivity input folder
    dictionary agglomerateVTKPropsDict_;
    word filename_;
    fileName relativePath_;
    const char* charFilename_;

    // Liquid connectivity list
    labelListList connectivityList_;

    // Init
    void init();

    //- Read connectivity list
    virtual labelListList readConnectivityList();
    
    //- Find particles in an agglomerate 
    virtual labelListList findParticlesIn();
   
    //- Locate model
    autoPtr<locateModelPost> locateModel_;	

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
            
    //- Calculate agglomerate velocity
    virtual vectorField agglomerateVel();    
    
    //- Calculate agglomerate position
    virtual vectorField agglomeratePos();   

    //- Calculate agglomerate radius
    virtual scalarField agglomerateRadius();
    
    //- Bridge IDs lenght
    int bridgeIDLen_;
    
public:

    //- Runtime type information
    TypeName("agglomerateCloud");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    agglomerateCloud
    (
       const dictionary& dict,
       const objectRegistry& obr,
       cfdemCloud& sm,
       const int& count,
       bool& locateAgglomerate
       //,
       //const labelList& particlesInSubDomain,
       //const int& nPInAgglomerate
    );

    //- Destructor
    virtual ~agglomerateCloud();

    // Member Functions
    //- Create agglomerate cloud
    virtual void createAgglomerate();
                
    //- Read filtered variable from file
    virtual void read();    

    //- Write stencils into the file
    virtual void write();
    
    // Access to agglomerate data
    int numberOfAgglomerates();
    vector velocity(int&);
    vector position(int&);
    scalar radius(int&);
    label type(int&);
    label cellIDs(int&);
    labelList particlesInAgglomerate(int& index);
    
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

#include "agglomerateCloudI.H"

#endif
