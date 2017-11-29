/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    mapping

Description
    Function object that calculate mapping

SourceFiles
    mapping.C

-----------------------------------------------------------------------------*/

#ifndef mapping_H
#define mapping_H

#include "dictionary.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "Ostream.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "volFields.H"

#include "boxFilter.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class fvMesh;
class regIOobject;
class boxFilter;

class mapping
{
    //- Disallow default bitwise copy construct
    mapping(const mapping&);

    //- Disallow default bitwise assignment
    void operator=(const mapping&);

protected:

    // Private data

    // Mesh
    const objectRegistry& obr_;
    
    // Filter types
    const boxFilter& filter_;
    
    // Read from file
    bool readFilteredVars_;
    bool readFromFile_;

    // Init
    void init();
	
public:

    //- Runtime type information
    TypeName("mapping");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    mapping
    (
       const objectRegistry& obr,
       const boxFilter& box,
       const bool& readFilteredVars
    );

    //- Destructor
    virtual ~mapping();

    // Member Functions

    //- Calculate mapping
    virtual volScalarField coarseningSca
    (
	const word& name,
	const int& filterWidth,
	const volScalarField& phi	
    );
       	
    virtual volScalarField coarseningScaFavre
    (
	const word& name,
	const int& filterWidth,
	const volScalarField& voidfraction,  	
	const volScalarField& phi	
    );
    
    virtual volVectorField coarseningVec
    (
	const word& name,
	const int& filterWidth,
	const volVectorField& phi	
    );
       	
    virtual volVectorField coarseningVecFavre
    (
	const word& name,
	const int& filterWidth,
	const volScalarField& voidfraction,  	
	const volVectorField& phi	
    );  

    //- Write stencils into the file
    virtual void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
