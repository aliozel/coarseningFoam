/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    convolution

Description
    Function object that calculate convolution

SourceFiles
    convolution.C

-----------------------------------------------------------------------------*/

#ifndef convolution_H
#define convolution_H

#include "dictionary.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "Ostream.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "volFields.H"

#include "constructFilter.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class fvMesh;
class regIOobject;
class constructFilter;

class convolution
{
    //- Disallow default bitwise copy construct
    convolution(const convolution&);

    //- Disallow default bitwise assignment
    void operator=(const convolution&);

protected:

    // Private data

    // Mesh
    const objectRegistry& obr_;
    
    // Filter types
    const constructFilter& filter_;
    
    // Read from file
    bool readFromFile_;
    bool readFilteredVars_;

    // Init
    void init();
	
public:

    //- Runtime type information
    TypeName("convolution");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    convolution
    (
       const objectRegistry& obr,
       const constructFilter& box,
       const bool& readFilteredVars
    );

    //- Destructor
    virtual ~convolution();

    // Member Functions

    //- Calculate convolution
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

    //- Calculate mean in filter volume 
    virtual volScalarField coarseningSum
    (
	const word& name,
	const int& filterWidth,
	const volScalarField& phi	
    );

    //- Write stencils into the file
    virtual void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
