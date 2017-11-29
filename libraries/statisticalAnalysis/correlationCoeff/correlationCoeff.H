/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Correlation coefficient

Description
    Function object that calculate correlation coefficient

SourceFiles
    correlationCoeff.C

-----------------------------------------------------------------------------*/

#ifndef correlationCoeff_H
#define correlationCoeff_H

#include "dictionary.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "Ostream.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "volFields.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class fvMesh;
class regIOobject;

class correlationCoeff
{
    //- Disallow default bitwise copy construct
    correlationCoeff(const correlationCoeff&);

    //- Disallow default bitwise assignment
    void operator=(const correlationCoeff&);

protected:

    // Private data

    // Mesh
    const objectRegistry& obr_;
    
    // Read from file
    bool readFromFile_;

    // Init
    void init();
	
public:

    //- Runtime type information
    TypeName("correlationCoeff");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    correlationCoeff
    (
       const objectRegistry& obr
    );

    //- Destructor
    virtual ~correlationCoeff();

    // Member Functions

    //- Calculate correlationCoeff
    virtual scalar calc
    (
	const volScalarField& phiA,
	const volScalarField& phiB		
    ); 

    virtual scalar calc
    (
	const scalarField& phiA,
	const scalarField& phiB		
    ); 	
    
    //- Read 
    virtual void read();    

    //- Write 
    virtual void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
