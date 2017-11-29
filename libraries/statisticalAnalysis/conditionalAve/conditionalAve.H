/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    conditionalAve

Description
    Function object that calculate conditionalAve

SourceFiles
    conditionalAve.C

-----------------------------------------------------------------------------*/

#ifndef conditionalAve_H
#define conditionalAve_H

#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Ostream.H"
#include "IOmanip.H"
#include "IOobjectList.H"
#include "labelIOList.H"
#include "volFields.H"

// swak4Foam ExpressionField
#include "expressionField.H"

// Drag autoPtr
#include "dragModelPost.H"

// Interpolation scheme
#include "interpolationCellPoint.H"

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

class conditionalAve
{
    //- Disallow default bitwise copy construct
    conditionalAve(const conditionalAve&);

    //- Disallow default bitwise assignment
    void operator=(const conditionalAve&);

protected:

    //- Private data

    //- Dictionary 
    const dictionary& dict_;	

    //- Mesh
    const objectRegistry& obr_;
    
    //- Total number of analyses
    label nVariable_;
    label nAveragingVariable_;
    label nTotalCase_;

    //- Conditional averaging parameters
    dictionary conditionalAveDict_;
    fileName outputPath_;   
    
    //- volScalarField variables 
    mutable wordList variableList;
    mutable wordList averagingVariableList;
    mutable labelList nBin;
    mutable scalarField minCondVar;
    mutable scalarField maxCondVar; 
    mutable Field <Field <scalar> > nReal;       
    mutable Field <Field <scalar> > condAveVar;
 	    
    //- Read from file
    bool readFromFile_;

    //- Boolean
    bool conditionalAveraging_;
    
    //- Init
    void init();

    //- Init scalarField
    void initUserDefine();
    
public:

    //- Runtime type information
    TypeName("conditionalAve");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    conditionalAve
    (
	 const dictionary& dict,
	 const dictionary& subDict,
	 const objectRegistry& obr,
	 const label& nVariable,
	 const label& nAveragingVariable,
	 const label& nTotalCase,
	 const bool& conditionalAveraging 
    );

    //- Destructor
    virtual ~conditionalAve();

    // Member Functions
    //- Calculate conditional averaging for volScalarField
    virtual void calc();

    //- Calculate conditional averaging for volScalarField
    virtual void calc
    (
	word& variableName,
	word& userDefineName,
	wordList& averagingVariableListCoarsening,
	int jj
    );
       
    //- Calculate conditional averaging for unresolved drag
    virtual void calcUserDefine
    (
    	 scalarField&,
	 autoPtr<dragModelPost>&,
         word&,
	 label& userLoopVar
    );

    //- Calculate conditional averaging for unresolved drag
    virtual void calcUserDefine
    (
    	 scalarField&,
	 autoPtr<dragModelPost>&,
         word&,
	 label& userLoopVar,
         List<List<scalar> >& particleWeights_,
         labelListList& neighboringCellIDs_,
	 int& maxCellPerPart_	 
    );
    	
    //- Write results into the file
    virtual void write(word&);
    
    //- Write results into the file
    virtual void writeUserDefine
    (
    	 word&, 
	 word&,
         label& userLoopVar
     );

    //- Calculate conditionalAve for volScalarField
    virtual void condAve
    (   
	 const fvMesh&         	mesh_,
	 const volScalarField& 	var_,
	 const volScalarField& 	condVar_,
	 const scalar&         	minCondVar_,
	 const scalar&         	maxCondVar_,
	 const label&          	nBin_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
    );	

    //- Calculate conditionalAve for user define scalar field
    virtual void condAveUserDefine
    (   
	 const scalarField& 	var_,
	 const scalarField& 	condVar_,
	 const scalar&         	minCondVar_,
	 const scalar&         	maxCondVar_,
	 const label&          	nBin_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
    );

    //- Output subroutine
    virtual void writeOutput
    (   
	 const fvMesh&       mesh_,
	 fileName&           outputPath_,
	 word&               varName_,
	 word&               condVarName_,
	 const scalar&       minCondVar_,
	 const scalar&       maxCondVar_,
	 const label&        nBin_,
	 scalarField&        condAveVar_,
	 scalarField&        nReal_,
	 word&		     userDefineName_	
    );                
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
