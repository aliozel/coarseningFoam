/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    multipleVarsConditionalAve

Description
    Function object that calculate multipleVarsConditionalAve

SourceFiles
    multipleVarsConditionalAve.C

-----------------------------------------------------------------------------*/

#ifndef multipleVarsConditionalAve_H
#define multipleVarsConditionalAve_H

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

class multipleVarsConditionalAve
{
    //- Disallow default bitwise copy construct
    multipleVarsConditionalAve(const multipleVarsConditionalAve&);

    //- Disallow default bitwise assignment
    void operator=(const multipleVarsConditionalAve&);

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
    dictionary multipleVarsConditionalAveDict_;
    fileName outputPath_;   
    
    //- volScalarField variables 
    mutable wordList variableList;
    mutable wordList averagingVariableList;
    mutable labelList nBin;
    mutable scalarField minCondVar;
    mutable scalarField maxCondVar; 
    mutable boolList logBin;
    mutable Field <Field <scalar> > nReal;       
    mutable Field <Field <scalar> > condAveVar;
 	    
    //- Read from file
    bool readFromFile_;
    
    //- Boolean
    bool multConditionalAveraging_;

    //- Init
    void init();

    //- Init scalarField
    void initUserDefine();
    
    //- Calculate multipleVarsConditionalAve for volScalarField
    virtual void condAve
    (   
	 const fvMesh&         	mesh_,
	 const volScalarField& 	var_,
	 const volScalarField& 	condVarI_,
	 const scalar&     	minCondVarI_,
	 const scalar&     	maxCondVarI_,
	 const label&       	nBinI_,
	 const volScalarField& 	condVarJ_,
	 const scalar&     	minCondVarJ_,
	 const scalar&     	maxCondVarJ_,
	 const label&       	nBinJ_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
    );	

    //- Calculate multipleVarsConditionalAve for user define scalar field
    virtual void condAveUserDefine
    (   
	 const scalarField& 	var_,
	 const scalarField& 	condVarI_,
	 const scalar&     	minCondVarI_,
	 const scalar&     	maxCondVarI_,
	 const label&       	nBinI_,
	 const scalarField& 	condVarJ_,
	 const scalar&     	minCondVarJ_,
	 const scalar&     	maxCondVarJ_,
	 const label&       	nBinJ_,
	 scalarField&    	condAveVar_,
	 scalarField&    	nReal_
    );

    //- Output subroutine
    virtual void writeOutput
    (   
	const fvMesh&       mesh_,
	fileName&           outputPath_,
	word&               varName_,
	word&               condVarNameI_,
	const scalar&       minCondVarI_,
	const scalar&       maxCondVarI_,
	const label&        nBinI_,
	word&               condVarNameJ_,
	const scalar&       minCondVarJ_,
	const scalar&       maxCondVarJ_,
	const label&        nBinJ_,
	scalarField&        condAveVar_,
	scalarField&        nReal_,
	word&		    userDefineName_
    );    

public:

    //- Runtime type information
    TypeName("multipleVarsConditionalAve");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    multipleVarsConditionalAve
    (
	 const dictionary& dict,
	 const dictionary& subDict,
	 const objectRegistry& obr,
	 const label& nVariable,
	 const label& nAveragingVariable,
	 const label& nTotalCase,
	 const bool& multConditionalAveraging  
    );

    //- Destructor
    virtual ~multipleVarsConditionalAve();

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
    	
    //- Write results into the file
    virtual void write(word&);
    
    //- Write results into the file
    virtual void writeUserDefine
    (
    	 word&, 
	 word&,
         label& userLoopVar
     );
            
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
