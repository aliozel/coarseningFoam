/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::GaussianFilter

Description
    Simple top-hat filter used in dynamic LES models.

    Implemented as a surface integral of the face interpolate of the field.

SourceFiles
    GaussianFilter.C

\*---------------------------------------------------------------------------*/

#ifndef GaussianFilter_H
#define GaussianFilter_H

#include "constructFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GaussianFilter Declaration
\*---------------------------------------------------------------------------*/

class GaussianFilter
:
    public constructFilter

{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        GaussianFilter(const GaussianFilter&);
        void operator=(const GaussianFilter&);

protected:

    // Private data

	// Read stencils from file
	bool readFromFile_;
	
	// Cell matrix
	mutable labelListList mA_;
	
	// Maximum number of cells
	mutable int mAx_;
	mutable int mAy_;
	mutable int mAz_;

	// Init
	void init();

public:

    //- Runtime type information
    TypeName("Gaussian");


    // Constructors

        //- Construct from components
        GaussianFilter(const fvMesh& mesh);

        //- Construct from IOdictionary
        GaussianFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~GaussianFilter()
    {}


    // Member Functions
    
	// Read from file
    	virtual bool readFromFile();

    	//- Call stencils
    	//virtual void stencils(const label&, const int&);
    	virtual labelList stencils
    	(
        	const int& filterWidth,
        	const label& cellID
    	) const;

        //- Call stencils with weights
        //virtual void stencils(const label&, const int&, scalarField& weights);
        virtual labelList stencils
        (
                const int& filterWidth,
                const label& cellID,
                scalarField& weights
        ) const;

    	//- Read stencils from file
    	virtual void read() const;

    	//- Write stencils into the file
    	virtual void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
