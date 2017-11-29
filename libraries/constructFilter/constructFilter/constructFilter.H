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
    Foam::constructFilter

Description
    Abstract class for LES filters

SourceFiles
    constructFilter.C
    newFilter.C

\*---------------------------------------------------------------------------*/

#ifndef constructFilter_H
#define constructFilter_H

#include "volFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class constructFilter Declaration
\*---------------------------------------------------------------------------*/

class constructFilter
{
    // Private data

        const fvMesh& mesh_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        constructFilter(const constructFilter&);
        void operator=(const constructFilter&);


public:

    //- Runtime type information
    TypeName("constructFilter");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            constructFilter,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& constructFilterDict
            ),
            (mesh, constructFilterDict)
        );


    // Constructors

        //- Construct from components
        constructFilter(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a reference to the selected LES filter
        static autoPtr<constructFilter> New
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~constructFilter()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

    	// Read from file
    	virtual bool readFromFile() = 0;

    	virtual labelList stencils
    	(
        	const int& filterWidth,
        	const label& cellID
    	) const = 0;

        virtual labelList stencils
        (
                const int& filterWidth,
                const label& cellID,
                scalarField& weights
        ) const = 0;

        //- Read stencils from file
        virtual void read() const = 0;

        //- Write stencils into the file
        virtual void write() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
