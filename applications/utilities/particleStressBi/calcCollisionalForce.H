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

InNamespace
    Foam

Description
    Function prototype for all simple post-processing functions
    e.g. filterinEulerianFieldDivPhi, filterinEulerianFieldMagU etc.

\*---------------------------------------------------------------------------*/

#ifndef calcCollisionalForce_H
#define calcCollisionalForce_H

#include "particleStressBi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{	

	void calcLiqIn(scalarField& liq, int& i, int& j, vectorField& liquidVol);
	
	void calcLiqBack(scalarField& liq, int& i, int& j, vectorField& liquidVol);

	void calcForce(	int& i, vector& velJ, vector& posJ, vector& omeJ, scalar& radJ, int& typeJ, 
			int& j, vector& velI, vector& posI, vector& omeI, scalar& radI, int& typeI, 	
			word& collisionModelI,
			scalarField& youngsModulus,
			scalarField& poissonsRatio,
			vectorField& coefficientRestitution,
			vectorField& coefficientFriction,
			scalar& k_n, 
			scalar& k_t, 
			scalar& gamma_n, 
			scalar& gamma_t,
			scalar& e_n,		  
			scalar& mu_f,
			vectorField& fcoll, 
			vectorField& ftan,
			const scalar& dt, 
			bool& tangential_history, 
			scalarField& liq, 
			bool& liquid_transfer, 
			vectorField& liquidVol, 
			scalar& surf_tension, 
			scalar& fluid_visc, 
			vectorField& fcap, 
			vectorField& fvisc,
			bool& first_touch,
			bool& cohesion,
			scalar& minimumDistanceVdW,
			vectorField& cohEnergyDens,
			vectorField& fcoh,
			const scalar& rhop,	   	    					        
			const bool& bidisperse,
			symmTensor& sigma_coll_JI,
			const bool& verbose,
			scalar& numberOfPairs		);

			
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
