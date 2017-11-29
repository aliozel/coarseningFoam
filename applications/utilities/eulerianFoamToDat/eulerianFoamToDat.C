/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHAN" " ILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    calcVoFMassFlowRate

Description
    Calculation of mass flow rate of VoF results

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"

#include "timeSelector.H"
#include "fvCFD.H"

#include "IOmanip.H"
#include "OFstream.H"

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    // Create output folder	
    OFstream* outputFile;
         	
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
	
	#include "readEulerianVar.H"

	outputFile =  new OFstream(mesh.time().path()/"dumpEulerian"+runTime.timeName()+".dat");
        Info << tab << "Writing Eulerian data in dump format " << mesh.time().path()/"dumpEulerian"+runTime.timeName()+".dat" << endl;
	
	//- Write number of cells
	*outputFile  << mesh.nCells() << endl;

        forAll(mesh.cells(),cellI)
	{
	   *outputFile << mesh.C()[cellI][0] << " "  
	               << mesh.C()[cellI][1] << " " 
		       << mesh.C()[cellI][2] << " "   
		       <<        U[cellI][0] << " " 
		       <<        U[cellI][1] << " " 
		       <<        U[cellI][2] << " " 
		       <<       Us[cellI][0] << " " 
		       <<       Us[cellI][1] << " " 
		       <<       Us[cellI][2] << " " 
		       << voidfraction[cellI] << endl;

	}  

    }
   
    delete outputFile;      	

}


// ************************************************************************* //
