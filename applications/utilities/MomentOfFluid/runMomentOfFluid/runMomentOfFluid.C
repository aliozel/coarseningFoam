/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    testMomentOfFluid

Description
    Driver for testing MomentOfFluid intersection algorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvMesh.H"
#include "argList.H"
#include "volFields.H"
#include "fvCFD.H"

#include "MomentOfFluid.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    runTime.setTime(timeDirs[0], 0);

#   include "createMesh.H"

    // Create fields
    volScalarField alpha
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField interfaceNormals
    (
        IOobject
        (
            "interfaceNormals",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector(0,0,0))
    );
    
    volVectorField interfaceNormalsGradAlpha
    (
        IOobject
        (
            "interfaceNormalsGradAlpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector(0,0,0))
    );      

    volScalarField kappaR
    (
        IOobject
        (
            "kappaR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("zero", dimensionSet(0,-1,0,0,0), 0)
    );
    
    volScalarField kappaG
    (
        IOobject
        (
            "kappaG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("zero", dimensionSet(0,-1,0,0,0), 0)
    );

    // Cell centers
    volVectorField refCentres(mesh.C());
    
    // Construct intersector
    MomentOfFluid mof(mesh);

    // Compute surfaces and output
    mof.constructInterface
    (
        alpha.internalField(),
        refCentres.internalField()
    );

    // Surface normals by reconstruction
    vectorField nByMofvec = mof.interfaceNormals();

    // Surface normals by grad(alpha)/|grad(alpha)|
    volVectorField gradAlpha(fvc::grad(alpha));    
    forAll(mesh.cells(),cellI)
    {
    
      interfaceNormals[cellI] =  nByMofvec[cellI];

      scalar magGradAlpha = mag(gradAlpha[cellI]);
      if(magGradAlpha>0) 
      {
      	interfaceNormalsGradAlpha[cellI] = gradAlpha[cellI]/magGradAlpha;
      }else
      {
        interfaceNormalsGradAlpha[cellI] = vector(0,0,0);
      }

    }
    interfaceNormals.correctBoundaryConditions();
    interfaceNormalsGradAlpha.correctBoundaryConditions();
    
    kappaR = fvc::div(interfaceNormals);
    kappaG = fvc::div(interfaceNormalsGradAlpha);

    // Create file output file
    fileName dirName(mesh.time().path());
    // Open stream for output
    OFstream file(dirName/"pdfNByReconstGradAlpha");
    file<< "# deltaBin " << tab << "nByREconstruction " << tab << "nByGradAlpha " << endl;

    // Pdf-nbin
    label nBin(100);

    scalar minKappaR(min(kappaR.internalField()));
    scalar maxKappaR(max(kappaR.internalField()));

    scalar minKappaG(min(kappaG.internalField()));
    scalar maxKappaG(max(kappaG.internalField()));

    scalar deltaBinR = ( maxKappaR - minKappaR ) / nBin;
    scalar deltaBinG = ( maxKappaG - minKappaG ) / nBin;
        
    scalarField nRealR(nBin,0);
    scalarField nRealG(nBin,0);
    
    forAll(mesh.cells(),cellI)
    {    
      label iBinR = floor( ( kappaR[cellI] - minKappaR ) / deltaBinR ) ;
      if( iBinR > -1 && iBinR < nBin) nRealR[iBinR] +=1.; 
      
      label iBinG = floor( ( kappaG[cellI] - minKappaG ) / deltaBinG ) ;
      if( iBinG > -1 && iBinG < nBin) nRealG[iBinG] +=1.;
    }  
    
    for(label ii=0; ii < nBin; ii++)
    {
    	file << minKappaR + (ii+0.5)*deltaBinR << tab << nRealR[ii]/sum(nRealR) << tab
	     << minKappaG + (ii+0.5)*deltaBinG << tab << nRealG[ii]/sum(nRealG) << endl;
    } 
    
    // Output VTK file
    mof.outputSurface();

   // Output
    interfaceNormals.write(); 
    interfaceNormalsGradAlpha.write();           	

    kappaR.write(); 
    kappaG.write();
    
}
