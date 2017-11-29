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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
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
    argList::validArgs.append("direction");
    argList::validArgs.append("nBin");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    const word direction = args[1];
    const word charNbin = args[2];
    
    int iStream(-1);
    int iSpan1(-1);
    int iSpan2(-1);
    if(direction=="x")
    {
        iStream = 0;
	iSpan1  = 1;
	iSpan2  = 2;
    }
    else if(direction=="y")
    {
        iStream = 1;
	iSpan1  = 0;
	iSpan2  = 2;
    }
    else if(direction=="z")
    {
        iStream = 2;
	iSpan1  = 0;
	iSpan2  = 1;
    }
    else
    {
        FatalError << " Direction is not defined " << abort(FatalError);
    }
    
    Info << "Streamwise direction: " << direction << endl;

    char* pEnd;
    int nBinSpan = strtod(charNbin.c_str(), &pEnd);
    int nBinStream;

    //pimpleControl pimple(mesh);
    const pointField& pp = mesh.points();

    // Min, max x-coordinates       
    scalar xMin = min(pp & vector(1,0,0));
    scalar xMax = max(pp & vector(1,0,0));

    // Min, max y-coordinates               
    scalar yMin = min(pp & vector(0,1,0));
    scalar yMax = max(pp & vector(0,1,0));

    // Min, max z-coordinates               
    scalar zMin = min(pp & vector(0,0,1));
    scalar zMax = max(pp & vector(0,0,1));

    Info << "Domain " << " x=[" << xMin << ":" << xMax << "]" << " y=[" << yMin << ":" << yMax << "]" << " z=[" << zMin << ":" << zMax << "]" << endl;
    
    scalar minCoordSpan(0);
    scalar maxCoordSpan(0);
    scalar minCoordStream(0);
    scalar maxCoordStream(0);
    scalar binSizeSpan(0);
    scalar binSizeStream(0);    
        
    if(iStream == 0)
    {
	if( ((zMax-zMin)-(yMax-yMin)) < 0 ) FatalError << " Domain is not squared section " << abort(FatalError);
	minCoordSpan = yMin;
	maxCoordSpan = yMax;
	minCoordStream = xMin;
	maxCoordStream = xMax;
    }
    else if(iStream == 1)
    {
	if( ((zMax-zMin)-(xMax-xMin)) < 0 ) FatalError << " Domain is not squared section " << abort(FatalError);
	minCoordSpan = zMin;
	maxCoordSpan = zMax;
	minCoordStream = yMin;
	maxCoordStream = yMax;
    }
    else if(iStream == 2)
    {
	if( ((yMax-yMin)-(xMax-xMin)) < 0 ) FatalError << " Domain is not squared section " << abort(FatalError);
	minCoordSpan = xMin;
	maxCoordSpan = xMax;
	minCoordStream = zMin;
	maxCoordStream = zMax;
    }
    else
    {
        FatalError << " Direction is not defined " << abort(FatalError);
    }
    
    nBinStream = nBinSpan * ( maxCoordStream-minCoordStream ) / (maxCoordSpan-minCoordSpan);
    
    Info << "Number of bins for spanwise direction: " << nBinSpan << endl;
    Info << "Number of bins for streamwise direction: " << nBinStream << endl;
    
    binSizeStream = ( maxCoordStream-minCoordStream )/nBinStream;   
    binSizeSpan = ( maxCoordSpan-minCoordSpan )/nBinSpan; 
    
    scalarField rdfStream(nBinStream,0);
    scalarField rdfSpan(nBinSpan,0);  
    
    scalarField nRealStream(nBinStream,0);
    scalarField nRealSpan(nBinSpan,0);  
           
    vector coordI(0,0,0);
    vector coordJ(0,0,0);    
    scalar alphaFlucI(0);
    scalar alphaFlucJ(0);
    scalar LenStream(0);
    scalar LenSpan(0);
    scalar LenSpan1(0);
    scalar LenSpan2(0);
        
    scalar alphaFluc2(0);
    scalar NrealAlphaFluc2(0);
    vector crossIJ(0,0,0);
    
    scalar alphaMean(0);
    scalar alphaalpha(0);
    scalarField alphaIalphaJStream(nBinStream,0);
    scalarField alphaIalphaJSpan(nBinSpan,0);  
    
    scalar SMALL(1.e-12);
     	
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
	
	mesh.readUpdate();
	#include "readEulerianVar.H"

        forAll(mesh.cells(),cellI)
	{
          coordI = mesh.C()[cellI];
	  alphaFlucI = alpha[cellI] - alphaMean;
	  alphaFluc2 += ( alpha[cellI] - alphaMean ) * ( alpha[cellI] - alphaMean );
	  alphaalpha += alpha[cellI]*alpha[cellI]; 
	  NrealAlphaFluc2 += 1.; 
	  
	  forAll(mesh.cells(),cellJ)
	  {
	    coordJ = mesh.C()[cellJ];
	    alphaFlucJ = alpha[cellJ] - alphaMean;

	    LenStream = Foam::sqrt( ( coordI[iStream]-coordJ[iStream] )
		                   *( coordI[iStream]-coordJ[iStream] ) );	     
	    LenSpan1  = Foam::sqrt( ( coordI[iSpan1] -coordJ[iSpan1]  )
		                   *( coordI[iSpan1] -coordJ[iSpan1]  ) );	
	    LenSpan2  = Foam::sqrt( ( coordI[iSpan2] -coordJ[iSpan2]  )
		                   *( coordI[iSpan2] -coordJ[iSpan2]  ) );
				      
	    // Streamwise direction
	    if( LenSpan1 <= SMALL && LenSpan2 <= SMALL && (LenStream-(maxCoordStream-minCoordStream)/2.) <= 0 )
	    {
	       // Conditionally averaging
	       label binSt                = floor(( LenStream - minCoordStream ) / binSizeStream );
	       rdfStream[binSt]          += alphaFlucI*alphaFlucJ;
	       alphaIalphaJStream[binSt] += alpha[cellI]*alpha[cellJ];
	       nRealStream[binSt]        += 1.;   		      
	    }
	    // Spanwise direction              
	    else if( (LenSpan1 <= SMALL && LenStream <= SMALL) || (LenSpan2 <= SMALL && LenStream <= SMALL) )
	    {
	          if( LenSpan1 <= SMALL ) LenSpan = LenSpan2; 
	          if( LenSpan2 <= SMALL ) LenSpan = LenSpan1; 
		  if( (LenSpan-(maxCoordSpan-minCoordSpan)/2.) <= 0 )
		  {		  
		    // Conditionally averaging
		    label binSp                = floor(( LenSpan  - minCoordSpan   ) / binSizeSpan   );
		    rdfSpan[binSp]            += alphaFlucI*alphaFlucJ;
    		    alphaIalphaJSpan[binSp]   += alpha[cellI]*alpha[cellJ];
		    nRealSpan[binSp]          += 1.; 
		  }  				   
	    }		 	      
	  }		
	}
    }

    // Create output folder	
    OFstream* outputFileStream;
    outputFileStream =  new OFstream(mesh.time().path()/"alpprime2rdfStream.dat");
    Info << tab << "Writing rdf along streamwise direction  into the file " << mesh.time().path()/"alpprime2rdfStream.dat" << endl;
    *outputFileStream  << "#nBin \t alpprime2rdf \t alp2rdf \t NrealStream \t alphaMean \t <alpprime^2> \t <alpp^2> \t nReal \t  "<< endl; 
    for(int ii=0; ii<nBinStream; ii++)
    {   
       *outputFileStream  << minCoordStream+(ii+0.5)*binSizeStream << tab 
	                  << rdfStream[ii]                         << tab
			  << alphaIalphaJStream[ii]		   << tab  
			  << nRealStream[ii]                       << tab
			  << alphaMean				   << tab			  
			  << alphaFluc2                            << tab
			  << alphaalpha                            << tab
			  << NrealAlphaFluc2                       << endl; 
    }
    
    OFstream* outputFileSpan;
    outputFileSpan =  new OFstream(mesh.time().path()/"alpprime2rdfSpan.dat");
    Info<< tab << "Writing rdf along spanwise direction into the file " << mesh.time().path()/"alpprime2rdfSpan.dat" << endl;
    *outputFileSpan  << "#nBin \t alpprime2rdf \t alp2rdf \t NrealSpan \t alphaMean \t <alpprime^2> \t <alpp^2> \t nReal \t  "<< endl;

    for(int ii=0; ii<nBinSpan; ii++)
    {   
       *outputFileSpan    << minCoordSpan+(ii+0.5)*binSizeSpan     << tab 
	                  << rdfSpan[ii]                           << tab
			  << alphaIalphaJSpan[ii]		   << tab  
			  << nRealSpan[ii]                         << tab
			  << alphaMean				   << tab
			  << alphaFluc2                            << tab
			  << alphaalpha                            << tab
			  << NrealAlphaFluc2                       << endl; 
    } 
        
    delete outputFileStream;      	
    delete outputFileSpan;      	

}


// ************************************************************************* //
