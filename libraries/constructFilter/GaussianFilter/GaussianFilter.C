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

\*---------------------------------------------------------------------------*/

#include "GaussianFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GaussianFilter, 0);
    addToRunTimeSelectionTable(constructFilter, GaussianFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussianFilter::GaussianFilter
(
    const fvMesh& mesh
)
:
    constructFilter(mesh)
{}


Foam::GaussianFilter::GaussianFilter(const fvMesh& mesh, const dictionary&)
:
    constructFilter(mesh),
    readFromFile_(false),
    mA_(0),
    mAx_(0),
    mAy_(0),
    mAz_(0)
{
    init();	
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GaussianFilter::init()
{
	const pointField& pp = mesh().points();

	// Min, max x-coordinates	
	scalar xMin = min(pp & vector(1,0,0));
	scalar xMax = max(pp & vector(1,0,0));

	// Min, max y-coordinates		
	scalar yMin = min(pp & vector(0,1,0));
	scalar yMax = max(pp & vector(0,1,0));

	// Min, max z-coordinates		
	scalar zMin = min(pp & vector(0,0,1));
	scalar zMax = max(pp & vector(0,0,1));	

	Info << "\nDomain " << " x=[" << xMin << ":" << xMax << "]" << " y=[" << yMin << ":" << yMax << "]" << " z=[" << zMin << ":" << zMax << "]" << endl;
		
	const scalar cell_volume_ref = mesh().V()[0];
	forAll(mesh().cells(),cellI)
	{
		if ( abs(cell_volume_ref-mesh().V()[cellI]) > SMALL)  FatalError<< "Non-uniform mesh !!! " << abort(FatalError);
	}

	const scalar delta = pow(cell_volume_ref,1./3.);
	Info << "Uniform mesh,  volume of a cell = " << cell_volume_ref << endl;
	Info << "Delta_x = Delta_y = Delta_z = " << delta << endl;
	
	// Create cell matrix
	labelList mACell(3,-1);
	mA_.resize(mesh().nCells(),mACell);
	int mAx(0); int mAy(0); int mAz(0);
	forAll(mesh().cells(),cellI)
	{
		
		//mAx = floor((mesh.C()[cellI][0]-xMin)/delta);
		//mAy = floor((mesh.C()[cellI][1]-yMin)/delta);
		//mAz = floor((mesh.C()[cellI][2]-zMin)/delta);			

		mAx = (mesh().C()[cellI][0]-xMin)/delta;
		mAy = (mesh().C()[cellI][1]-yMin)/delta;
		mAz = (mesh().C()[cellI][2]-zMin)/delta;	

		mA_[cellI][0] = mAx;
		mA_[cellI][1] = mAy;
		mA_[cellI][2] = mAz;

		//if ( mAx > mAx_ ) mAx_ = mAx; 
		//if ( mAy > mAy_ ) mAy_ = mAy; 
		//if ( mAz > mAz_ ) mAz_ = mAz;			

		mAx_ = mAx; 
		mAy_ = mAy; 
		mAz_ = mAz;
		
		//Info << " mA_  = " << mA_[cellI] << endl; 	
	} 
	
	Info << "Box filter; mAx = " << mAx_ << " mAy = " << mAx_ << " mAz = " << mAz_ << endl;
}

bool Foam::GaussianFilter::readFromFile()
{
     Info << "Filter stencils will be reading from the file " << endl;
     return	readFromFile_;
}

Foam::labelList Foam::GaussianFilter::stencils
(
    const int& filterWidth,
    const label& cellID
)const
{
    labelList stencilList(filterWidth*filterWidth*filterWidth,-1);
    
    //Info << " mAx_ " << mAx_ << " mAy_ " << mAx_ << " mAz_ " << mAz_ << endl;

    if(!readFromFile_) 
    {						
	 label nStencil = 0;
	 int filter = floor(filterWidth/2) - 1;
	 for ( label kk = -(filter+1); kk <= filter+1; kk++)
	 {
		 for ( label jj = -(filter+1); jj <= filter+1; jj++)
		 {
			 for ( label ii = -(filter+1); ii <= filter+1; ii++)
			 {									
				 stencilList[nStencil] =   ((mA_[cellID][0] + ii + mAx_ + 1) % (mAx_+1))
					         	 + ((mA_[cellID][1] + jj + mAy_ + 1) % (mAy_+1))*(mAx_+1)
						     	 + ((mA_[cellID][2] + kk + mAz_ + 1) % (mAz_+1))*(mAx_+1)*(mAy_+1) ;					

				 nStencil++;		
			 }	
		 }
	 }
    }else
    {
	    // Do nothing		
    }

    //Info << " stencilList = " << stencilList << endl;
    	
    return stencilList;		
}

Foam::labelList Foam::GaussianFilter::stencils
(   
    const int& filterWidth,
    const label& cellID,
    scalarField& weights
)const
{   
    labelList stencilList(filterWidth*filterWidth*filterWidth,-1);
    
    //Info << " mAx_ " << mAx_ << " mAy_ " << mAx_ << " mAz_ " << mAz_ << endl;
    
    if(!readFromFile_)
    {    
         label nStencil = 0;
         int filter = floor(filterWidth/2) - 1;
         for ( label kk = -(filter+1); kk <= filter+1; kk++)
         {       
                 for ( label jj = -(filter+1); jj <= filter+1; jj++)
                 {       
                         for ( label ii = -(filter+1); ii <= filter+1; ii++)
                         {       
                                 stencilList[nStencil] =   ((mA_[cellID][0] + ii + mAx_ + 1) % (mAx_+1))
                                                         + ((mA_[cellID][1] + jj + mAy_ + 1) % (mAy_+1))*(mAx_+1)
                                                         + ((mA_[cellID][2] + kk + mAz_ + 1) % (mAz_+1))*(mAx_+1)*(mAy_+1) ;

                    		 // Find distances between particle and neighbouring cells                                   
                    		 scalar dist_s =  mag( mesh().C()[cellID] - mesh().C()[stencilList[nStencil]])
                                                / Foam::pow(mesh().V()[stencilList[nStencil]],1./3.);
				 // Compute weights
                    		 if(dist_s <= 0.5)
                    		 {
                            		weights[nStencil] =  1./4.*pow(dist_s,4)-5./8.*pow(dist_s,2)+115./192.;
                    		 }
                    		 else if (dist_s > 0.5 && dist_s <= 1.5)
                    		 {
                            		weights[nStencil] = -1./6.*pow(dist_s,4)+5./6.*pow(dist_s,3)-5./4.*pow(dist_s,2)+5./24.*dist_s+55./96.;
                    		 }
                    		 else if (dist_s > 1.5 && dist_s <= 2.5)
                    		 {
                            		weights[nStencil] =  pow(2.5-dist_s,4)/24.;
                    		 }
                    		 else
                    		 {
                            	 	weights[nStencil] =  0;
                   		 }
                                
                                 nStencil++;
                         }
                 }
         }
    }else
    {       
            // Do nothing               
    }
    
    //Info << " stencilList = " << stencilList << endl;
    
    return stencilList;
}

void Foam::GaussianFilter::write() const
{
    if(!readFromFile_)
        {

    }

}

void Foam::GaussianFilter::read() const
{
    if(readFromFile_)
        {

    }

}
// ************************************************************************* //
