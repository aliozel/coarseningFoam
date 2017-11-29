/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Calculate convolution

-----------------------------------------------------------------------------*/

#include "convolution.H"

namespace Foam {
    defineTypeNameAndDebug(convolution,0);
}

Foam::convolution::convolution
(
    const objectRegistry& obr,
    const constructFilter& box,
    const bool& readFilteredVars
)
:
    	obr_(obr),
	filter_(box),
	readFromFile_(false),
	readFilteredVars_(readFilteredVars)
{
	// Init
	init();	
}

Foam::convolution::~convolution()
{}

void Foam::convolution::init()
{
	// Read filtered variable from files
        if(readFilteredVars_)  readFromFile_ = true; 
	if(Pstream::parRun()) readFromFile_ = true;
	if(readFromFile_) Info << "Filtered variables is reading from files " <<endl;
}

Foam::volScalarField Foam::convolution::coarseningSca
(
	const word& name,
	const int& filterWidth,
	const volScalarField& phi
)
{

	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	char charfPhi[100];
	sprintf(charfPhi, "%dX%dX%d",filterWidth,filterWidth,filterWidth);
	
	volScalarField fphi
	(
	    IOobject
	    (
		name+charfPhi,
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ, //READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    phi
        );
			
	if(!readFromFile_) 
	{
		Info << tab << "Filtering " << name << " filter size = " << filterWidth << "X" << filterWidth << "X" << filterWidth << endl;
		forAll(mesh.cells(),cellI)
		{
			scalar totalVol(0);
			scalar fVar(0);
			//const labelList fCell(filter_.stencils(filterWidth,cellI)); 
			scalarField weights(filterWidth*filterWidth*filterWidth,1.);
			const labelList fCell(filter_.stencils(filterWidth,cellI,weights));
			
			//Info << "cellI = " << cellI << " fCell = " << fCell << endl;

			forAll(fCell,filterCellI)
			{
				label cID = fCell[filterCellI]; 
			
				scalar cWeight = weights[filterCellI];			

				totalVol +=   cWeight
					    * mesh.V()[cID];

				fVar +=    cWeight
					*  mesh.V()[cID]
					*      phi[cID]; 

			}
			
			if ( totalVol > 0 )
			{
				fVar /= totalVol;
			}
			else 
			{
				fVar = scalar(0);
			}
			
			fphi[cellI] = fVar;

		} 
		
		fphi.write();
		Info << tab << "Writing " << name+charfPhi << endl;		
		
	}else
	{
		Info << tab << "Reading " << name+charfPhi << endl;
		//fphi.read();
		IOobject fphimHeader
		(
        		name+charfPhi,
        		mesh.time().timeName(),
        		mesh,
        		IOobject::MUST_READ
		);
		volScalarField fphim(fphimHeader,mesh);
		fphi = fphim;
	}

	// For conditional averaging
	volScalarField fphiCond
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ, //READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fphi
        );
	fphiCond.write();

	return fphi;		
}
	
Foam::volScalarField Foam::convolution::coarseningScaFavre
(
	const word& name,
	const int& filterWidth,
	const volScalarField& voidfraction,  	
	const volScalarField& phi
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	char charfPhi[100];
	sprintf(charfPhi, "%dX%dX%d",filterWidth,filterWidth,filterWidth);
	
	volScalarField fphi
	(
	    IOobject
	    (
		name+charfPhi,
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ, //READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    phi
        );
		
	if(!readFromFile_) 
	{
		Info << tab << "Filtering " << name << " filter size = " << filterWidth << "X" << filterWidth << "X" << filterWidth << endl;

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction(0);
			scalar fVar(0);
			//const labelList fCell(filter_.stencils(filterWidth,cellI)); 
                        scalarField weights(filterWidth*filterWidth*filterWidth,1.);
                        const labelList fCell(filter_.stencils(filterWidth,cellI,weights));

			forAll(fCell,filterCellI)
			{
				label cID = fCell[filterCellI]; 

                                scalar cWeight = weights[filterCellI];

				barvoidfraction +=   cWeight
						   * voidfraction[cID] 
						   *     mesh.V()[cID];

				fVar +=   cWeight
					* voidfraction[cID] 
					*     mesh.V()[cID]
					*          phi[cID]; 

			}
			
			if ( barvoidfraction > 0 )
			{
				fVar /= barvoidfraction;
			}
			else 
			{
				fVar = scalar(0);
			}
			
			fphi[cellI] = fVar;
		} 

		Info << tab << "Writing " << name+charfPhi << endl;		
		fphi.write();
		
	}else
	{
		Info << tab << "Reading " << name+charfPhi << endl;
		//fphi.read();
		IOobject fphimHeader
		(
        		name+charfPhi,
        		mesh.time().timeName(),
        		mesh,
        		IOobject::MUST_READ
		);
		volScalarField fphim(fphimHeader,mesh);
		fphi = fphim;
	}

        // For conditional averaging
        volScalarField fphiCond
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ, //READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fphi
        );
        fphiCond.write();

	return fphi;		
}	

Foam::volVectorField Foam::convolution::coarseningVec
(
	const word& name,
	const int& filterWidth,
	const volVectorField& phi
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	char charfPhi[100];
	sprintf(charfPhi, "%dX%dX%d",filterWidth,filterWidth,filterWidth);
	
	volVectorField fphi
	(
	    IOobject
	    (
		name+charfPhi,
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ, //READ_IF_PRESENT,
		IOobject::AUTO_WRITE		
	    ),
	    phi
        );
	
	//Info << " Dimensions coarsening Vec" << phi.dimensions() << endl;
		
	if(!readFromFile_) 
	{
		Info << tab << "Filtering " << name << " filter size = " << filterWidth << "X" << filterWidth << "X" << filterWidth << endl;
		
		forAll(mesh.cells(),cellI)
		{
			scalar totalVol(0);
			vector fVar(0,0,0);
			//const labelList fCell(filter_.stencils(filterWidth,cellI)); 
                        scalarField weights(filterWidth*filterWidth*filterWidth,1.);
                        const labelList fCell(filter_.stencils(filterWidth,cellI,weights));

			forAll(fCell,filterCellI)
			{
				label cID = fCell[filterCellI]; 

                                scalar cWeight = weights[filterCellI];

				totalVol +=   cWeight
					    * mesh.V()[cID];

				fVar +=   cWeight
					* mesh.V()[cID]
					*      phi[cID]; 

			}
			
			if ( totalVol > 0 )
			{
				fVar /= totalVol;
			}
			else 
			{
				fVar = vector(0,0,0);
			}
			
			fphi[cellI].x() = fVar[0];
			fphi[cellI].y() = fVar[1];
			fphi[cellI].z() = fVar[2];
		} 

		Info << tab << "Writing " << name+charfPhi << endl;		
		fphi.write();
		
	}else
	{
		Info << tab << "Reading " << name+charfPhi << endl;
		//fphi.read();
		IOobject fphimHeader
		(
        		name+charfPhi,
        		mesh.time().timeName(),
        		mesh,
        		IOobject::MUST_READ
		);
		volVectorField fphim(fphimHeader,mesh);
		fphi = fphim;
	}

        // For conditional averaging
        volVectorField fphiCond
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ, //READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fphi
        );
        fphiCond.write();

	return fphi;		
}
	
Foam::volVectorField Foam::convolution::coarseningVecFavre
(
	const word& name,
	const int& filterWidth,
	const volScalarField& voidfraction,  	
	const volVectorField& phi
)
{
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	char charfPhi[100];
	sprintf(charfPhi, "%dX%dX%d",filterWidth,filterWidth,filterWidth);
	
	volVectorField fphi
	(
	    IOobject
	    (
		name+charfPhi,
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ, //READ_IF_PRESENT,
		IOobject::AUTO_WRITE		
	    ),
	    phi
        );
		
	if(!readFromFile_)
	{
		Info << tab << "Filtering " << name << " filter size = " << filterWidth << "X" << filterWidth << "X" << filterWidth << endl;

		forAll(mesh.cells(),cellI)
		{
			scalar barvoidfraction(0);
			vector fVar(0,0,0);
			//const labelList fCell(filter_.stencils(filterWidth,cellI)); 
			scalarField weights(filterWidth*filterWidth*filterWidth,1.);
                        const labelList fCell(filter_.stencils(filterWidth,cellI,weights));

			forAll(fCell,filterCellI)
			{
				label cID = fCell[filterCellI]; 

                                scalar cWeight = weights[filterCellI];

				barvoidfraction +=   cWeight 
						   * voidfraction[cID] 
						   *     mesh.V()[cID];

				fVar +=   cWeight
					* voidfraction[cID] 
					*     mesh.V()[cID]
					*          phi[cID]; 

			}
			
			if ( barvoidfraction > 0 )
			{
				fVar /= barvoidfraction;
			}
			else 
			{
				fVar = vector(0,0,0);
			}
			
			fphi[cellI].x() = fVar[0];
			fphi[cellI].y() = fVar[1];
			fphi[cellI].z() = fVar[2];
		} 

		Info << tab << "Writing " << name+charfPhi << endl;		
		fphi.write();

	}else
	{
		Info << tab << "Reading " << name+charfPhi << endl;
		//fphi.read();
		IOobject fphimHeader
		(
        		name+charfPhi,
        		mesh.time().timeName(),
        		mesh,
        		IOobject::MUST_READ
		);
		volVectorField fphim(fphimHeader,mesh);
		fphi = fphim;
	}

        // For conditional averaging
        volVectorField fphiCond
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ, //READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fphi
        );
        fphiCond.write();

	return fphi;		
}	

Foam::volScalarField Foam::convolution::coarseningSum
(
	const word& name,
	const int& filterWidth,
	const volScalarField& phi
)
{

	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	char charfPhi[100];
	sprintf(charfPhi, "%dX%dX%d",filterWidth,filterWidth,filterWidth);
	
	volScalarField fphi
	(
	    IOobject
	    (
		name+charfPhi,
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ, //READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    phi
        );
			
	forAll(mesh.cells(),cellI)
	{
		scalar fVar(0);
		//const labelList fCell(filter_.stencils(filterWidth,cellI)); 
                scalarField weights(filterWidth*filterWidth*filterWidth,1.);
                const labelList fCell(filter_.stencils(filterWidth,cellI,weights));

		scalar sumWeight(0);
		forAll(fCell,filterCellI)
		{
			label cID = fCell[filterCellI]; 
                        scalar cWeight = weights[filterCellI];
			sumWeight += cWeight;
			fVar += phi[cID]; 
		}
		fphi[cellI] = fVar/sumWeight;		
	} 

	fphi.write();
	Info << tab << "Writing " << name+charfPhi << endl;		

	return fphi;		
}

void Foam::convolution::write()
{
    if(readFromFile_) 
    {
	//Info << " Filtered variable is reading from file " << endl;
    }

}


// ************************************************************************* //
