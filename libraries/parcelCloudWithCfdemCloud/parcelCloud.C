/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Create parcel 

-----------------------------------------------------------------------------*/

#include "parcelCloud.H"

namespace Foam {
    defineTypeNameAndDebug(parcelCloud,0);
}

Foam::parcelCloud::parcelCloud
(
    const dictionary& dict,
    const objectRegistry& obr,
    cfdemCloud& sm,
    const labelList& particlesInSubDomain,
    const int& nPInParcel
)
:
	dict_(dict),
	obr_(obr),
	sm_(sm),
	particlesInSubDomain_(particlesInSubDomain),
	nPInParcel_(nPInParcel),
	particlesIn_(particlesInSubDomain.size(),labelList(nPInParcel,-1)),
	velocities_(particlesInSubDomain.size(),vector(0,0,0)),
	positions_(particlesInSubDomain.size(),vector(0,0,0)),
	radius_(particlesInSubDomain.size(),scalar(0)),
	types_(particlesInSubDomain.size(),-1),
	cellIDs_(particlesInSubDomain.size(),-1),
	readFromFile_(false),
	useAllParticles_(false),
	nPParcel_(0),
	locateModel_
    	(
          locateModelPost::New
          (
            dict_,
            obr_
          ) 
        ),
        turbulenceModelType_(dict.lookup("turbulenceModelType")),
        turbulence_
	(
            #if defined(version21) || defined(version16ext)
        	#ifdef comp
                    refCast<const fvMesh>(obr).lookupObject<compressible::turbulenceModel>
        	#else
                    refCast<const fvMesh>(obr).lookupObject<incompressible::turbulenceModel>
        	#endif
            #elif defined(version15)
        	refCast<const fvMesh>(obr).lookupObject<incompressible::RASModel>
            #endif
            (
        	turbulenceModelType_
            )
	)
{
	// Init
	init();
}

Foam::parcelCloud::~parcelCloud()
{}

void Foam::parcelCloud::init()
{
	// Use all particles to create parcel or pick some		
        if(dict_.found("useAllParticles"))
        {
                Pout << tab << "Using all particles to construct parcels " << endl;
                useAllParticles_ = true;
        }
	
}

Foam::labelListList Foam::parcelCloud::findParticlesIn()
{
	if(useAllParticles_)
	{
		particlesIn_.resize(particlesInSubDomain_.size(),labelList(nPInParcel_,-1));
	}else
	{
		particlesIn_.resize(particlesInSubDomain_.size()/nPInParcel_,labelList(nPInParcel_,-1));
	}
	
	// Number of particles in a parcel
	int k = nPInParcel_; 

	// Dimensions, exact OR approximate  
	int dim =3; double eps = 0;

	// Number of points
	int nPts;
	nPts = particlesInSubDomain_.size();

	Pout << tab << "Number of particles in a parcel = " << nPInParcel_ << endl;		

	// Domain Min/Max
	const fvMesh& mesh = refCast<const fvMesh>(obr_);
	const pointField& pp = mesh.points();	
	
	// Min, max x-coordinates	
	scalar minX = Foam::min(pp & vector(1,0,0));
	scalar maxX = Foam::max(pp & vector(1,0,0));

	// Min, max y-coordinates		
	scalar minY = Foam::min(pp & vector(0,1,0));
	scalar maxY = Foam::max(pp & vector(0,1,0));

	// Min, max z-coordinates		
	scalar minZ = Foam::min(pp & vector(0,0,1));
	scalar maxZ = Foam::max(pp & vector(0,0,1));

	// Squared radius
	const scalar sqRad = pow( (maxX - minX) * (maxX - minX) 
				 +(maxY - minY) * (maxY - minY)
				 +(maxZ - minZ) * (maxZ - minZ), 0.5);

	Pout << tab << "Squared radius = " << sqRad << endl;

	// Data points
	ANNpointArray dataPts;
	// Query points
	ANNpoint queryPt;

	ANNidxArray	nnIdx;          // 	near neighbour indices
	ANNdistArray 	dists;		//	near neighbour distances
	ANNkd_tree* 	kdTree;		//	search structure

	Pout << tab << "Created kdTree variables " << endl;

	// Allocate 
	queryPt = annAllocPt(dim);
	dataPts = annAllocPts(nPts, dim);
	nnIdx = new ANNidx[k];
	dists = new ANNdist[k];

	Pout << tab << "Allocated kdTree variables " << endl;

        labelList particleCreateParcelList(particlesInSubDomain_.size());

        for(int ii =0; ii < particlesInSubDomain_.size(); ii++)
	{
			label particleGlobalID = particlesInSubDomain_[ii];
			dataPts[ii][0] = sm_.position(particleGlobalID).x();
			dataPts[ii][1] = sm_.position(particleGlobalID).y();
			dataPts[ii][2] = sm_.position(particleGlobalID).z();
                	particleCreateParcelList[ii] = particleGlobalID;
	}

        Pout << tab << "Creating kdTree..." << endl;
        kdTree = new ANNkd_tree(dataPts, nPts, dim);

	Pout << tab << "Entering particle loops to create parcel" << endl; 

    	label parcelI = 0;
    	forAll(particlesInSubDomain_,ii)
	{				
	    	label particleGlobalID = particlesInSubDomain_[ii];

        	if ( particleCreateParcelList[ii] > -1 )
        	{
        	    queryPt[0] = sm_.position(particleGlobalID).x();
        	    queryPt[1] = sm_.position(particleGlobalID).y();
        	    queryPt[2] = sm_.position(particleGlobalID).z();

        	    kdTree->annkFRSearch(
                                	    queryPt,				// query point					
                                	    sqRad,				// squared radius
                                	    k,                  		// number of the near neighbours to return
                                	    nnIdx,				// nearest neighbor array
                                	    dists,				// dist to near neighbours
                                	    eps			);

        	    int partSum = 0;
		    int i = 0;
			
		    while( i < k && partSum < nPInParcel_ )
		    {
						
			if ( particleCreateParcelList[nnIdx[i]] != -1 )
                	{
			    //Info  << " parcelI " << parcelI << " partSum " << partSum << " Neighbour part " 
			    //      << nnIdx[i] << " particlesInSubDomain " << particlesInSubDomain_[nnIdx[i]] 
			    //      << " particlesIn " << particlesIn_[parcelI][partSum] << endl;
			    
                	    if (!useAllParticles_) particleCreateParcelList[nnIdx[i]] = -1 ;
                	    particlesIn_[parcelI][partSum] = particlesInSubDomain_[nnIdx[i]];
                	    partSum++;
			    if( partSum == nPInParcel_ ) parcelI++;
                	};
			i++;				
        	    };

        	}

    	}
	
	// Resize number of parcel 
	nPParcel_ = parcelI;
	if (!useAllParticles_)
	{
	  particlesIn_.resize(nPParcel_);
	  velocities_.resize(nPParcel_);
	  positions_.resize(nPParcel_);
	  radius_.resize(nPParcel_);
	  types_.resize(nPParcel_);
	  cellIDs_.resize(nPParcel_);
	}  
		
        int totalNpParcel = parcelI;
        if(Pstream::parRun())
        {
            reduce(totalNpParcel, sumOp<int>());
        }

        Info << tab << "Total number of parcels = " << totalNpParcel << endl;
		
	return particlesIn_;
	
}


void Foam::parcelCloud::createParcel()
{
    if(!readFromFile_) 
    {
	// Find closest particles
	findParticlesIn();
	// Parcel radius
	parcelRadius();
	// Calculate parcel velocity
	parcelVel();
	// Locate parcel
	parcelPos();	
	// Find parcel cell
        locateM().findCell(NULL,positions_,cellIDs_,particlesInSubDomain_.size(),radius_);
    }else
    {
    	// Do nothing
    }	
}

vectorField Foam::parcelCloud::parcelVel()
{
    if(!readFromFile_) 
    {
	for(label parI=0; parI < particlesIn_.size(); parI++)
	{
		for(label partIn=0; partIn < nPInParcel_; partIn++)
		{
			velocities_[parI] += sm_.velocity(particlesIn_[parI][partIn]);
		}
		velocities_[parI] /= nPInParcel_;
	}	
    }else
    {
    	// Do nothing
    }
    return velocities_;
}

vectorField Foam::parcelCloud::parcelPos()
{
    if(!readFromFile_) 
    {
	for(label parI=0; parI < particlesIn_.size(); parI++)
	{
		for(label partIn=0; partIn < nPInParcel_; partIn++)
		{
			positions_[parI] += sm_.position(particlesIn_[parI][partIn]);
		}
		positions_[parI] /= nPInParcel_;
	}	
    }else
    {
    	// Do nothing
    }
    return positions_;
}

scalarField Foam::parcelCloud::parcelRadius()
{
    if(!readFromFile_) 
    {
	    for(int index = 0;index < nPParcel_; ++index)
	    {
        	radius_[index] = sm_.radius(index);
	    }
    }else
    {
     	// Do nothing
    }		    
    return radius_;
}

void Foam::parcelCloud::read()
{
    if(readFromFile_) 
    {
	// Do nothing
    }

}

void Foam::parcelCloud::write()
{
    if(!readFromFile_) 
    {
	// Do nothing
    }

}

// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

int Foam::parcelCloud::numberOfParticles()
{
    return nPParcel_;		
}

vector Foam::parcelCloud::velocity(int& index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = velocities_[index][i];
    return vel;
}

vector Foam::parcelCloud::position(int& index)
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positions_[index][i];
    return pos;
}

scalar Foam::parcelCloud::radius(int& index)
{
    scalar rad;
    rad = radius_[index];
    return rad; 	
}

label Foam::parcelCloud::type(int& index)
{
    label ty = static_cast<label>(types_[index]);
    return ty; 	
}

label Foam::parcelCloud::cellIDs(int& index)
{
    label cellI = cellIDs_[index];
    return cellI;
}

labelList Foam::parcelCloud::particlesInParcel(int& index)
{
    labelList partsIn = particlesIn_[index]; 
    return partsIn;	
}

// ************************************************************************* //
