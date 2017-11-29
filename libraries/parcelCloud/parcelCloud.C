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
	const int & partNP,
	List<vector>& partPositions,
	List<vector>& partVelocities,
	List<scalar>& partRadii,	
	const int& nPInParcel
)
:
	dict_(dict),
	obr_(obr),
	partNP_(partNP),
	partPositions_(partPositions),
	partVelocities_(partVelocities),
	partRadii_(partRadii),
	nPInParcel_(nPInParcel),
	particlesIn_(partNP_,labelList(nPInParcel,-1)),
	velocities_(partNP_,vector(0,0,0)),
	positions_(partNP_,vector(0,0,0)),
	radius_(partNP_,scalar(0)),
	types_(partNP_,-1),
	cellIDs_(partNP_,-1),
	useAllParticles_(false),
	nPParcel_(0),
	treeSearch_(dict_.lookup("treeSearch")),
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
	),
	parcelCoarseningDict_(dict_)
{
	// Init
	init();
}

Foam::parcelCloud::~parcelCloud()
{}

void Foam::parcelCloud::init()
{
	// Use all particles to create parcel or pick some		
        if(parcelCoarseningDict_.found("useAllParticles"))
        {
                Pout << tab << "Using all particles to construct parcels " << endl;
                useAllParticles_ = true;
        }	
}

Foam::labelListList Foam::parcelCloud::findParticlesIn()
{
	if(useAllParticles_)
	{
		particlesIn_.resize(partNP_,labelList(nPInParcel_,-1));
	}else
	{
		particlesIn_.resize(partNP_/nPInParcel_,labelList(nPInParcel_,-1));
	}
	
	// Number of particles in a parcel
	int k = nPInParcel_; 

	// Dimensions, exact OR approximate  
	int dim =3; double eps = 0;

	// Number of points
	int nPts;
	nPts = partNP_;

	Pout << tab << "Number of particles in a parcel = " << nPInParcel_ << endl;		

	/*
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
	//const scalar sqRad = pow( (maxX - minX) * (maxX - minX) 
	 scalar sqRad = pow( (maxX - minX) * (maxX - minX) 
			    +(maxY - minY) * (maxY - minY)
			    +(maxZ - minZ) * (maxZ - minZ), 0.5);
	
	*/
	scalar sqRad = 0; 
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

        labelList particleCreateParcelList(partNP_);

        for(int ii =0; ii < partNP_; ii++)
	{
			dataPts[ii][0] = partPositions_[ii][0];
			dataPts[ii][1] = partPositions_[ii][1];
			dataPts[ii][2] = partPositions_[ii][2];
                	particleCreateParcelList[ii] = ii;
	}

        Pout << tab << "Creating kdTree..." << endl;
        kdTree = new ANNkd_tree(dataPts, nPts, dim);

	Pout << tab << "Entering particle loops to create parcel" << endl; 
        Pout << tab << "Number of particles " << partNP_ << " to create parcels " << endl; 

    	label parcelI = 0;
    	for(int ii =0; ii < partNP_; ii++)
	{				

        	if ( particleCreateParcelList[ii] > -1 )
        	{
        	    queryPt[0] = partPositions_[ii][0];
        	    queryPt[1] = partPositions_[ii][1];
        	    queryPt[2] = partPositions_[ii][2];

        	    
		    /*
		    kdTree->annkFRSearch(
                                	    queryPt,				// query point					
                                	    sqRad,				// squared radius
                                	    k,                  		// number of the near neighbours to return
                                	    nnIdx,				// nearest neighbor array
                                	    dists,				// dist to near neighbours
                                	    eps			);

        	    */
		
		    	
		    kdTree->annkSearch(
                                	    queryPt,				// query point					
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
			    /*
			    if( ii%100000 == 0 )
			    {
			       
			       Info  << tab << " parcelI " << parcelI << " partSum " << partSum << " Neighbour part " 
			             << nnIdx[i] << " particlesInSubDomain " << partNP_ //particlesInSubDomain_[nnIdx[i]] 
			             << " particlesIn " << particlesIn_[parcelI][partSum] << endl;			    
			    }
			    */
                	    if (!useAllParticles_) particleCreateParcelList[nnIdx[i]] = -1 ;
                	    particlesIn_[parcelI][partSum] = nnIdx[i];
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

     // Find closest particles
     findParticlesIn();	

     // Parcel radius
     parcelRadius();

     // Calculate parcel velocity
     parcelVel();

     // Locate parcel
     parcelPos();

     // Find parcel cell
     parcelLocate();

     // Write parcels into
     // write();	
}

void Foam::parcelCloud::readParcel()
{
     // Read parcels from file
     read();	
}

vectorField Foam::parcelCloud::parcelVel()
{

    for(label parI=0; parI < particlesIn_.size(); parI++)
    {
	    for(label partIn=0; partIn < nPInParcel_; partIn++)
	    {
		    velocities_[parI][0] += partVelocities_[particlesIn_[parI][partIn]][0];
		    velocities_[parI][1] += partVelocities_[particlesIn_[parI][partIn]][1];
		    velocities_[parI][2] += partVelocities_[particlesIn_[parI][partIn]][2];

		    //Info << parI << " " << partIn << " " << particlesIn_[parI][partIn] << endl;
	    }
	    velocities_[parI] /= nPInParcel_;
    }	

    return velocities_;
}

vectorField Foam::parcelCloud::parcelPos()
{

    for(label parI=0; parI < particlesIn_.size(); parI++)
    {
	    for(label partIn=0; partIn < nPInParcel_; partIn++)
	    {
		    positions_[parI][0] += partPositions_[particlesIn_[parI][partIn]][0];
		    positions_[parI][1] += partPositions_[particlesIn_[parI][partIn]][1];
		    positions_[parI][2] += partPositions_[particlesIn_[parI][partIn]][2];
	    }
	    positions_[parI] /= nPInParcel_;
    }	

    return positions_;
}

scalarField Foam::parcelCloud::parcelRadius()
{

    for(int index = 0;index < nPParcel_; ++index)
    {
        radius_[index] = partRadii_[index];
    }
	    
    return radius_;
}

labelList Foam::parcelCloud::parcelLocate()
{
    meshSearch searchEngine(refCast<const fvMesh>(obr_),polyMesh::FACEPLANES);
    
    for(int index = 0;index < nPParcel_; ++index)
    {
	  if(radius_[index] > SMALL)
	  {
              // create position vector
              vector pos = positions_[index];
              // find cell
              cellIDs_[index] = searchEngine.findCell(pos,cellIDs_[index],treeSearch_);
	      //Info << " index = " << index << " cellIDs = " << cellIDs_[index] << endl;
	  }
    }
	    
    return cellIDs_;
}

//- Write parcels into folders
void Foam::parcelCloud::write()
{
    
    Info << tab << "Parcels are writing into files" << endl;
    //- Create name for parcel 			  
    char charName[100];
    sprintf(charName, "nP%d",nPInParcel_);
    word parcelName(charName);
			  
    IOField<vector> parcelVel
    ( 
      IOobject
      (
	"parcelVel"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::NO_READ,
	IOobject::AUTO_WRITE 
      ),
      velocities_
    );

    IOField<vector> parcelPos
    ( 
      IOobject
      (
	"parcelPos"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::NO_READ,
	IOobject::AUTO_WRITE 
      ),
      positions_
    );

    IOField<scalar> parcelRad
    ( 
      IOobject
      (
	"parcelRad"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::NO_READ,
	IOobject::AUTO_WRITE 
      ),
      radius_
    );
    
    IOList<label> parcelCellID
    ( 
      IOobject
      (
	"parcelCellID"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::NO_READ,
	IOobject::AUTO_WRITE 
      ),
      cellIDs_
    ); 

    IOList< List < label > > particlesInParcel
    ( 
      IOobject
      (
	"particlesInParcel"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::NO_READ,
	IOobject::AUTO_WRITE 
      ),
      particlesIn_
    ); 
                   
    //parcelVel.write();
    //parcelPos.write();
    //parcelRad.write();
    //parcelCellID.write();
    //particlesInParcel.write();
     
}

void Foam::parcelCloud::read()
{
    //Info << tab << "Parcels are reading files" << endl;
    //- Create name for parcel 			  
    char charName[100];
    sprintf(charName, "nP%d",nPInParcel_);
    word parcelName(charName);
    
    IOField<vector> parcelVel
    ( 
      IOobject
      (
	"parcelVel"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::MUST_READ,
	IOobject::NO_WRITE 
      )
    );

    IOField<vector> parcelPos
    ( 
      IOobject
      (
	"parcelPos"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::MUST_READ,
	IOobject::NO_WRITE 
      )
    );

    IOField<scalar> parcelRad
    ( 
      IOobject
      (
	"parcelRad"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::MUST_READ,
	IOobject::NO_WRITE 
      )
    );
    
    IOList<label> parcelCellID
    ( 
      IOobject
      (
	"parcelCellID"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::MUST_READ,
	IOobject::NO_WRITE 
      )
    );

    IOList< List < label > > particlesInParcel
    ( 
      IOobject
      (
	"particlesInParcel"+parcelName,
	refCast<const fvMesh>(obr_).time().timeName(),
	refCast<const fvMesh>(obr_),
	IOobject::MUST_READ,
	IOobject::NO_WRITE 
      ),
      particlesIn_
    ); 
            
    // Set number of parcels
    nPParcel_ = parcelVel.size();
    
    //- Resize vectors     
    velocities_.resize(nPParcel_);
    positions_.resize(nPParcel_);
    radius_.resize(nPParcel_);
    cellIDs_.resize(nPParcel_);
    particlesIn_.resize(nPParcel_);
       
    // Set variables
    velocities_ = parcelVel;
    positions_ = parcelPos;
    radius_ = parcelRad;
    cellIDs_ = parcelCellID;          
    particlesIn_ = particlesInParcel;
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
