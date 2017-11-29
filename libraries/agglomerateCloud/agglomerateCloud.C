/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
		 Create agglomerate

-----------------------------------------------------------------------------*/

#include "agglomerateCloud.H"

namespace Foam {
    defineTypeNameAndDebug(agglomerateCloud,0);
}

Foam::agglomerateCloud::agglomerateCloud
(
    const dictionary& dict,
    const objectRegistry& obr,
    cfdemCloud& sm,
    const int& count,
    bool& locateAgglomerate
)
:
    dict_(dict),
    obr_(obr),
    sm_(sm),
    count_(count),
    locateAgglomerate_(locateAgglomerate),
    nPInAgglomerate_(sm_.numberOfParticles(),-1),
    particlesIn_(sm_.numberOfParticles(),labelList(nPInAgglomerate_,-1)),
    velocities_(sm_.numberOfParticles(),vector(0,0,0)),
    positions_(sm_.numberOfParticles(),vector(0,0,0)),
    radius_(sm_.numberOfParticles(),scalar(0)),
    types_(sm_.numberOfParticles(),-1),
    cellIDs_(sm_.numberOfParticles(),-1),
    readFromFile_(false),
    nPAgglomerate_(0),
    agglomerateVTKPropsDict_(dict_.subDict("oneWayVTKProps")),
    filename_(agglomerateVTKPropsDict_.lookup("agglomerateFilename")),
    relativePath_(agglomerateVTKPropsDict_.lookup("agglomerateRelativePath")),
    connectivityList_(sm_.numberOfParticles(),labelList(10,-1)),    
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
    ),
    bridgeIDLen_(10)
{
    // Init
    init();
}

Foam::agglomerateCloud::~agglomerateCloud()
{}

void Foam::agglomerateCloud::init()
{
    // Set size label list
    nPInAgglomerate_.setSize(sm_.numberOfParticles(),0);

    // Define agglomerate dictionary
    string HH=string(filename_);
    charFilename_=HH.c_str();
    Info << tab << "Relative path for agglomerate data: " << relativePath_ << endl;
}

Foam::labelListList Foam::agglomerateCloud::readConnectivityList()
{
    // get path to particle connectivity files
    char index[100];
    sprintf(index, charFilename_,count_);
    fileName H(sm_.mesh().time().path()/relativePath_/index);
    Info << tab << "Opening file: " << H << endl;	

    // Set file pointer
    string HH=string(H);
    const char * paricleFilePath=HH.c_str();
    ifstream* inputPtr;
    inputPtr = new ifstream(paricleFilePath);

    // Connectivity list
    int ** connectivityList;
	
    // Read data header
    string just_read = " ";
    for(int index = 0;index < 38; ++index)
    {
    	*inputPtr >> just_read;   // ITEM: ATOMS id f_connectivityList[0:9]
    	//Info      << " just_read = " << just_read << endl; 	
    }
    // Read connectivity list
    for(int index = 0;index < sm_.numberOfParticles(); ++index)
    {
        *inputPtr >> connectivityList_[index][0] >> connectivityList_[index][1] >> connectivityList_[index][2]
	          >> connectivityList_[index][3] >> connectivityList_[index][4] >> connectivityList_[index][5]
		  >> connectivityList_[index][6] >> connectivityList_[index][7] >> connectivityList_[index][8]
		  >> connectivityList_[index][9] ;   		  
	// Info << "connectivityList_["<< index <<"][0] " <<  connectivityList_[index][0] << endl;     
    }
    
    return connectivityList_;
} 

Foam::labelListList Foam::agglomerateCloud::findParticlesIn()
{			
    //label agglomerateI = 0;	

    // Sorted IDs
    labelList sortedID(sm_.numberOfParticles());

    // Sorted position
    labelList sortedPosition(sm_.numberOfParticles());

    // Bridge IDs
    // connectivityList_

    // Agglomerate IDs
    labelList agglomID(sm_.numberOfParticles());

    // Agglomerate Sizes
    labelList agglomerateSizes(floor(sm_.numberOfParticles()/2));

    label swapCount(1);
    label agglomerateCount(0);
    label subID;
    label swappedParticleID;
    label temp; 

    for(int index = 0;index < sm_.numberOfParticles(); ++index)
    {
       label hostID = sortedID[index];
       // Check if particle is in an agglomerate
       if(agglomID[hostID+1]>=0) 
       {
	  swapCount = swapCount - 1;
	  // Traverse contacts
	  for(int conn = 1; conn < bridgeIDLen_; ++conn)
	  {
	     if(connectivityList_[hostID+1][conn] != -1) 	//  look only at contacts
	     {
		subID = connectivityList_[hostID+1][conn];
		if(agglomID[subID+1] == -1)			// Not declared in an agglomerate
		{
		   agglomID[subID+1] = agglomID[hostID+1];
		   agglomerateSizes[agglomID[hostID+1]] = agglomerateSizes[agglomID[hostID+1]+1]+1; // add 1 to aglomerate size
		   // Perform swap
		   swappedParticleID = sortedID[index+swapCount];
		   sortedID[index+swapCount] = subID;
		   sortedID[sortedPosition[subID+1]] = swappedParticleID; 
		   temp = sortedPosition[swappedParticleID+1];
		   sortedPosition[swappedParticleID+1] = sortedPosition[subID+1];
		   sortedPosition[subID+1] = temp;
		   swapCount = swapCount + 1;
		} 
	     }
	  }	
       } 
       // if particle is not in an agglomerate
       else	
       {
	  label firstContactSwitch(0);
	  for(int conn = 1; conn < bridgeIDLen_; ++conn)
	  {
	     if(connectivityList_[hostID+1][conn] != -1) //  look only at contacts
	     {
		if(firstContactSwitch == 0)
		{
		   agglomerateCount++;
		   agglomerateSizes[agglomerateCount] = agglomerateSizes[agglomerateCount] + 1;
		   agglomID[hostID+1] = agglomerateCount - 1;
		   firstContactSwitch = -1;  
		}
		subID = connectivityList_[hostID+1][conn];
		if(agglomID[subID+1] == -1)
		{
		   agglomID[subID+1] = agglomID[subID+1];
		   agglomerateSizes[agglomID[hostID+1]] = agglomerateSizes[agglomID[hostID+1]+1]+1; // add 1 to aglomerate size
		   // Perform swap
		   swappedParticleID = sortedID[index+swapCount];
		   sortedID[index+swapCount] = subID;
		   sortedID[sortedPosition[subID+1]] = swappedParticleID; 
		   temp = sortedPosition[swappedParticleID+1];
		   sortedPosition[swappedParticleID+1] = sortedPosition[subID+1];
		   sortedPosition[subID+1] = temp;
		   swapCount = swapCount + 1;		       
		}
	     }		   
	  }		   
       }	 
    }

    // Resize number of agglomerate
    //nPAgglomerate_ = agglomerateI;
    particlesIn_.resize(nPAgglomerate_);
    velocities_.resize(nPAgglomerate_);
    positions_.resize(nPAgglomerate_);
    radius_.resize(nPAgglomerate_);
    types_.resize(nPAgglomerate_);
    cellIDs_.resize(nPAgglomerate_); 

    return particlesIn_;	 
} 

void Foam::agglomerateCloud::createAgglomerate()
{
    if(!readFromFile_) 
    {
	// Read connectivity list
	readConnectivityList();
	// Find closest particles
	findParticlesIn();
	// Agglomerate radius
	agglomerateRadius();
	// Calculate Agglomerate velocity
	agglomerateVel();
	// Locate Agglomerate
	agglomeratePos();	
	// Find Agglomerate cell
        if(locateAgglomerate_)
		//locateM().findCell(NULL,positions_,cellIDs_,particlesInSubDomain_.size(),radius_);
		locateM().findCell(NULL,positions_,cellIDs_,sm_.numberOfParticles(),radius_);    
    }else
    {
    	// Do nothing
    }	
}

vectorField Foam::agglomerateCloud::agglomerateVel()
{
    if(!readFromFile_) 
    {
	for(label parI=0; parI < particlesIn_.size(); parI++)
	{
		for(label partIn=0; partIn < nPInAgglomerate_[parI]; partIn++)
		{
			velocities_[parI] += sm_.velocity(particlesIn_[parI][partIn]);
		}
		velocities_[parI] /= nPInAgglomerate_[parI];
	}	
    }else
    {
    	// Do nothing
    }
    return velocities_;
}

vectorField Foam::agglomerateCloud::agglomeratePos()
{
    if(!readFromFile_) 
    {
	for(label parI=0; parI < particlesIn_.size(); parI++)
	{
		for(label partIn=0; partIn < nPInAgglomerate_[parI]; partIn++)
		{
			positions_[parI] += sm_.position(particlesIn_[parI][partIn]);
		}
		positions_[parI] /= nPInAgglomerate_[parI];
	}	
    }else
    {
    	// Do nothing
    }
    return positions_;
}

scalarField Foam::agglomerateCloud::agglomerateRadius()
{
    if(!readFromFile_) 
    {
	    for(int index = 0;index < nPAgglomerate_; ++index)
	    {
        	radius_[index] = sm_.radius(index);
	    }
    }else
    {
     	// Do nothing
    }		    
    return radius_;
}

void Foam::agglomerateCloud::read()
{
    if(readFromFile_) 
    {
	// Do nothing
    }

}

void Foam::agglomerateCloud::write()
{
    if(!readFromFile_) 
    {
	// Do nothing
    }

}

// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

int Foam::agglomerateCloud::numberOfAgglomerates()
{
    return nPAgglomerate_;		
}

vector Foam::agglomerateCloud::velocity(int& index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = velocities_[index][i];
    return vel;
}

vector Foam::agglomerateCloud::position(int& index)
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positions_[index][i];
    return pos;
}

scalar Foam::agglomerateCloud::radius(int& index)
{
    scalar rad;
    rad = radius_[index];
    return rad; 	
}

label Foam::agglomerateCloud::type(int& index)
{
    label ty = static_cast<label>(types_[index]);
    return ty; 	
}

label Foam::agglomerateCloud::cellIDs(int& index)
{
    label cellI = cellIDs_[index];
    return cellI;
}

labelList Foam::agglomerateCloud::particlesInAgglomerate(int& index)
{
    labelList partsIn = particlesIn_[index]; 
    return partsIn;	
}

// ************************************************************************* //
