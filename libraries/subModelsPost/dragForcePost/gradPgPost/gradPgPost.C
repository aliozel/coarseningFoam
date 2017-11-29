/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    Carrier phase pressure gradient force

Class
     Carrier phase pressure gradient force

SourceFiles
    gradPgPost.C

-----------------------------------------------------------------------------*/

#include "error.H"

#include "gradPgPost.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradPgPost, 0);

addToRunTimeSelectionTable
(
    dragModelPost,
    gradPgPost,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gradPgPost::gradPgPost
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    dragModelPost(dict,obr),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    interpolation_(false)
{
    if (propsDict_.found("interpolation")) 
    {
        interpolation_=true;
        Info << tab << "Using interpolated value of gradPg" << endl;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradPgPost::~gradPgPost()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradPgPost::setForce
(
    const volVectorField& gradPg_,
    cfdemCloud& particleCloud_
) const
{
    // Resize arrays
    //force_.resize(particleCloud_.numberOfParticles());
    gradPgForce_.resize(particleCloud_.numberOfParticles());
    
    vector gradP(0,0,0);
    scalar ds(0);
    scalar Vs(0);
    vector position(0,0,0);
    vector force(0,0,0);
    label cellI(-1);

    interpolationCellPoint<vector> gradPgInterpolator(gradPg_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
       cellI = particleCloud_.cellIDs()[index][0];

       if (cellI > -1) // particle Found
       {
           position = particleCloud_.position(index);

           if(interpolation_) 
           {
               gradP = gradPgInterpolator.interpolate(position,cellI);
           }else
           {
               gradP = gradPg_[cellI];
           }

           ds = 2*particleCloud_.radius(index);
           Vs = ds*ds*ds*M_PI/6;

           // calc particle's pressure gradient force
           force = -Vs*gradP;

           if(verbose_ && index >=0 && index <2)
           {
               //Pout << "periodicPressure " << endl; 
	       Pout << "index = " << index << endl;
               Pout << "gradP = " << gradP << endl;
               Pout << "force = " << force << endl;
           }
	   
	   // Drag force
      	   for(int j=0;j<3;j++) gradPgForce_[index][j] = force[j];
           //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

       }
    }
}

void gradPgPost::setForceParcel
(
    const volVectorField& gradPg_,
    parcelCloud& parcelCloud_
) const
{
    // Resize arrays
    //force_.resize(parcelCloud_.numberOfParticles());
    gradPgForce_.resize(parcelCloud_.numberOfParticles());
    
    vector gradP(0,0,0);
    scalar ds(0);
    scalar Vs(0);
    vector position(0,0,0);
    vector force(0,0,0);
    label cellI(-1);

    interpolationCellPoint<vector> gradPgInterpolator(gradPg_);

    for(int index = 0;index <  parcelCloud_.numberOfParticles(); index++)
    {
       cellI = parcelCloud_.cellIDs(index);

       position = parcelCloud_.position(index);

       if(interpolation_) 
       {
           gradP = gradPgInterpolator.interpolate(position,cellI);
       }else
       {
           gradP = gradPg_[cellI];
       }

       ds = 2*parcelCloud_.radius(index);
       Vs = ds*ds*ds*M_PI/6;

       // calc particle's pressure gradient force
       force = -Vs*gradP;

       if(verbose_ && index >=0 && index <2)
       {
           //Pout << "periodicPressure " << endl; 
	   Pout << "index = " << index << " gradP = " << force << endl;
       }

       // Drag force
       for(int j=0;j<3;j++) gradPgForce_[index][j] = force[j];
       //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

   }
}

void gradPgPost::setForce
(
    const volVectorField& gradPg_,
    const int & partNP_,
    labelList& partCellIDs_,
    scalarList& partRadii_	
) const
{
    // Resize arrays
    //force_.resize(particleCloud_.numberOfParticles());
    gradPgForce_.resize(partNP_);
    
    vector gradP(0,0,0);
    scalar ds(0);
    scalar Vs(0);
    vector position(0,0,0);
    vector force(0,0,0);
    label cellI(-1);

    interpolationCellPoint<vector> gradPgInterpolator(gradPg_);

    for(int index = 0;index < partNP_; index++)
    {
       cellI = partCellIDs_[index];

       if (cellI > -1) // particle Found
       {

           if(interpolation_) 
           {
               gradP = gradPgInterpolator.interpolate(position,cellI);
           }else
           {
               gradP = gradPg_[cellI];
           }

           ds = partRadii_[index];
           Vs = ds*ds*ds*M_PI/6;

           // calc particle's pressure gradient force
           force = -Vs*gradP;

           if(verbose_ && index >=0 && index <2)
           {
               //Pout << "periodicPressure " << endl; 
	       Pout << "index = " << index << endl;
               Pout << "gradP = " << gradP << endl;
               Pout << "force = " << force << endl;
           }
	   
	   // Drag force
      	   for(int j=0;j<3;j++) gradPgForceEulerian_[index][j] = force[j];
           //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

       }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
