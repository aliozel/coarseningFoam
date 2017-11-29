/*------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description 
    WenYu drag model

Class
    WenYuDragPost

SourceFiles
    WenYuDragPost.C

-----------------------------------------------------------------------------*/

#include "error.H"
#include "WenYuDragPost.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WenYuDragPost, 0);

addToRunTimeSelectionTable
(
    dragModelPost,
    WenYuDragPost,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
WenYuDragPost::WenYuDragPost
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    dragModelPost(dict,obr),
    propsDict_(dict.subDict(typeName + "Props")),
    //rhoParticle_(readScalar(propsDict_.lookup("rhoParticle"))),    
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(refCast<const fvMesh>(obr).lookupObject<volScalarField> (densityFieldName_)),
    verbose_(false),
    interpolation_(false),
    maxCellPerPart_(0),
    taupFuncPhi_(false),
    taupConst_(false)    
{
    rhoParticle_ = readScalar(propsDict_.lookup("rhoParticle"));
    if (propsDict_.found("interpolation")) 
    {
        interpolation_=true;
        Info << tab << "Using interpolated value of U" << endl;
    }
    if (propsDict_.found("verbose")) 
    {
        verbose_=true;
    }
    if (propsDict_.found("maxCellPerPart")) 
    {
        maxCellPerPart_ = readScalar(propsDict_.lookup("maxCellPerPart"));    
        Info << tab << "Max. number of cells per particle " << maxCellPerPart_ 
		    << " for mapping & interpolation" << endl;
    }
    if (propsDict_.found("taupFuncPhi")) 
    {
        taupFuncPhi_= true;
	taupConst_ = false;
        Info << tab << "Taup is only function of phi " << endl;	
    }
    if (propsDict_.found("taupConst")) 
    {
        taupConst_= true;
	taupFuncPhi_= true; 
        Info << tab << "Taup is only function of phi" << endl;	
    }
    Switch weighthing(false);
    dict.lookup("weighthing") >> weighthing;
    if(!weighthing) 
    {
        maxCellPerPart_ = 0; //1;    
        Info << tab << "Max. number of cells per particle " << maxCellPerPart_ 
		    << " for mapping & interpolation + weighting will not be used" << endl;	
    }          
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

WenYuDragPost::~WenYuDragPost()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WenYuDragPost::setForce
(
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_
) const
{
    // Resize arrays
    dragForce_.resize(particleCloud_.numberOfParticles());
    taup_.resize(particleCloud_.numberOfParticles());
    Ufxp_.resize(particleCloud_.numberOfParticles());
    cellIDs_.resize(particleCloud_.numberOfParticles());
    velocities_.resize(particleCloud_.numberOfParticles());
    positions_.resize(particleCloud_.numberOfParticles());	
    voidfractionxp_.resize(particleCloud_.numberOfParticles());
	
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label  cellI=0;

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar alps(0);
    scalar betaP(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    
    //- Number of particle for drag autoPtr 
    nP_ = particleCloud_.numberOfParticles();
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {	   
      cellI = particleCloud_.cellIDs()[index][0];
      cellIDs_[index] = cellI; 

      if (cellI > -1) // particle Found
      {
          position = particleCloud_.position(index);
	  positions_[index] = position;
	  if(interpolation_)
          {	            
              voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
              Ufluid = UInterpolator_.interpolate(position,cellI);
              //Ensure interpolated void fraction to be meaningful
              // Info << " --> voidfraction: " << voidfraction << endl;
              if(voidfraction>1.00) voidfraction = 1.00;
              if(voidfraction<0.40) voidfraction = 0.40;
          }
	  else
          {
	      voidfraction = particleCloud_.voidfraction(index);
              Ufluid = U_[cellI];
          }

	  // Fluid velocity at particle position
	  Ufxp_[index] = Ufluid;
	  // Fluid volume fraction at particle position
	  voidfractionxp_[index] = voidfraction;

          Up = particleCloud_.velocity(index);
	  velocities_[index] = Up;
          Ur = Ufluid-Up;
          ds = 2*particleCloud_.radius(index);

	  nuf = nufField[cellI];
          rho = rho_[cellI];
          magUr = mag(Ur);
          Volp = ds*ds*ds*M_PI/6;
          alps = 1-voidfraction+SMALL;

          if (magUr > 0)
          {
              // calc particle Re number
		Rep = voidfraction*ds*magUr/(nuf+SMALL);

              // calc CD
	      if (Rep < 1000)
	      {
		CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	      }else
	      {
        	CD = 0.44*pow(voidfraction,-1.65); 
	      }  

	      // calc drag coefficient 
	      betaP = 3./4.*rho*CD*magUr/ds;  	

              // Relaxation time
	      taup_[index] = rhoParticle_/betaP;

	      // calc particle's drag
              drag = Volp*betaP*Ur;
	  }

          if(verbose_ && index <1 )
          {
	      Pout << tab << "WenYu Drag"  << endl;
	      Pout << tab << "index = " << index << endl;
              Pout << tab << "Up = " << Up << endl;
              Pout << tab << "Ur = " << Ur << endl;
              Pout << tab << "dp = " << ds << endl;
              Pout << tab << "dparcel = " << 2*particleCloud_.radius(index) << endl;
              Pout << tab << "rho = " << rho << endl;
              Pout << tab << "nuf = " << nuf << endl;
              Pout << tab << "voidfraction = " << voidfraction << endl;
              Pout << tab << "Rep = " << Rep << endl;
              Pout << tab << "Volp = " << Volp << endl;
              Pout << tab << "alps = " << alps << endl;
              Pout << tab << "CD = " << CD << endl;
	      Pout << tab << "betaP = " << betaP << endl;		   		    
              Pout << tab << "drag = " << drag << endl;
	      Pout << "" << endl;
          }
      }

      // Drag force
      for(int j=0;j<3;j++) dragForce_[index][j] = drag[j];		  

    }

} 

void WenYuDragPost::setForceEulerian
(
	const volScalarField& voidfraction_,
	const volVectorField& U_,
	cfdemCloud& particleCloud_,
	const volVectorField& Us_
) const
{
    // Resize arrays
    dragForceEulerian_.resize(U_.mesh().cells().size());
    taupEulerian_.resize(U_.mesh().cells().size());
	
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector Ur(0,0,0);
    vector drag(0,0,0);
    scalar voidfraction(0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar alps(0);

    forAll(U_.mesh().cells(),cellI)   	    
    {	   
    
      voidfraction = voidfraction_[cellI];
      Ur = U_[cellI]-Us_[cellI];
      ds = 2*particleCloud_.radius(0); // WARNING BE CAREFUL FOR BIDISPERSE CASE
      nuf = nufField[cellI];
      rho = rho_[cellI];
      magUr = mag(Ur);
      alps = 1-voidfraction+SMALL;

      if (magUr > 0)
      {
          // calc particle Re number
	  Rep = voidfraction_[cellI]*ds*magUr/(nuf+SMALL);

          // calc CD
	  if (Rep < 1000)
	  {
	    CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	  }else
	  {
            CD = 0.44*pow(voidfraction,-1.65); 
	  }  

          // Relaxation time
	  taupEulerian_[cellI] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);

	  // calc particle's drag
          drag = alps*rhoParticle_/taupEulerian_[cellI]*Ur;
      }

      if(verbose_ && cellI <1 )
      {
	  Pout << tab << "WenYu Drag"  << endl;
	  Pout << tab << "index = " << cellI << endl;
          Pout << tab << "Up = " << Us_[cellI] << endl;
          Pout << tab << "Ur = " << Ur << endl;
          Pout << tab << "dp = " << ds << endl;
          Pout << tab << "rho = " << rho << endl;
          Pout << tab << "nuf = " << nuf << endl;
          Pout << tab << "voidfraction = " << voidfraction << endl;
          Pout << tab << "Rep = " << Rep << endl;
          Pout << tab << "alps = " << alps << endl;
          Pout << tab << "CD = " << CD << endl;
          Pout << tab << "drag = " << drag << endl;
	  Pout << "" << endl;
      }

      // Drag force
      for(int j=0;j<3;j++) dragForceEulerian_[cellI][j] = drag[j];		  

    }

}

// Overwrite function for not particle cloud
void WenYuDragPost::setForce
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    const int & partNP_,
    labelList& partCellIDs_,
    vectorField& partPositions_,
    vectorField& partVelocities_,
    scalarList& partRadii_
) const
{
    // Resize arrays
    dragForce_.resize(partNP_);
    taup_.resize(partNP_);
    Ufxp_.resize(partNP_);
    cellIDs_.resize(partNP_);
    velocities_.resize(partNP_);
    positions_.resize(partNP_);	
    voidfractionxp_.resize(partNP_);
	
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI(-1);

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar alps(0);
    scalar betaP(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    
    //- Number of particle for drag autoPtr 
    nP_ = partNP_;
    
    for(int index = 0;index <  partNP_; index++)
    {	   
      cellI = partCellIDs_[index];
      cellIDs_[index] = cellI; 

      if (cellI > -1) // particle Found
      {
          for(int idir=0; idir<3; idir++) position[idir] = partPositions_[index][idir];
	  positions_[index] = position;
	  if(interpolation_)
          {	            
              voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
              Ufluid = UInterpolator_.interpolate(position,cellI);
              //Ensure interpolated void fraction to be meaningful
              // Info << " --> voidfraction: " << voidfraction << endl;
              if(voidfraction>1.00) voidfraction = 1.00;
              if(voidfraction<0.40) voidfraction = 0.40;
          }
	  else
          {
	      voidfraction = voidfraction_[cellI];
              Ufluid = U_[cellI];
          }

	  // Fluid velocity at particle position
	  Ufxp_[index] = Ufluid;
	  // Fluid volume fraction at particle position
	  voidfractionxp_[index] = voidfraction;

          for(int idir=0; idir<3; idir++) Up[idir] = partVelocities_[index][idir];
	  velocities_[index] = Up;
          Ur = Ufluid-Up;
          ds = partRadii_[index];

	  nuf = nufField[cellI];
          rho = rho_[cellI];
          magUr = mag(Ur);
          Volp = ds*ds*ds*M_PI/6;
          alps = 1-voidfraction+SMALL;

          if (magUr > 0)
          {
              // calc particle Re number
		Rep = voidfraction*ds*magUr/(nuf+SMALL);

              // calc CD
	      if (Rep < 1000)
	      {
		CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	      }else
	      {
        	CD = 0.44*pow(voidfraction,-1.65); 
	      }  

	      // calc drag coefficient 
	      betaP = 3./4.*rho*CD*magUr/ds;  	

              // Relaxation time
	      taup_[index] = rhoParticle_/betaP;

	      // calc particle's drag
              drag = Volp*betaP*Ur;
	  }

          if(verbose_ && index <1 )
          {
	      Pout << tab << "WenYu Drag"  << endl;
	      Pout << tab << "index = " << index << endl;
              Pout << tab << "Up = " << Up << endl;
              Pout << tab << "Ur = " << Ur << endl;
              Pout << tab << "dp = " << ds << endl;
              Pout << tab << "rho = " << rho << endl;
              Pout << tab << "nuf = " << nuf << endl;
              Pout << tab << "voidfraction = " << voidfraction << endl;
              Pout << tab << "Rep = " << Rep << endl;
              Pout << tab << "Volp = " << Volp << endl;
              Pout << tab << "alps = " << alps << endl;
              Pout << tab << "CD = " << CD << endl;
	      Pout << tab << "betaP = " << betaP << endl;		   		    
              Pout << tab << "drag = " << drag << endl;
	      Pout << "" << endl;
          }
      }

      // Drag force
      for(int j=0;j<3;j++) dragForce_[index][j] = drag[j];		  

    }

} 

void WenYuDragPost::setForceEulerian
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    scalarField& partRadii_,
    const volVectorField& Us_
) const
{
    // Resize arrays
    dragForceEulerian_.resize(U_.mesh().cells().size());
    taupEulerian_.resize(U_.mesh().cells().size());

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector Ur(0,0,0);
    vector drag(0,0,0);
    scalar voidfraction(0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar alps(0);

    forAll(U_.mesh().cells(),cellI)   	    
    {	   
    
      voidfraction = voidfraction_[cellI];
      Ur = U_[cellI]-Us_[cellI];
      ds = 2.*partRadii_[0]; // WARNING BE CAREFUL FOR BIDISPERSE CASE
      nuf = nufField[cellI];
      rho = rho_[cellI];
      magUr = mag(Ur);
      alps = 1-voidfraction+SMALL;

      if (magUr > 0)
      {
          // calc particle Re number
	  Rep = voidfraction_[cellI]*ds*magUr/(nuf+SMALL);

          // calc CD
	  if (Rep < 1000)
	  {
	    CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	  }else
	  {
            CD = 0.44*pow(voidfraction,-1.65); 
	  }  

          // Relaxation time
	  taupEulerian_[cellI] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);
	  
	  // calc particle's drag
          drag = alps*rhoParticle_/taupEulerian_[cellI]*Ur;
      }

      if(verbose_ && cellI <1 )
      {
	  Pout << tab << "WenYu Drag-scalarField radii"  << endl;
	  Pout << tab << "index = " << cellI << endl;
          Pout << tab << "Up = " << Us_[cellI] << endl;
          Pout << tab << "Ur = " << Ur << endl;
          Pout << tab << "dp = " << ds << endl;
          Pout << tab << "rho = " << rho << endl;
          Pout << tab << "nuf = " << nuf << endl;
          Pout << tab << "voidfraction = " << voidfraction << endl;
          Pout << tab << "Rep = " << Rep << endl;
          Pout << tab << "alps = " << alps << endl;
          Pout << tab << "CD = " << CD << endl;
          Pout << tab << "drag = " << drag << endl;
	  Pout << "" << endl;
      }

      // Drag force
      for(int j=0;j<3;j++) dragForceEulerian_[cellI][j] = drag[j];		  

    }

}

//-GradPg force
void WenYuDragPost::setForce
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
               Pout << tab << "Gas pressure gradient"  << endl;
	       Pout << "index = " << index << endl;
               Pout << "gradP = " << gradP << endl;
               Pout << "force = " << force << endl;
	       Pout << "" << endl;
           }
	   
	   // Drag force
      	   for(int j=0;j<3;j++) gradPgForce_[index][j] = force[j];
           //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

       }
    }
}

//-New
//-functions
//-
//-
//-
//-
//-
//-

// Overwrite function for List<vector>
void WenYuDragPost::setForce
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    const int & partNP_,
    labelList& partCellIDs_,
    List<vector>& partPositions_,
    List<vector>& partVelocities_,
    List<scalar>& partRadii_
) const
{

    // Resize arrays
    dragForce_.resize(partNP_);
    taup_.resize(partNP_);
    Ufxp_.resize(partNP_);
    cellIDs_.resize(partNP_);
    velocities_.resize(partNP_);
    positions_.resize(partNP_);	
    voidfractionxp_.resize(partNP_);
	
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI(-1);

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar alps(0);
    scalar mp(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    
    //- Number of particle for drag autoPtr 
    nP_ = partNP_;
    

    for(int index = 0;index <  partNP_; index++)
    {	   
      cellI = partCellIDs_[index];
      cellIDs_[index] = cellI; 

      if (cellI > -1) // particle Found
      {
          for(int idir=0; idir<3; idir++) position[idir] = partPositions_[index][idir];
	  positions_[index] = position;
	  if(interpolation_)
          {	            
              voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
              Ufluid = UInterpolator_.interpolate(position,cellI);
              //Ensure interpolated void fraction to be meaningful
              // Info << " --> voidfraction: " << voidfraction << endl;
              if(voidfraction>1.00) voidfraction = 1.00;
              if(voidfraction<0.40) voidfraction = 0.40;
          }
	  else
          {
	      voidfraction = voidfraction_[cellI];
              Ufluid = U_[cellI];
          }

	  // Fluid velocity at particle position
	  Ufxp_[index] = Ufluid;
	  // Fluid volume fraction at particle position
	  voidfractionxp_[index] = voidfraction;

          for(int idir=0; idir<3; idir++) Up[idir] = partVelocities_[index][idir];
	  velocities_[index] = Up;
          Ur = Ufluid-Up;
          ds = 2.*partRadii_[index];

	  nuf = nufField[cellI];
          rho = rho_[cellI];
          magUr = mag(Ur);
          Volp = ds*ds*ds*M_PI/6;
	  mp = rhoParticle_*Volp;
          alps = 1-voidfraction+SMALL;

          if (magUr > 0)
          {
              // calc particle Re number
              Rep = voidfraction*ds*magUr/(nuf+SMALL);

              // calc CD
	      if (Rep < 1000)
	      {
		CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	      }else
	      {
        	CD = 0.44*pow(voidfraction,-1.65); 
	      }  

	      //- Non linearity test
	      if(taupFuncPhi_) CD = 24./Rep*pow(voidfraction,-1.65);
	  	
              // Relaxation time
	      taup_[index] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);

	      // Test correlation alppUr
	      if(taupConst_) taup_[index] = 0.025;

	      // calc particle's drag
              drag = mp/taup_[index]*Ur;
	     
	  }

          if(verbose_ && index <1 )
          {
	      Pout << tab << "WenYu Drag-Lists"  << endl;
	      Pout << tab << "index = " << index << endl;
	      Pout << tab << "cellI = " << cellI << endl;
              Pout << tab << "Up = " << Up << endl;
              Pout << tab << "Ur = " << Ur << endl;
              Pout << tab << "dp = " << ds << endl;
              Pout << tab << "rho = " << rho << endl;
              Pout << tab << "nuf = " << nuf << endl;
              Pout << tab << "voidfraction = " << voidfraction << endl;
              Pout << tab << "Rep = " << Rep << endl;
              Pout << tab << "Volp = " << Volp << endl;
              Pout << tab << "alps = " << alps << endl;
              Pout << tab << "CD = " << CD << endl;
	      Pout << tab << "taup = " << taup_[index] << endl;		   		    
              Pout << tab << "drag = " << drag << endl;
	      Pout << "" << endl;
          }
      }

      // Drag force
      for(int j=0;j<3;j++) dragForce_[index][j] = drag[j];		  

    }

} 

void WenYuDragPost::setForce
(
    const volVectorField& gradPg_,
    const int & partNP_,
    labelList& partCellIDs_,
    List<vector>& partPositions_,
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

       for(int idir=0; idir<3; idir++) position[idir] = partPositions_[index][idir];
       
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

           if(verbose_ && index ==0 )
           {
               Pout << tab << "Gas pressure gradient " << endl; 
	       Pout << tab << "index = " << index << endl;
               Pout << tab << "gradP = " << gradP << endl;
	       Pout << "" << endl;
           }
	   
	   // Drag force
      	   for(int j=0;j<3;j++) gradPgForce_[index][j] = force[j];
           //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

       }
    }
}

void WenYuDragPost::setForceEulerian
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    List<scalar>& partRadii_,
    const volVectorField& Us_
) const
{
    // Resize arrays
    nP_ = U_.mesh().cells().size();
    dragForceEulerian_.resize(U_.mesh().cells().size());
    taupEulerian_.resize(U_.mesh().cells().size());
    positions_.resize(U_.mesh().cells().size());
    cellIDs_.resize(U_.mesh().cells().size());

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector Ur(0,0,0);
    vector drag(0,0,0);
    scalar voidfraction(0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar alps(0);
    
    forAll(U_.mesh().cells(),cellI)   	    
    {	       
      cellIDs_[cellI] = cellI;
      positions_[cellI] = U_.mesh().C()[cellI];
      voidfraction = voidfraction_[cellI];
      
      // Limit voidfraction
      if(voidfraction>1.00) voidfraction = 1.00;
      if(voidfraction<0.40) voidfraction = 0.40;
      
      Ur = U_[cellI]-Us_[cellI];
      ds = 2.*partRadii_[0]; // WARNING BE CAREFUL FOR BIDISPERSE CASE
      nuf = nufField[cellI];
      rho = rho_[cellI];
      magUr = mag(Ur);     
      alps = 1.-voidfraction; //+SMALL;

      if (magUr > 0)
      {
          // calc particle Re number
	  Rep = voidfraction*ds*magUr/(nuf+SMALL);

          // calc CD
	  if (Rep < 1000)
	  {
	    CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	  }else
	  {
            CD = 0.44*pow(voidfraction,-1.65); 
	  }  

	  //- Non linearity test
	  if(taupFuncPhi_) CD = 24./Rep*pow(voidfraction,-1.65);

          // Relaxation time
	  taupEulerian_[cellI] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);
	  
	  // Test correlation alppUr
	  if(taupConst_) taupEulerian_[cellI] = 0.025;
	  
	  // calc particle's drag
          drag = alps*rhoParticle_/taupEulerian_[cellI]*Ur;
	  
      }

      if(verbose_ && cellI <1 )
      {
	  Pout << tab << "WenYu Drag-ListRadii"  << endl;
	  Pout << tab << "index = " << cellI << endl;
          Pout << tab << "Up = " << Us_[cellI] << endl;
          Pout << tab << "Ur = " << Ur << endl;
          Pout << tab << "dp = " << ds << endl;
          Pout << tab << "rho = " << rho << endl;
          Pout << tab << "nuf = " << nuf << endl;
          Pout << tab << "voidfraction = " << voidfraction << endl;
          Pout << tab << "Rep = " << Rep << endl;
          Pout << tab << "alps = " << alps << endl;
          Pout << tab << "CD = " << CD << endl;
	  Pout << tab << "taupEulerian = " << taupEulerian_[cellI] << endl;
          Pout << tab << "drag = " << drag << endl;
	  Pout << "" << endl;
      }
      	
      // Drag force
      for(int j=0;j<3;j++) dragForceEulerian_[cellI][j] = drag[j];		  

    }

}

//-Set force Eulerian gas pressure term
void WenYuDragPost::setForceEulerian
(
    const volScalarField& voidfraction_, 
    const volVectorField& gradPg_	
) const
{
    // Resize arrays
    gradPgForceEulerian_.resize(voidfraction_.mesh().cells().size());
    
    vector force(0,0,0);

    forAll(voidfraction_.mesh().cells(),cellI)   	    
    {	   
    
      // calc particle's pressure gradient force
      force = -voidfraction_[cellI]*gradPg_[cellI];
      
      if(verbose_ && cellI == 0 )
      {          
          Pout << tab << "Gas pressure gradient"  << endl;
	  Pout << tab << "voidfraction = " << voidfraction_[cellI] << endl;
	  Pout << tab << "gradP = " << gradPg_[cellI] << endl;
	  Pout << "" << endl;
      }

      // Drag force
      for(int j=0;j<3;j++) gradPgForceEulerian_[cellI][j] = force[j];

    }
}

void WenYuDragPost::setForceParcel
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    parcelCloud& parcelCloud_
) const
{    
    // Info << " parcelCloud_.numberOfParticles() = " << parcelCloud_.numberOfParticles() << endl;
    // Resize arrays
    dragForce_.resize(parcelCloud_.numberOfParticles());
    taup_.resize(parcelCloud_.numberOfParticles());
    Ufxp_.resize(parcelCloud_.numberOfParticles());
    cellIDs_.resize(parcelCloud_.numberOfParticles());
    velocities_.resize(parcelCloud_.numberOfParticles());
    positions_.resize(parcelCloud_.numberOfParticles());    	
    voidfractionxp_.resize(parcelCloud_.numberOfParticles());

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = parcelCloud_.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = parcelCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI(-1);

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar mp(0);    
    scalar alps(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    //- Number of particle for drag autoPtr 
    nP_ = parcelCloud_.numberOfParticles();
       
    for(int index = 0;index <  parcelCloud_.numberOfParticles(); index++)
    {
       cellI = parcelCloud_.cellIDs(index);
       cellIDs_[index] = cellI; 
       if (cellI > -1) // particle Found
       {
	   position = parcelCloud_.position(index);  
	   positions_[index] = position;          
	   if(interpolation_)
           {	       
               voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
               Ufluid = UInterpolator_.interpolate(position,cellI);
               //Ensure interpolated void fraction to be meaningful
               // Info << " --> voidfraction: " << voidfraction << endl;
               if(voidfraction>1.00) voidfraction = 1.00;
               if(voidfraction<0.40) voidfraction = 0.40;
           }
	   else
           {
	       voidfraction = voidfraction_[cellI];
               Ufluid = U_[cellI];
           }

	   // Fluid velocity at particle position
	   Ufxp_[index] = Ufluid;
	   // Fluid volume fraction at particle position
	   voidfractionxp_[index] = voidfraction;

           Up = parcelCloud_.velocity(index);
	   velocities_[index] = Up;
           Ur = Ufluid-Up;
           ds = 2*parcelCloud_.radius(index);

	   nuf = nufField[cellI];
           rho = rho_[cellI];
           magUr = mag(Ur);
           Volp = ds*ds*ds*M_PI/6;
	   mp = rhoParticle_*Volp;
           alps = 1-voidfraction;

           if (magUr > 0)
           {
               // calc particle Re number
		 Rep = voidfraction*ds*magUr/(nuf+SMALL);

               // calc CD
	       if (Rep < 1000)
	       {
		 CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	       }else
	       {
        	 CD = 0.44*pow(voidfraction,-1.65); 
	       }  

	       //- Non linearity test
	       if(taupFuncPhi_) CD = 24./Rep*pow(voidfraction,-1.65); 

               // Relaxation time
	       taup_[index] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);

	       // Test correlation alppUr
	       if(taupConst_) taup_[index] = 0.025;
	      
	       // calc particle's drag
               drag = mp/taup_[index]*Ur;

	   }

           if(verbose_ && index <1 )
           {
	       Pout << tab << "WenYu Drag"  << endl;
	       Pout << tab << "index = " << index << endl;
               Pout << tab << "Up = " << Up << endl;
               Pout << tab << "Ur = " << Ur << endl;
               Pout << tab << "dp = " << ds << endl;
	       int nPIn = parcelCloud_.particlesInParcel(index).size();
	       Pout << tab << "nPIn = " << nPIn << endl;
               Pout << tab << "dparcel = " << pow(nPIn,1./3.)*ds << endl;
               Pout << tab << "rho = " << rho << endl;
               Pout << tab << "nuf = " << nuf << endl;
               Pout << tab << "voidfraction = " << voidfraction << endl;
               Pout << tab << "Rep = " << Rep << endl;
               Pout << tab << "Volp = " << Volp << endl;
               Pout << tab << "alps = " << alps << endl;
               Pout << tab << "CD = " << CD << endl;
	       Pout << tab << "taup = " << taup_[index] << endl;		   		    
               Pout << tab << "drag = " << drag << endl;
	       Pout << tab << "cell ID = " << cellI << endl;
	       Pout << "" << endl;
           }
       }

       // Drag force
       for(int j=0;j<3;j++) dragForce_[index][j] = drag[j];		  
	   
    }
    
} 

void WenYuDragPost::setForceParcel
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

//- New Functions using same weighting functions
//-
//-
//-
//-
//-
//-
//-
//-
//-
//-
//-
//-
//-
//- below functions I am using for mapping

// Overwrite function for List<vector>
void WenYuDragPost::setForce
(
    const volScalarField& voidfraction_,
    const volVectorField& U_,
    const int & partNP_,
    labelList& partCellIDs_,
    List<vector>& partPositions_,
    List<vector>& partVelocities_,
    List<scalar>& partRadii_,
    List<List<scalar> >& particleWeights_,
    labelListList& neighboringCellIDs_
) const
{

    // Resize arrays
    dragForce_.resize(partNP_);
    taup_.resize(partNP_);
    Ufxp_.resize(partNP_);
    cellIDs_.resize(partNP_);
    velocities_.resize(partNP_);
    positions_.resize(partNP_);	
    voidfractionxp_.resize(partNP_);
	
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = turbulence().mu()/rho_;
    	//const volScalarField nufField = refCast<const fvMesh>(obr_).turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = turbulence().nu();
	//const volScalarField& nufField = refCast<const fvMesh>(obr_).turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI(-1);
    label cellJ(-1);

    vector Up(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar CD(0);
    scalar Volp(0);
    scalar mp(0);    
    scalar alps(0);
    
    scalar weightP(0);
    scalar sumWeightP(0);
    vector voidfractionUfluid(0,0,0); 
    
    //- Number of particle for drag autoPtr 
    nP_ = partNP_;
    
    for(int index = 0;index <  partNP_; index++)
    {	   
      cellI = partCellIDs_[index];
      cellIDs_[index] = cellI; 

      sumWeightP = 0.;	
      if (cellI > -1) // particle Found
      {
          for(int idir=0; idir<3; idir++) position[idir] = partPositions_[index][idir];
	  positions_[index] = position;

  	  // Centering cell
	  weightP = particleWeights_[index][0];
	  sumWeightP += weightP;
          voidfraction = weightP*voidfraction_[cellI];
	  voidfractionUfluid = weightP*U_[cellI];; //weightP*voidfraction_[cellI]*U_[cellI];
	  
	  if(index==0) Info << tab << "index = " << index << " weightP = " << weightP;

	  // Neighboring cells    
	 
	  for(int subCell=1 ; subCell < maxCellPerPart_; subCell++)
	  {
              cellJ = neighboringCellIDs_[index][subCell];
              if (cellJ > -1)
              {
        	  weightP = particleWeights_[index][subCell];
		  sumWeightP += weightP;
		  
		  if(index==0) Info << tab << " weightP-"<< subCell <<" = " << weightP;
        	  voidfraction += weightP*voidfraction_[cellJ];
	          voidfractionUfluid += weightP*U_[cellJ]; //weightP*voidfraction_[cellJ]*U_[cellJ];
              }
	  }
	  if(index==0) Info << endl;
	  
	  // Normalize Uf & voidfraction
	  Ufluid = voidfractionUfluid/sumWeightP; //voidfractionUfluid/voidfraction;
          voidfraction = voidfraction/sumWeightP;

	  // Limit
	  if(voidfraction>1.00) voidfraction = 1.00;
          if(voidfraction<0.40) voidfraction = 0.40;
	   
	  // Fluid velocity at particle position
	  Ufxp_[index] = Ufluid;
	  

	  
	  // Fluid volume fraction at particle position
	  voidfractionxp_[index] = voidfraction;

          for(int idir=0; idir<3; idir++) Up[idir] = partVelocities_[index][idir];
	  velocities_[index] = Up;
          Ur = Ufluid-Up;
          ds = 2.*partRadii_[index];

	  nuf = nufField[cellI];
          rho = rho_[cellI];
          magUr = mag(Ur);
          Volp = ds*ds*ds*M_PI/6;
	  mp = Volp*rhoParticle_;
          alps = 1-voidfraction;

          if (magUr > 0)
          {
              // calc particle Re number
		Rep = voidfraction*ds*magUr/(nuf+SMALL);

              // calc CD
	      if (Rep < 1000)
	      {
		CD = 24./Rep*(1.+0.15*pow(Rep,0.687))*pow(voidfraction,-1.65);
	      }else
	      {
        	CD = 0.44*pow(voidfraction,-1.65); 
	      }  

	      //- Non linearity test
	      if(taupFuncPhi_) CD = 24./Rep*pow(voidfraction,-1.65);
	       
              // Relaxation time
	      taup_[index] = 4./3.*rhoParticle_/rho*ds/(CD*magUr);

	      // Test correlation alppUr
	      if(taupConst_) taup_[index] = 0.025;
	      
	      // calc particle's drag
              drag = mp/taup_[index]*Ur;

	  }

          if(verbose_ && index <1 )
          {
	      Pout << tab << "WenYu Drag-Lists"  << endl;
	      Pout << tab << "index = " << index << endl;
	      Pout << tab << "cellI = " << cellI << endl;
              Pout << tab << "Up = " << Up << endl;
              Pout << tab << "Ur = " << Ur << endl;
              Pout << tab << "dp = " << ds << endl;
              Pout << tab << "rho = " << rho << endl;
              Pout << tab << "nuf = " << nuf << endl;
              Pout << tab << "voidfraction = " << voidfraction << endl;
              Pout << tab << "Rep = " << Rep << endl;
              Pout << tab << "Volp = " << Volp << endl;
              Pout << tab << "alps = " << alps << endl;
              Pout << tab << "CD = " << CD << endl;
	      Pout << tab << "taup = " << taup_[index] << endl;		   		    
              Pout << tab << "drag = " << drag << endl;
              Pout << tab << "cell ID = " << cellI << endl;
	      Pout << "" << endl;
          }
      }

      // Drag force
      for(int j=0;j<3;j++) dragForce_[index][j] = drag[j];		  

    }

} 

void WenYuDragPost::setForce
(
    const volVectorField& gradPg_,
    const int & partNP_,
    labelList& partCellIDs_,
    List<vector>& partPositions_,
    scalarList& partRadii_,
    List<List<scalar> >& particleWeights_,
    labelListList& neighboringCellIDs_	
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
    label cellJ(-1);
    
    scalar weightP(0.);
    scalar sumWeightP(0.);

    for(int index = 0;index < partNP_; index++)
    {
       cellI = partCellIDs_[index];
       
       for(int idir=0; idir<3; idir++) position[idir] = partPositions_[index][idir];
       
       sumWeightP = 0.;	
       if (cellI > -1) // particle Found
       {

  	   // Centering cell
	   weightP = particleWeights_[index][0];
	   sumWeightP += weightP;
           gradP = weightP*gradPg_[cellI];

	   // Neighboring cells    
	   for(int subCell=1 ; subCell < maxCellPerPart_; subCell++)
	   {
               cellJ = neighboringCellIDs_[index][subCell];
               if (cellJ > -1)
               {
        	   weightP = particleWeights_[index][subCell];
	   	   sumWeightP += weightP;
        	   gradP += weightP*gradPg_[cellJ];
               }
	   }

           gradP /= sumWeightP; 
           ds = partRadii_[index];
           Vs = ds*ds*ds*M_PI/6;

           // calc particle's pressure gradient force
           force = -Vs*gradP;

           if(verbose_ && index ==0 )
           {
               Pout << tab << "Gas pressure gradient " << endl; 
	       Pout << tab << "index = " << index << endl;
               Pout << tab << "gradP = " << gradP << endl;
	       Pout << "" << endl;
           }
	   
	   // Drag force
      	   for(int j=0;j<3;j++) gradPgForce_[index][j] = force[j];
           //for(int j=0;j<3;j++) force_[index][j] += force[j];		  

       }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
