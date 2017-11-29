// Collision parameters
word collisionModelI("hertz");

scalarField youngsModulus(nParticleClass,scalar(0));
scalarField poissonsRatio(nParticleClass,scalar(0));

vectorField coefficientRestitution(nParticleClass,vector(0,0,0));
vectorField coefficientFriction(nParticleClass,vector(0,0,0));

scalar k_n(0);
scalar k_t(0);
scalar gamma_n(0); 
scalar gamma_t(0);
scalar mu_f(0);
scalar e_n(0);
scalar e_t(0); 
scalar meff(0);

scalar dt(0);

bool tangential_history=false;
if( calcCollision )
{
	if(particleStressPropsDict.found("tangential_history")) 	
	{
		tangential_history = true;
		Info << "Tangential history effect is active " << endl;
	}
}

bool liquid_transfer=false;
if(particleStressPropsDict.found("liquid_transfer")) liquid_transfer = true;	
scalar surf_tension(0);
scalar fluid_visc(0);

bool cohesion=false;
if(particleStressPropsDict.found("cohesion")) cohesion = true;	
scalar minimumDistanceVdW(0);
vectorField cohEnergyDens(nParticleClass,vector(0,0,0));

// Liquid bridge form
bool first_touch(false); 

