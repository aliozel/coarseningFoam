// Init 
scalarField liq(particlesInSubDomain.size());

vectorField fcoll(particlesInSubDomain.size(),vector(0,0,0));
vectorField ftan(particlesInSubDomain.size(),vector(0,0,0));
vectorField fcap(particlesInSubDomain.size(),vector(0,0,0));
vectorField fvisc(particlesInSubDomain.size(),vector(0,0,0));
vectorField fcoh(particlesInSubDomain.size(),vector(0,0,0));
vectorField liquidVol(particlesInSubDomain.size(),vector(0,0,0));

// Liquid bridge form
bool first_touch(false); 

