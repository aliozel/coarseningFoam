// Uf@p weighted baralpp
interpolationCellPoint<vector> barAlppUfxp(barAlppUf);
// tildeUs@xp
interpolationCellPoint<vector> tildeUsxp(tildeUs);

//vectorField uSlipLag(nP[me],vector(0,0,0));				  
vectorField gij(nP[me],vector(0,0,0));
vectorField alppVdi(nP[me],vector(0,0,0));
//scalarField unresolvedDragUSlip(nP[me],0);

// Unresolved part of particle drag force
//vectorField normalizedParticleDragCoeff(nP[me],vector(0,0,0));
IOobject normalizedParticleDragCoeffIO
(
	"normalizedParticleDragCoeff"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedParticleDragCoeff(normalizedParticleDragCoeffIO,nP[me]);

// Unresolved part of particle drag force numerator
IOobject normalizedParticleDragCoeffNumIO
(
	"normalizedParticleDragCoeffNum"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedParticleDragCoeffNum(normalizedParticleDragCoeffNumIO,nP[me]);

// Unresolved part of particle drag force denominator
IOobject normalizedParticleDragCoeffDomIO
(
	"normalizedParticleDragCoeffDom"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedParticleDragCoeffDom(normalizedParticleDragCoeffDomIO,nP[me]);

IOobject unresolvedDragIO
(
	"unresolvedDrag"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> unresolvedDrag(unresolvedDragIO,nP[me]);

IOobject uSlipLagIO
(
	"uSlipLag"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> uSlipLag(uSlipLagIO,nP[me]);

IOobject resolvedUSlipLagIO
(
	"resolvedUSlipLag"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> resolvedUSlipLag(resolvedUSlipLagIO,nP[me]);

IOobject gradBarAlppxpIO
(
	"gradBarAlppxp"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> gradBarAlppxp(gradBarAlppxpIO,nP[me]);
interpolationCellPoint<vector> IntGradBarAlppxp(gradBarAlpp);

IOobject resolvedTaupIO
(
	"resolvedTaup"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<scalar> resolvedTaup(resolvedTaupIO,nP[me]);

IOobject diffResolvedUUsxpIO
(
	"diffResolvedUUsxp"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> diffResolvedUUsxp(diffResolvedUUsxpIO,nP[me]);

/* IOobject gradTildeUfzxpIO
(
	"gradTildeUfzxp"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> gradTildeUfzxp(gradTildeUfzxpIO,nP[me]);
interpolationCellPoint<vector> IntGradTildeUfzxp(gradTildeUfz);	 */	    

IOobject barForceEulerianxpIO
(
	"barForceEulerianxp"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> barForceEulerianxp(barForceEulerianxpIO,nP[me]);
interpolationCellPoint<vector> IntbarForceEulerianxp(barForceEulerian);	

IOobject resolvedLagDragIO
(
	"resolvedLagDrag"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> resolvedLagDrag(resolvedLagDragIO,nP[me]);


//-Eulerian variables
IOobject normalizedEulerianDragCoeffIO
(
	"normalizedEulerianDragCoeff"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedEulerianDragCoeff(normalizedEulerianDragCoeffIO,mesh.cells().size());

IOobject normalizedEulerianDragCoeffNumIO
(
	"normalizedEulerianDragCoeffNum"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedEulerianDragCoeffNum(normalizedEulerianDragCoeffNumIO,mesh.cells().size());

IOobject normalizedEulerianDragCoeffDomIO
(
	"normalizedEulerianDragCoeffDom"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<vector> normalizedEulerianDragCoeffDom(normalizedEulerianDragCoeffDomIO,mesh.cells().size());

IOobject baralppxpIO
(
	"baralppxp"+filterWidthName,
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE,
	false
);		    		    
IOField<scalar> baralppxp(baralppxpIO,nP[me]);
