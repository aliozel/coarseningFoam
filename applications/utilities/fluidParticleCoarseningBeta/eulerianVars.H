//-Eulerian variables
/*
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
*/
//-Eulerian variables
volVectorField normalizedEulerianDragCoeff
(
    IOobject
    (
        "normalizedEulerianDragCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedEulerianDragCoeffNum
(
    IOobject
    (
        "normalizedEulerianDragCoeffNum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);
volVectorField normalizedEulerianDragCoeffDom
(
    IOobject
    (
        "normalizedEulerianDragCoeffDom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "zero", dimensionSet(1,-2,-2,0,0), vector(0,0,0) )
);

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
