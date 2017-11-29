// Define "Eulerian" coordination number in a cell {numberOfContacts|_{cell}/nP|_{cell}}		
volScalarField coordNumberEul
(
    IOobject
    (
        "coordNumberEul",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0,0,0,0,0),scalar(0))
);

// Define "Eulerian" coordination number in a cell {numberOfContacts|_{cell}/(nP|_{cell}-nP|{without contacts})}		
volScalarField coordNumberEulWithoutZeroContact
(
    IOobject
    (
        "coordNumberEulWithoutZeroContact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0,0,0,0,0),scalar(0))
);

// Define "Eulerian" coordination number in a cell {numberOfContacts|_{cell}/nP|_{cell}}		
volScalarField numberOfParticlesInCell
(
    IOobject
    (
        "numberOfParticlesInCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0,0,0,0,0),scalar(0))
);

// Define "Eulerian" coordination number in a cell {numberOfContacts|_{cell}/(nP|_{cell}-nP|{without contacts})}		
volScalarField numberOfParticlesInCellWithoutColl
(
    IOobject
    (
        "numberOfParticlesInCellWithoutColl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
//        IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0,0,0,0,0),scalar(0))
);
