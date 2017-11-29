// Read gas velocity 
Info<< tab <<"Reading U" << endl;
//U.read();
IOobject UHeader
(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volVectorField U(UHeader,mesh);

// Read particle velocity 
Info<< tab <<"Reading Us" << endl;
IOobject UsHeader
(
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volVectorField Us(UsHeader,mesh);

// Read gas volume fraction
Info<< tab <<"Reading gas volume fraction" << endl;
IOobject voidfractionHeader
(
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField voidfraction(voidfractionHeader,mesh);
