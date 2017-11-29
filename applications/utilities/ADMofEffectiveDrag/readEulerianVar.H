// Read voidfraction
Info<< tab <<"Reading voidfraction" << endl;
IOobject mappedVoidfractionmHeader
(
        "mappedVoidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField mappedVoidfractionm(mappedVoidfractionmHeader,mesh);

// Read gas velocity 
Info<< tab <<"Reading U" << endl;
IOobject UmHeader
(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volVectorField Um(UmHeader,mesh);

// Read solid velocity 
Info<< tab <<"Reading Us" << endl;
IOobject UsmHeader
(
        "mappedUs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volVectorField Usm(UsmHeader,mesh);

// End of reading
volScalarField alpp = 1.-mappedVoidfractionm;
volVectorField U = Um;
volVectorField Us = Usm;
