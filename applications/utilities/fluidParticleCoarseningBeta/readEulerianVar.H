// Read gas velocity 
Info<< tab <<"Reading U" << endl;
//U.read();
IOobject UmHeader
(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volVectorField Um(UmHeader,mesh);

// Read gas volume fraction
Info<< tab <<"Reading gas volume fraction" << endl;
//voidfraction.read();
IOobject voidfractionmHeader
(
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField voidfractionm(voidfractionmHeader,mesh);

// Read gas pressure gradient
Info<< tab <<"Reading Pg" << endl;
//p.read();
IOobject pmHeader
(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField pm(pmHeader,mesh);

//- Pressure gradient
//volVectorField gradPg(fvc::grad(p));

// Read gas pressure gradient
Info<< tab <<"Reading gas density" << endl;
//rho.read();
IOobject rhomHeader
(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField rhom(rhomHeader,mesh);

// Read interfacial mom. transfer coefficient
Info<< tab <<"Reading Ksl" << endl;
IOobject KslmHeader
(
        "Ksl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
volScalarField Kslm(KslmHeader,mesh);

U = Um;
voidfraction = voidfractionm;
p = pm;
//Calculate gas pressure
//volVectorField gradPg(fvc::grad(p));
gradPg = fvc::grad(p);
rho = rhom;
Ksl = Kslm;


