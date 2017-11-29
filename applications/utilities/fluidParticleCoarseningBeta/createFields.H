    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ); 

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "\nCreating dummy density field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
        
    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
     (
	 IOobject
	 (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
	 ),
	 linearInterpolate(U) & mesh.Sf()
     );

//
// - Filtered variables
//
    volVectorField tildeUf
    (
        IOobject
        (
            "tildeUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );
    
    volScalarField barAlpf
    (
        IOobject
        (
            "barAlpf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
    );

    volScalarField barAlpp
    (
        IOobject
        (
            "barAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
    );

    volVectorField barAlppUf
    (
        IOobject
        (
            "barAlppUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );

    volScalarField barAlppAlpp
    (
        IOobject
        (
            "barAlppAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))   
    );

    volVectorField tildeUs
    (
        IOobject
        (
            "tildeUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );

    volVectorField tildeAlppUr
    (
        IOobject
        (
            "tildeAlppUr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0), vector(0,0,0))
    );

    volVectorField dragFpEulerian
    (
	IOobject
	(
            "dragFpEulerian",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );    

    // Create filtered Eulerian force	
    volVectorField barForceEulerian
    (
	IOobject
	(
            "barForceEulerian",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
        dimensionedVector("0", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );

    volVectorField gradPg
    (
	IOobject
	(
            "gradPg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
        dimensionedVector("0", dimensionSet(0,1,-2,0,0),vector(0,0,0)) 
    );

    volVectorField tildeGradPg
    (
	IOobject
	(
            "tildeGradPg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh,
        dimensionedVector("0", dimensionSet(0,1,-2,0,0),vector(0,0,0)) 
    );

    volScalarField mappedVoidfraction
    (
        IOobject
        (
            "mappedVoidfraction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0,0,0,0,0),scalar(1.)) 
    );
    
    volVectorField mappedUs
    (
        IOobject
        (
            "mappedUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(0,1,-1,0,0),vector(0,0,0)) 
    );

    volScalarField weightField
    (
        IOobject
        (
            "weightField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0,0,0,0,0),scalar(0)) 
    );

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    //- alppFluc2
    volScalarField flucAlppFlucAlpp
    (
        IOobject
        (
            "flucAlppFlucAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))   
    );

    //- barAlppFluc2
    volScalarField barFlucAlppFlucAlpp
    (
        IOobject
        (
            "barFlucAlppFlucAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))   
    );

    //- gradBarAlppGradBarAlpp
    volScalarField gradBarAlppGradBarAlpp
    (
        IOobject
        (
            "gradBarAlppGradBarAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0, -2, 0, 0, 0), scalar(0))   
    );    
