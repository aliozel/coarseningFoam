//
// - Filtered variables
//
    volVectorField tildeUf
    (
        IOobject
        (
            "tildeUf"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField barAlpf
    (
        IOobject
        (
            "barAlpf"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField barAlpp
    (
        IOobject
        (
            "barAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField barAlppUf
    (
        IOobject
        (
            "barAlppUf"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField barAlppAlpp
    (
        IOobject
        (
            "barAlppAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh  
    );

    volVectorField tildeUs
    (
        IOobject
        (
            "tildeUs"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
    );

    volVectorField tildeAlppUr
    (
        IOobject
        (
            "tildeAlppUr"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
    );
    
    // Create filtered Eulerian force	
    volVectorField barForceEulerian
    (
	IOobject
	(
            "barForceEulerian"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	mesh
    );

    volVectorField tildeGradPg
    (
	IOobject
	(
            "tildeGradPg"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	mesh
    );

    //- barAlppFluc2
    volScalarField barFlucAlppFlucAlpp
    (
        IOobject
        (
            "barFlucAlppFlucAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh  
    );

    volScalarField barAlpp
    (
        IOobject
        (
            "barAlpp"+filterWidthName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
