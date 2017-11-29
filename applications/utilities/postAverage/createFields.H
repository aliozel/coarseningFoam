    Info<< "Reading field U\n" << endl;
    volVectorField particleEulerVelocity
    (
        IOobject
        (
            "particleEulerVelocity",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );
