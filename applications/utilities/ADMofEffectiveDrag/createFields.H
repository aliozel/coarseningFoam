volScalarField alpp
(
    IOobject
    (
        "alpp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);

volScalarField barAlppAlpp
(
    IOobject
    (
        "barAlppAlpp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), scalar(0))
);



