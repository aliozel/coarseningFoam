// Test

  // Calculate Eulerian resolved drag force
Info << "\n" << endl;
Info << tab << "Calculating Eulerian resolved drag force..." << endl;
resolvedDragFEulerian().setForceEulerian(barAlpf,tildeUf,radii[me],tildeUs);

volScalarField OneAlpf = 1/mappedVoidfraction;
volScalarField barOneAlpf = convKernel.coarseningSca("barOneAlpf",filterWidth[fWidth],OneAlpf);
volScalarField OneAlpfPrime = OneAlpf - barOneAlpf;
volVectorField FPrime = dragFpEulerian - barForceEulerian;

//Y1
volVectorField barOneAlpfF = convKernel.coarseningVec("barOneAlpfF",filterWidth[fWidth],OneAlpf*dragFpEulerian);
volVectorField Y1
(
    IOobject
    (
        "Y1"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barOneAlpfF
);
Y1.write();

//Y2
volVectorField barOneAlpfBarF = barOneAlpf*barForceEulerian;
volVectorField Y2
(
    IOobject
    (
        "Y2"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barOneAlpfBarF
);
Y2.write();

//Y3
volVectorField gradBarOneAlpf = fvc::grad(barOneAlpf);
volTensorField gradBarF = fvc::grad(barForceEulerian);
volVectorField gradBarOneAlpfGradBarF = gradBarF&gradBarOneAlpf;
volVectorField Y3
(
    IOobject
    (
        "Y3"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradBarOneAlpfGradBarF
);
Y3.write();  

//Y4
volVectorField laplacianBarF = fvc::laplacian(barForceEulerian);
volVectorField barOneAlpfLaplacianBarF = barOneAlpf*laplacianBarF;
volVectorField Y4
(
    IOobject
    (
        "Y4"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barOneAlpfLaplacianBarF
);
Y4.write();  

//Y5
volScalarField laplacianBarOneAlpf = fvc::laplacian(barOneAlpf);
volVectorField laplacianBarOneAlpfBarF = laplacianBarOneAlpf*barForceEulerian;
volVectorField Y5
(
    IOobject
    (
        "Y5"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    laplacianBarOneAlpfBarF
);
Y5.write();

//Y6
volVectorField barOneAlpfPrimeFPrime = convKernel.coarseningVec("barOneAlpfPrimeFPrime",filterWidth[fWidth],OneAlpfPrime*FPrime);
volVectorField Y6
(
    IOobject
    (
        "Y6"+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    barOneAlpfPrimeFPrime
);
Y6.write();
