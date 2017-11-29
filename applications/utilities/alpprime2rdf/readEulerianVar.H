// Read volume fraction of gas
IOobject colorFuncHeader
(
	"voidfraction",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ
);

Info<< "        Reading voidfraction" << endl;
volScalarField voidfraction(colorFuncHeader,mesh);

volScalarField alpha = 1.-voidfraction;
alphaMean = fvc::domainIntegrate(alpha*mesh.V()).value()/fvc::domainIntegrate(mesh.V()).value();
Info << tab << "Domain-average solid vol. frac = " << alphaMean << endl;
