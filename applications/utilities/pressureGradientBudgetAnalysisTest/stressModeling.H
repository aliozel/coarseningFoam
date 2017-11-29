// Total stress
volSymmTensorField totalParticleStress
(
    IOobject
    (
	"totalParticleStress",
	mesh.time().timeName(),
	mesh,
	IOobject::NO_READ, //READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
);
totalParticleStress = sigmaKin - sigmaColl;
