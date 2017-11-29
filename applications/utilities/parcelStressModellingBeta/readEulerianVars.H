// Read gas volume fraction
/* 
IOobject voidfractionmHeader
(
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
);
Info<< " 	Reading voidfraction" << endl;
volScalarField voidfractionm(voidfractionmHeader,mesh);
 
// Read Eulerian particle velocity
IOobject Usmheader
(
   "Us",
   runTime.timeName(),
   mesh,
   IOobject::MUST_READ	
);

Info<< " 	Reading Us" << endl;
volVectorField Usm(Usmheader,mesh);
*/

// Read gas velocity
IOobject Umheader
(
   "U",
   runTime.timeName(),
   mesh,
   IOobject::MUST_READ	
);

Info<< " 	Reading U" << endl;
volVectorField Um(Umheader,mesh);

// Read particle collisional stresses
IOobject sigmaCollmHeader
(
    "sigmaColl",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);

Info<< " 	Reading sigmaColl" << endl;
volSymmTensorField sigmaCollm(sigmaCollmHeader,mesh);

//voidfraction = voidfractionm;
//Us = Usm;
//sigmaColl = sigmaCollm;
U=Um;

forAll(mesh.cells(),cellI)
{
   sigmaColl[cellI].xx() = sigmaCollm[cellI][0];
   sigmaColl[cellI].yy() = sigmaCollm[cellI][1];
   sigmaColl[cellI].zz() = sigmaCollm[cellI][2];
   
   sigmaColl[cellI].xy() = sigmaCollm[cellI][3];
   sigmaColl[cellI].xz() = sigmaCollm[cellI][4];      
   sigmaColl[cellI].yz() = sigmaCollm[cellI][5];         
}	

