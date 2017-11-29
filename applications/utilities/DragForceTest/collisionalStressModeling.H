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

forAll(mesh.cells(),cellI)
{
   sigmaColl[cellI].xx() = sigmaCollm[cellI][0];
   sigmaColl[cellI].yy() = sigmaCollm[cellI][1];
   sigmaColl[cellI].zz() = sigmaCollm[cellI][2];
   
   sigmaColl[cellI].xy() = sigmaCollm[cellI][3];
   sigmaColl[cellI].xz() = sigmaCollm[cellI][4];      
   sigmaColl[cellI].yz() = sigmaCollm[cellI][5];         
}	

