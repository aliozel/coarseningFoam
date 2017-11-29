    Info<< "Reading field p\n" << endl;
    /*
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
    */ 
    
    Info<< "Reading gas velocity field U" << endl;
    //Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
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
     
//===============================
// particle interaction modelling
//===============================
    /*
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
        //dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );
    */
    /*
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
    */
    /*
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
        //dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );
    */
    /*
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
    */
    
    /*
    Info<< "\nExplicit 0term fzero" << endl;
    volVectorField fzero
    (
        IOobject
        (
            "fzero",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );
    */
    
    /*  
    Info<< "\nResidual 0term fresidual" << endl;
    volVectorField fresidual
    (
        IOobject
        (
            "fresidual",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedVector("residual", dimensionSet(1,-2,-2,0,0),vector(0,0,0)) 
    );
    */    
//===============================

//#   include "createPhi.H"
//#ifndef createPhi_H
//#define createPhi_H
/*
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
 */
//#endif

    /*
    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
    
    scalar rhoParticle(0);	
    bool fullyPeriodicFLow(false);

    bool mixtureVelocityCorrectionPeriodicFlow(false);
    if (pimpleDict.found("mixtureVelocityCorrectionPeriodicFlow"))
    {
        mixtureVelocityCorrectionPeriodicFlow = true;
        Info << "Mixture velocity correction is activated for periodic flow " << endl;
    } 
    
    // Boundary patches, find coupled patches
    label nCoupledPatches(0);
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].type() == "cyclic" )
        {
        	nCoupledPatches++;
	}
    }

    if ( nCoupledPatches == 6 ) 
    {
	Info << "Number of coupled patches = " << nCoupledPatches << endl;
        fullyPeriodicFLow = true; 
	Info << "Fully periodic flow " << endl;
	rhoParticle = readScalar(pimpleDict.lookup("rhoParticle"));
    }	

    bool correctCouplingError(false);
    if (pimpleDict.found("correctCouplingError"))
    {
        correctCouplingError = true;
	Info << "Momentum coupling error correction is active " << endl;
    }
   
    if ( nCoupledPatches < 6 )
    {
	mixtureVelocityCorrectionPeriodicFlow = false;
    }	
     
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
    */

    // Filtered variables 
    volScalarField barAlpp
    (
        IOobject
        (
            "barAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ, //MUST_READ,
            IOobject::AUTO_WRITE
        ),
        //1-voidfraction //mesh
	mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(0))	
    );

    volScalarField barAlppAlpp
    (
        IOobject
        (
            "barAlppAlpp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ, //MUST_READ,
            IOobject::AUTO_WRITE
        ),
        //1-voidfraction //mesh
	mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(0))	
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
        //(1.-voidfraction)*U 
	mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0))
    );

    volVectorField tildetildeUs
    (
        IOobject
        (
            "tildetildeUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //(1.-voidfraction)*U 
	mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0))
    );

    volSymmTensorField tildeUsUs
    (
        IOobject
        (
            "tildeUsUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //(1.-voidfraction)*U 
	mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
    );

   

    volSymmTensorField tildetUstUs
    (
        IOobject
        (
            "tildetUstUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //(1.-voidfraction)*U 
	mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
    );

    volVectorField tildeUsPrime
    (
        IOobject
        (
            "tildeUsPrime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //(1.-voidfraction)*U 
	mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0))
    );

    volSymmTensorField diffTildeeAlppUsUs
    (
        IOobject
        (
            "diffTildeeAlppUsUs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
    );

    

    volSymmTensorField sigmaColl
    (
        IOobject
        (
            "sigmaColl",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedSymmTensor("zero",dimensionSet(1,-1,-2,0,0),symmTensor(0,0,0,0,0,0))
    );

    /*  
    volTensorField sigmaColl
    (
        IOobject
        (
            "sigmaColl",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedTensor("zero",dimensionSet(1,-1,-2,0,0),tensor(0,0,0,0,0,0,0,0,0))
    );
    */

   // Mean velocities in each cell
   volVectorField meanVelCell
   (
        IOobject
        (
           "meanUp",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
         ),
       mesh,
       dimensionedVector( "zero", dimensionSet(0,1,-1,0,0), vector(0,0,0) )
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
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(1.)) 
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
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0)) 
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
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0),scalar(0)) 
    );

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
        //(1.-voidfraction)*U 
	mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0),vector(0,0,0))
    );

    volSymmTensorField tildetUstUsHatStar
    (
        IOobject
        (
            "tildetUstUsHatStar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //(1.-voidfraction)*U 
	mesh,
	dimensionedSymmTensor("zero",dimensionSet(0,2,-2,0,0),symmTensor(0,0,0,0,0,0))
    );

    volScalarField mappedCoordNumber
    (
        IOobject
        (
            "mappedCoordNumber",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,1,-1,0,0),scalar(0))
    );

