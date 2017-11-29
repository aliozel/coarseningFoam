// Particle pressure
volScalarField Pp
(
    IOobject
    (
        "Pp",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

// Liggths definition (sigmaColl: -r:F_cont)
//Pp = 1./3. * tr( tildeSigmaKin + sigmaColl ) ;
Pp = 1./3. * tr( tildeSigmaKin - sigmaColl ) ;
// Write into the results folder
Pp.write();

dimensionedScalar PpLogsmall("zero", dimensionSet(1,-1,-2,0,0), 1e-32);
dimensionedScalar PpONE("zero", dimensionSet(1,-1,-2,0,0), 1.);
// Particle pressure
volScalarField PpLog
(
    IOobject
    (
        "PpLog",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
 );

// Liggths definition (sigmaColl: -r:F_cont)
//Pp = 1./3. * tr( tildeSigmaKin + sigmaColl ) ;
PpLog = log10(max(Pp,PpLogsmall)/PpONE) ;
// Write into the results folder
PpLog.write();

// Particle pressure
volScalarField Ppcontact
(
    IOobject
    (
        "Ppcontact",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

// Liggths definition (sigmaColl: -r:F_cont)
//Pp = 1./3. * tr( tildeSigmaKin + sigmaColl ) ;
Ppcontact = 1./3. * tr( - sigmaColl ) ;
// Write into the results folder
Ppcontact.write();

// Particle pressure
volScalarField Ppkinetic
(
    IOobject
    (
        "Ppkinetic",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

// Liggths definition (sigmaColl: -r:F_cont)
//Pp = 1./3. * tr( tildeSigmaKin + sigmaColl ) ;
Ppkinetic = 1./3. * tr( tildeSigmaKin ) ;
// Write into the results folder
Ppkinetic.write();

volScalarField Ppxy
(
    IOobject
    (
        "Ppxy",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

//Ppxy = 1./2. *(tildeSigmaKin.component(tensor::XX) + sigmaColl.component(tensor::XX) +tildeSigmaKin.component(tensor::YY) + sigmaColl.component(tensor::YY));
forAll(mesh.cells(),cellI)
{
  Ppxy[cellI] = 1./2. *(   tildeSigmaKin[cellI].xx() - sigmaColl[cellI].xx()
                         + tildeSigmaKin[cellI].yy() - sigmaColl[cellI].yy() );
}			 

// Write into the results folder
Ppxy.write();


volScalarField Ppxx
(
    IOobject
    (
        "Ppxx",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

volScalarField Ppyy
(
    IOobject
    (
        "Ppyy",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

volScalarField Ppzz
(
    IOobject
    (
        "Ppzz",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );




volScalarField PpxxPp
(
    IOobject
    (
        "PpxxPp",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
 );


forAll(mesh.cells(),cellI)
{
  Ppxx[cellI] = tildeSigmaKin[cellI].xx() - sigmaColl[cellI].xx();
                    
}	
Ppxx.write();		 
PpxxPp = Ppxx / max(Pp,PpLogsmall);
// Write into the results folder
PpxxPp.write();

volScalarField PpyyPp
(
    IOobject
    (
        "PpyyPp",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
 );


forAll(mesh.cells(),cellI)
{
  Ppyy[cellI] = tildeSigmaKin[cellI].yy() - sigmaColl[cellI].yy();
                    
}	
Ppyy.write();		 		 
PpyyPp = Ppyy / max(Pp,PpLogsmall);
// Write into the results folder
PpyyPp.write();

volScalarField PpzzPp
(
    IOobject
    (
        "PpzzPp",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
 );


forAll(mesh.cells(),cellI)
{
  Ppzz[cellI] = tildeSigmaKin[cellI].zz() - sigmaColl[cellI].zz();
                    
}	
Ppzz.write();		 		 
PpzzPp = Ppzz / max(Pp,PpLogsmall);
// Write into the results folder
PpzzPp.write();

/*
volScalarField Ppzz
(
    IOobject
    (
        "Ppzz",//+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
 );

Ppzz = tildeSigmaKin.component(tensor::ZZ) + sigmaColl.component(tensor::ZZ) ;


// Write into the results folder
Ppzz.write();
*/
// Calculate the particulate pressure gradient
volVectorField gradPp(fvc::grad(Pp));

// Particle shear stress
volTensorField sigmap
(
    IOobject
    (
        "sigmap", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
);

//sigmap = ( tildeSigmaKin + sigmaColl ) - tensor(I) * Pp;
//sigmap = ( tildeSigmaKin - sigmaColl ) - tensor(I) * Pp;
sigmap = tensor(I) * Pp - ( tildeSigmaKin - sigmaColl );
// Write into the results folder
sigmap.write();

// Particle viscosity
volScalarField mup
(
    IOobject
    (
        "mup", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

// Particle shear stresses
volTensorField S("S",fvc::grad(mappedUs) + fvc::grad(mappedUs)().T()-2./3.*tr(fvc::grad(mappedUs))*tensor(I));

dimensionedScalar SSsmall("zero", dimensionSet(0,0,-2,0,0,0,0), SMALL);

//mup = ( sigmap && S ) / ( max ( S && S, SSsmall ) );
mup = ( sigmap && S ) / ( max ( S && S, SSsmall ) );
// Limit by zero
//mup.max(0.);
// Write into the results folder
mup.write(); 

// Particle viscosity
volScalarField mupVS
(
    IOobject
    (
        "mupVS", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);


dimensionedScalar SSsmallVS("zero", dimensionSet(0,0,-2,0,0,0,0), 1.0e-37);

//mup = ( sigmap && S ) / ( max ( S && S, SSsmall ) );
mupVS = ( sigmap && S ) / ( max ( S && S, SSsmallVS ) );
mupVS.write(); 

// Particle viscosity/sqrt(p)
dimensionedScalar PpSmall("zero", dimensionSet(1,-1,-2,0,0,0,0), 1.e-06);
volScalarField mupSqrtPp
(
    IOobject
    (
        "mupSqrtPp", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mup/sqrt((mag(Pp)+PpSmall)*rhoParticle)/dp
); 

// Pp/(rhoParticle*(S:S)*dp^2)
volScalarField normPp
(
    IOobject
    (
        "normPp", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Pp/(rhoParticle*( (S && S) + SSsmall )*dp*dp)
); 

// Pp/(rhoParticle*(S:S)*dp^2)
dimensionedScalar mupSmall("zero", dimensionSet(1,-1,-2,0,0,0,0), 1.e-06);
volScalarField normmup
(
    IOobject
    (
        "normmup", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mup/(rhoParticle*sqrt( (S && S) + SSsmall )*dp*dp)
); 

///ADDITION BY YILE GU
//////////////////////////////////////////////////////
// Particle hat shear stress
volTensorField sigmapHAT
(
    IOobject
    (
        "sigmapHAT", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
);

sigmapHAT = sigmap;

volTensorField SHAT("SHAT",S);

forAll(mesh.cells(),cellI)
{
    sigmapHAT[cellI].xx() = 0.;
    SHAT[cellI].xx() = 0.;
    sigmapHAT[cellI].yy() = 0.;
    SHAT[cellI].yy() = 0.;
    sigmapHAT[cellI].zz() = 0.;
    SHAT[cellI].zz() = 0.;
}

dimensionedScalar mupLogsmall("zero", dimensionSet(1,-1,-1,0,0), 1e-32);
dimensionedScalar mupONE("zero", dimensionSet(1,-1,-1,0,0), 1.);

volScalarField mupLog
(
    IOobject
    (
        "mupLog", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
);

mupLog = log10(max(mup,mupLogsmall)/mupONE);
mupLog.write(); 

volScalarField mup2
(
    IOobject
    (
        "mup2", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mup2 = ( sigmapHAT && SHAT ) / ( max ( SHAT && SHAT, SSsmall ) );
// Limit by zero
//mup2.max(0.);
// Write into the results folder
mup2.write(); 

volScalarField mup3
(
    IOobject
    (
        "mup3", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mup3 = sqrt(( sigmap && sigmap ) / ( max ( S && S, SSsmall ) ));
// Limit by zero
//mup3.max(0.);
// Write into the results folder
mup3.write(); 

volScalarField mup4
(
    IOobject
    (
        "mup4", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mup4 = sqrt(( sigmapHAT && sigmapHAT ) / ( max ( SHAT && SHAT, SSsmall ) ));
// Limit by zero
//mup4.max(0.);
// Write into the results folder
mup4.write(); 

volScalarField mup5
(
    IOobject
    (
        "mup5", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);
dimensionedScalar SigmaSsmall("zero", dimensionSet(1,-1,-3,0,0,0,0), SMALL);
mup5 = ( sigmap && sigmap ) / ( ( sigmap && S ) + SigmaSsmall );
// Limit by zero
//mup5.max(0.);
// Write into the results folder
mup5.write(); 

volScalarField mup6
(
    IOobject
    (
        "mup6", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-1,0,0), scalar(0) )
);

mup6 = ( sigmapHAT && sigmapHAT ) / ( ( sigmapHAT && SHAT ) + SigmaSsmall );
// Limit by zero
//mup6.max(0.);
// Write into the results folder
mup6.write(); 

volScalarField shearstress
(
    IOobject
    (
        "shearstress", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

shearstress =  sqrt(sigmapHAT && sigmapHAT) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
shearstress.write(); 

volScalarField pressureuslip
(
    IOobject
    (
        "pressureuslip", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-3,0,0,0), scalar(0) )
);

pressureuslip =  Pp/(mag(U-mappedUs)*mag(U-mappedUs)) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
pressureuslip.write(); 

volScalarField sigmagradus
(
    IOobject
    (
        "sigmagradus", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-3,0,0), scalar(0) )
);

sigmagradus = -( tildeSigmaKin - sigmaColl )&&fvc::grad(mappedUs);
// Write into the results folder
sigmagradus.write(); 

dimensionedScalar sigmagradusLogsmall("zero", dimensionSet(1,-1,-3,0,0), 1e-32);
dimensionedScalar sigmagradusONE("zero", dimensionSet(1,-1,-3,0,0), 1.);
volScalarField sigmagradusLog
(
    IOobject
    (
        "sigmagradusLog", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
);

sigmagradusLog = log10(max(sigmagradus,sigmagradusLogsmall)/sigmagradusONE);
// Write into the results folder
sigmagradusLog.write(); 

volScalarField uslipuslip
(
    IOobject
    (
        "uslipuslip", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,1,-1,0,0), scalar(0) )
);

uslipuslip = mag(U-mappedUs);
// Limit by zero
//mup6.max(0.);
// Write into the results folder
uslipuslip.write(); 

volScalarField shearstress2
(
    IOobject
    (
        "shearstress2", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

shearstress2 =  sqrt(sigmap && sigmap) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
shearstress2.write(); 

// Particle shear stress
volTensorField sigmap_kin
(
    IOobject
    (
        "sigmap_kin", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimensionSet(1,-1,-2,0,0), tensor(0,0,0,0,0,0,0,0,0) )
);

//sigmap = ( tildeSigmaKin + sigmaColl ) - tensor(I) * Pp;
//sigmap = ( tildeSigmaKin - sigmaColl ) - tensor(I) * Pp;
sigmap_kin = tensor(I) * Ppkinetic - ( tildeSigmaKin  );
// Write into the results folder
sigmap_kin.write();

volScalarField shearstress2_kin
(
    IOobject
    (
        "shearstress2_kin", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

shearstress2_kin =  sqrt(sigmap_kin && sigmap_kin) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
shearstress2_kin.write(); 

volScalarField shearstress_xy
(
    IOobject
    (
        "shearstress_xy", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

shearstress_xy =  sqrt(sigmap.component(tensor::XY)*sigmap.component(tensor::XY)) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
shearstress_xy.write(); 

volScalarField shearstress_xzyz
(
    IOobject
    (
        "shearstress_xzyz", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-2,0,0), scalar(0) )
);

shearstress_xzyz =  sqrt((1./2.*(sigmap.component(tensor::XZ)+sigmap.component(tensor::YZ)))*(1./2.*(sigmap.component(tensor::XZ)+sigmap.component(tensor::YZ)))) ;
// Limit by zero
//mup6.max(0.);
// Write into the results folder
shearstress_xzyz.write(); 

volScalarField divU
(
    IOobject
    (
        "divU", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

divU =  fvc::div(mappedUs);
// Limit by zero
//mup6.max(0.);
// Write into the results folder
divU.write(); 

volScalarField ShearRate
(
    IOobject
    (
        "ShearRate", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

ShearRate =  sqrt( 2.0*(0.5*fvc::grad(mappedUs) + 0.5*fvc::grad(mappedUs)().T()-1.0/3.0*tr(fvc::grad(mappedUs))*tensor(1.,0.,0.,0.,1.,0.,0.,0.,1.)) &&
                       (0.5*fvc::grad(mappedUs) + 0.5*fvc::grad(mappedUs)().T()-1.0/3.0*tr(fvc::grad(mappedUs))*tensor(1.,0.,0.,0.,1.,0.,0.,0.,1.))     );
ShearRate.write(); 

volScalarField ShearRate2
(
    IOobject
    (
        "ShearRate2", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

ShearRate2 =  sqrt( 2.0 * (S/2.0) && (S/2.0));
// Limit by zero
//mup6.max(0.);
// Write into the results folder
ShearRate2.write();

volScalarField ShearRateMagGradU
(
    IOobject
    (
        "ShearRateMagGradU", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,-1,0,0), scalar(0) )
);

ShearRateMagGradU = Foam::sqrt(2.0)*mag(symm(fvc::grad(mappedUs)));
ShearRateMagGradU.write(); 

volScalarField Phi
(
    IOobject
    (
        "Phi", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(0,0,0,0,0), scalar(0) )
);

Phi = 1.-mappedVoidfraction;
Phi.write(); 

scalar partDia = radii[0][0];
volScalarField granularTemp = Ppkinetic/max((1-mappedVoidfraction),1e-12)/rhoParticle;
volScalarField g0 = (1.-Phi/2.)/Foam::pow(1.-Phi,3.) + 0.58*Phi*Phi/Foam::pow(Foam::sqrt(Foam::pow((0.61-Phi),2.)),1.5);

volScalarField GammaSlip
(
    IOobject
    (
        "GammaSlip", //+filterWidthName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "zero", dimensionSet(1,-1,-3,0,0), scalar(0) )
);

//- Some tricks about dimensions; dummy "1"s with dimemsion
dimensionedScalar granularTempInv("granularTempInv",granularTemp.dimensions(), 1.);
dimensionedScalar dimGammaSlip("dimGammaSlip",dimensionSet(1,-3,-1,0,0), 1.);

GammaSlip = dimGammaSlip*81.*Phi*1.8e-5*1.8e-5*uslipuslip*uslipuslip/(g0*rhoParticle*partDia*partDia*partDia*Foam::sqrt(pi*max(granularTemp,granularTempInv*1e-12)/granularTempInv));
GammaSlip.write(); 
