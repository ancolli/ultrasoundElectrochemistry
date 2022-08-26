#include "createRDeltaT.H"

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading field p_Re\n" << endl;
volScalarField p_Re
(
    IOobject
    (
        "p_Re",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_Im\n" << endl;
volScalarField p_Im
(
    IOobject
    (
        "p_Im",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing absolute acoustic pressure\n" << endl;
volScalarField ap
(
    IOobject
    (
        "ap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,//MUST_READ,//
        IOobject::AUTO_WRITE
    ),
    //mesh
    sqrt((pow(p_Re, 2) + pow(p_Im, 2))/2)
);

dimensionedScalar omega_f
(
    "omega_f",
    inv(dimTime),
    transportProperties.lookup("omega_f")
);

dimensionedScalar c
(
   "c",
   dimVelocity,
   transportProperties.lookup("c")
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

dimensionedScalar A
(
    "A",
    dimPressure*inv(dimArea),
    transportProperties.lookup("A")
);

dimensionedScalar B
(
    "B",
    dimPressure*inv(dimArea),
    transportProperties.lookup("B")
);

volScalarField k_r_sq
(
   IOobject
   (
      "k_r_sq",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   sqr(omega_f/c)  
);

volScalarField k_i_sq
(
   IOobject
   (
      "k_i_sq",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("k_i_sq", dimensionSet(0,-2,0,0,0,0,0), 0.0)
);

Info<< "Reading field p\n" << endl;
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

Info<< "Reading field U\n" << endl;
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

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimVelocity, vector(0,0,0))//fvc::grad(ap)/(omega_f)//rho*
);

surfaceScalarField phiUa
(
    IOobject
    (
        "phiUa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()
);

Info<< "Reading field C\n" << endl;
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField Flux
    (
        IOobject
        (
            "Flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE//AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", C.dimensions()/dimTime*dimLength, 0.0)
    );

Info<< "Reading field Deff\n" << endl;

    volScalarField Deff
    (
        IOobject
        (
            "Deff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,//MUST_READ,/
            IOobject::AUTO_WRITE//NO_WRITE
        ),
        mesh,
        dimensionedScalar("Deff",dimensionSet(0,2,-1,0,0,0,0),0.0)
    );

 Info<< "Reading diffusivity D\n" << endl;

    

    dimensionedScalar Dmol
    (
        "Dmol",
   	dimVelocity*dimLength,
   	transportProperties.lookup("Dmol")
    );
    
    Info<< "Reading turbulent Schmidt number\n" << endl;

    dimensionedScalar Sct
    (
    	"Sct",
   	dimless,
        transportProperties.lookup("Sct")
    );

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
