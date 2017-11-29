// Post-processing dictionary 
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingPropertiesADMofScalarVariance",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Define couplingProperties dictionary
const dictionary dict(couplingProperties);

// Define postprocessing dictionary
const dictionary postProcessingDict(dict.subDict("postProcessingProps"));					

// Fluid coarsening
dictionary fluidCoarseningDict(postProcessingDict.subDict("fluidCoarsening"));
labelList filterWidth(labelList(fluidCoarseningDict.lookup("filterWidth")));

// ADM dictionary
dictionary ADMDict(postProcessingDict.subDict("ADMDict"));
label ADMOrder(readLabel(ADMDict.lookup("order")));
Info << "ADM will be performed with " << ADMOrder << "th order" << endl;

