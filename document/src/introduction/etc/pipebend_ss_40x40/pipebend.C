#include "DynLinkTools/StandardIncludes.h"

/// Element Information Creation Block
static const UInt sBSetIds[] = {1};

static const enum XFemElemType  sXfemBlockElemType  = NON_XFEM_ELEMENT;
static const UInt sNumDofsPerNodeMax                = 3;

/// Time Integration options

static const Real sEpsilon               = 1.0e-5 ;
static const Real smthng_rad             = 0.075;
static const UInt sNumTimeIts            = 3;

//----------------------------------------------------------------------------------------------------
//
//           THIS IS WHERE MOST OF THE CHANGES WILL OCCUR
//
//----------------------------------------------------------------------------------------------------
extern "C"
    void DLModel ( const char* aModelFile, ModelData* aModelDataStruct, Problem *problem )
{
    // Local variables
    ParseModelHelper parseModelHelper = ParseModelHelper();

    int PhaseIds[]   = {-1, 1, 0, 0, 0, 0};


    UInt modelId        = aModelDataStruct->GetModelId();
    UInt meshId         = aModelDataStruct->GetMeshId();
    UInt numSolnWriters = aModelDataStruct->mSolutionData.size();

/// MODEL --------------------------------------------------------
    Model* model = problem->CreateModel ( aModelDataStruct,modelId,meshId,0,0,0,0,
                                          0,0,0,numSolnWriters);
    Mesh* msh = model->GetMesh();

    parseModelHelper.SetNumLs(model,1);
    parseModelHelper.SetMainPhaseFlags(model,PhaseIds,2);

/// Functions
    // Physical design variables function creation
    UInt numE[3] = {0,0,0};
    UInt *IndicesList=NULL, *IdsList=NULL;
    UInt numPropFunctions = parseModelHelper.GetBaseSetElemsInfoByID(msh,sBSetIds[0],IndicesList,IdsList);
    std::vector<UInt>       ElemIndicesVec(IndicesList,IndicesList+numPropFunctions);
    std::vector<PropFunction*>  EmodFunctionPtrsVec(numPropFunctions,NULL);
    std::vector<PropFunction*>  RhoFunctionPtrsVec(numPropFunctions,NULL);
    Real funcVal            = 1.0;
    enum FunctionVarType funcVarType = FUNC_VAR_DESVAR;
    for ( UInt in = 0; in < numPropFunctions; ++in )
    {
        EmodFunctionPtrsVec[in] = model->CreatePropFunction(                 in+1,FUNCTION_CONSTANT,&funcVarType, numE,&funcVal);
	RhoFunctionPtrsVec[in]  = model->CreatePropFunction(numPropFunctions+in+1,FUNCTION_CONSTANT,&funcVarType, numE,&funcVal);
    }


	std::vector<UInt> NpIds;
	NpIds.push_back(10);

	model->CreateNodeProperty ( NpIds[0],NP_STRUC );


/// Element Properties block
	std::vector<UInt>    BlkEpIds;
	BlkEpIds.push_back(20);
    model->CreateElementProperty ( BlkEpIds[0],EP_FLUID );
	std::vector<UInt>    BdyBlkEpIds;
	BdyBlkEpIds.push_back(21);
    ElementProperty*  elemProp = model->CreateElementProperty( BdyBlkEpIds[0],EP_FLUID);

    elemProp->SetDirective( (UInt) ED_FLUID_BNDY_INTEG, BI_PRES );

/// Materials Block
	std::vector<UInt> BlkMatIds;
	BlkMatIds.push_back(30);
    Material* material = model->CreateMaterial ( BlkMatIds[0], MATERIAL_FLUID_IDEAL);
    material->InitializeProperty( MP_FLUID_ICMP_RYNLDS , .1, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_DNSTY , 1.0, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_GSV,      1.0e-6, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_ALPHA,    0.2, numPropFunctions, &ElemIndicesVec[0], numPropFunctions, &EmodFunctionPtrsVec[0] );
    material->InitializeProperty( MP_FLUID_ICMP_ATHRESH, 0.0001, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_SWITCHBW,0.0001, 1, NULL, 1, NULL );
    material->InitializeProperty(MP_FLUID_ICMP_VISC, 1./.1, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_SOLID_DENSITY,0.2, numPropFunctions, &ElemIndicesVec[0], numPropFunctions, &RhoFunctionPtrsVec[0] );
//    UInt  SideMatIds[] = {31,31,31,31,31,31};
    std::vector<UInt> SideMatIds;
    SideMatIds.push_back(31);    
    material = model->CreateMaterial ( 31, MATERIAL_FLUID_IDEAL);
    material->InitializeProperty( MP_FLUID_ICMP_RYNLDS , .1, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_DNSTY , 1.0, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_GSV,      1.0e-6, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_ALPHA,    0.2, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_ATHRESH, 0.0001, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_SWITCHBW,0.0001, 1, NULL, 1, NULL );
    material->InitializeProperty(MP_FLUID_ICMP_VISC, 1./.1, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_FLUID_ICMP_SOLID_DENSITY,0.2, 1, NULL, 1, NULL );

    /// Nodes
    if (sXfemBlockElemType == NON_XFEM_ELEMENT)
        parseModelHelper.CreateAllNodes(model,NpIds);
    else
        parseModelHelper.CreateAllNodes(model,NpIds);

/// Block Elements
    enum ElementType BlockElemType = FLUID_ICMP_QUAD4;
    parseModelHelper.CreateElementsByBaseSetId(model,sBSetIds[0],NON_COUPLED_ELEMENT,NON_IC_ELEMENT,sXfemBlockElemType,&BlockElemType,BlkEpIds,BlkMatIds);

//    UInt SideEpIds[] = {21,21,21,21,21,21};
std::vector<UInt> SideEpIds;
SideEpIds.push_back(21);
    UInt threeB = 3;

    enum ElementType SideElemType = FLUID_ICMP_SURF;
    parseModelHelper.CreateElementsBySideSetId(model,threeB,NON_COUPLED_ELEMENT,NON_IC_ELEMENT,NON_XFEM_ELEMENT,&SideElemType,SideEpIds,SideMatIds);

/// Finalize Model
    parseModelHelper.FinalizeModel(model);

///MPC Block
    model->SetMpcDofs();

    UInt InflowSSId = 1, ChannelTopBottomSSId = 3, CylSSId = 2;

/// Dirichlet Conditions
    parseModelHelper.ApplyNodeConditionToSSetById( model, ChannelTopBottomSSId, DIRICHLET, REGULAR_DOF, RHO, 1.0, 0);

    parseModelHelper.ApplyNodeConditionToSSetById( model, CylSSId, DIRICHLET, REGULAR_DOF, RHOVX, 0.0, 0);
    parseModelHelper.ApplyNodeConditionToSSetById( model, CylSSId, DIRICHLET, REGULAR_DOF, RHOVY, 0.0, 0);

    parseModelHelper.ApplyNodeConditionToSSetById( model, InflowSSId, DIRICHLET, REGULAR_DOF, RHOVX, 1.0, 0);
//    parseModelHelper.ApplyNodeConditionToSSetById( model, InflowSSId, DIRICHLET, REGULAR_DOF, RHOVX, 1.0, 0); // don't use the BC Function above
    parseModelHelper.ApplyNodeConditionToSSetById( model, InflowSSId, DIRICHLET, REGULAR_DOF, RHOVY, 0.0, 0);

/// Initial Conditions Block
    parseModelHelper.ApplyNodeConditionToAll( model, IDISP, REGULAR_DOF, RHOVX, 0.0001, 0);
    parseModelHelper.ApplyNodeConditionToAll( model, IDISP, REGULAR_DOF, RHOVY, 0.0, 0);
    parseModelHelper.ApplyNodeConditionToAll( model, IDISP, REGULAR_DOF, RHO, 1.0, 0);

    parseModelHelper.ApplyNodeConditionToSSetById( model, InflowSSId, IDISP, REGULAR_DOF, RHOVX, 0., 0);
    parseModelHelper.ApplyNodeConditionToSSetById( model, InflowSSId, IDISP, REGULAR_DOF, RHOVY, 0.0, 0);
}


extern "C"
    ProblemData* InputData ( Problem *problem )
{
#ifdef PARALLEL
    // Get the run suffix to enable viewing with paraview
    int myPID = 0, numProcs=1, strSize=30;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPID), MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int numDigits = std::floor(std::log10(numProcs))+1;
    char *exoSuffix=(char*) alloca(sizeof(char)*strSize), formatStr[] = ".%i.%0\0\0\0\0\0\0\0\0", *nemFolder=(char*) alloca(sizeof(char)*strSize);
    std::strncpy(exoSuffix,"\0",strSize), std::strncpy(nemFolder,"\0",strSize);
    std::sprintf(formatStr+6,"%ii",numDigits);
    std::sprintf(exoSuffix,formatStr,numProcs,myPID);
    std::sprintf(nemFolder,"./nem-%i/",numProcs);
#endif

    // Local variables
    ParseInputHelper parseInputHelper = ParseInputHelper();

    // Create the ProblemData,MeshData,SolverData, and ModelData
    ProblemData* pd        = problem->GetProblemData();
    MeshData*   meshData   = pd->CreateMeshDataStruct();
    SolverData* solverData = pd->CreateSolverDataStruct();
    ModelData*  modelData  = pd->CreateModelDataStruct();
    meshData->mMeshFileData           = new MeshFileData();
    solverData->mTimeSolverData       = new TimeSolverData();
    solverData->mNonlinearSolverDataList.push_back(new NonlinearSolverData());
    solverData->mLinearSolverDataList.push_back(new LinearSolverData());
    modelData->mModelFileData         = new ModelFileData();
    modelData->mFluidModelData        = new FluidModelData();
    SolutionData* solutionData        = new SolutionData();
    modelData->mSolutionData.push_back( solutionData );


///Set Problem and Mesh
    meshData->SetMeshId ( 1 );
    meshData->mMeshFileData->mMeshFileType = MESH_EXODUS;
    const char* InputMeshName             = "pipebend_40x40.g";
#ifdef PARALLEL
    std::sprintf(meshData->mMeshFileData->mMeshFileBase, "%s%s"  , nemFolder, InputMeshName );
    std::sprintf(meshData->mMeshFileData->mMeshFile    , "%s%s%s", nemFolder, InputMeshName, exoSuffix );
#else
    std::sprintf(meshData->mMeshFileData->mMeshFileBase, "%s", InputMeshName );
    std::sprintf(meshData->mMeshFileData->mMeshFile    , "%s", InputMeshName );
#endif


/// Time Integration options
    solverData->SetSolverId ( 1 );
    solverData->mTimeSolverData->mType           = TIME_SOLVER_FEMDOC_BDF;
    solverData->mTimeSolverData->mBdfOrder       = 1;
    solverData->mTimeSolverData->mTime           = 0.0;
    solverData->mTimeSolverData->mTimeIts        = sNumTimeIts;
    solverData->mTimeSolverData->mTimeItList[0]  = solverData->mTimeSolverData->mTimeIts;
    solverData->mTimeSolverData->mDt             = 10000.;
    solverData->mTimeSolverData->mDtList[0]      = solverData->mTimeSolverData->mDt;
//    solverData->mTimeSolverData->mMaxTime         = solverData->mTimeSolverData->mTimeIts*solverData->mTimeSolverData->mDt;
    solverData->mTimeSolverData->mBeta           = 0.25;           // Newmark Alpha
    solverData->mTimeSolverData->mGamma          = 0.5;            // Newmark Delta
    solverData->mTimeSolverData->mAlphaF         = 0.0;
    solverData->mTimeSolverData->mAlphaM         = 0.0;
    solverData->mTimeSolverData->mTotResNorm     = 1.0e-12;
    solverData->mTimeSolverData->mTotResNormDrop = 1.0e-12;
    solverData->mTimeSolverData->mNondimDivisor  = 1.0;

    solverData->mTimeSolverData->mJacCheck = false;
    solverData->mTimeSolverData->mJacFree = false;
    solverData->mTimeSolverData->mNaNCheck = true;
    
            solverData->mTimeSolverData->mReinitByStandaloneFile = true;
    const char* tempPath = "tempDir/";
    solverData->mTimeSolverData->mSolVecDir = new char[200];
    strcpy(solverData->mTimeSolverData->mSolVecDir,tempPath);


/// Newton solver and linear solver options
    solverData->mNonlinearSolverDataList[0]->mType             = NONLINEAR_SOLVER_FEMDOC_NEWTON;
    solverData->mNonlinearSolverDataList[0]->mMaxItsList[0]    = 15;
    solverData->mNonlinearSolverDataList[0]->mRelaxation       = 1.0;
    solverData->mNonlinearSolverDataList[0]->mRlxList[0]       = solverData->mNonlinearSolverDataList[0]->mRelaxation;
    solverData->mNonlinearSolverDataList[0]->mTotResNormDrop   = 1.0e-6;
    solverData->mNonlinearSolverDataList[0]->mTotResNorm       = 1.0e-20;


/// Linear solver options
    solverData->mLinearSolverDataList[0]->mSparseMatrixType    = SPARSE_MATRIX_TRILINOS_EPETRA_FECRS;
#ifdef PARALLEL
    solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AZTEC;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_CG;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_CGS;
    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_BICGSTAB;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_GMRES;
    solverData->mLinearSolverDataList[0]->mPreconditioner      = LINSOL_PREC_AZTEC_DD_ILUT;
    solverData->mLinearSolverDataList[0]->mKrylovSpace         = 30;
    solverData->mLinearSolverDataList[0]->mMaxIts              = 1000;
    solverData->mLinearSolverDataList[0]->mDropTol             = 1e-20;
    solverData->mLinearSolverDataList[0]->mGraphFill           = 0;
    solverData->mLinearSolverDataList[0]->mIlutFill            = 10.0;
    solverData->mLinearSolverDataList[0]->mPolyOrder           = 4;
    solverData->mLinearSolverDataList[0]->mAlphaThresh         = 0.00000;
    solverData->mLinearSolverDataList[0]->mRhoThresh           = 0.0;
    solverData->mLinearSolverDataList[0]->mConvergence         = LINSOL_CONV_R0;
    solverData->mLinearSolverDataList[0]->mDiagnostics         = LINSOL_DIAGNOSTICS_NONE;
    solverData->mLinearSolverDataList[0]->mEpsilon             = 1.0e-10;
    solverData->mLinearSolverDataList[0]->mOverlap             = 1;
    solverData->mLinearSolverDataList[0]->mReorder             = 1;

//    solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AMESOS;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_MUMPS;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_UMFPACK;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_KLU;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_SUPERLUDIST;
#else
//    solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AZTEC;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_CG;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_CGS;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_BICGSTAB;
////    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AZTEC_GMRES;
//    solverData->mLinearSolverDataList[0]->mPreconditioner      = LINSOL_PREC_AZTEC_DD_ILUT;
//    solverData->mLinearSolverDataList[0]->mKrylovSpace         = 30;
//    solverData->mLinearSolverDataList[0]->mMaxIts              = 30;
//    solverData->mLinearSolverDataList[0]->mDropTol             = 1e-20;
//    solverData->mLinearSolverDataList[0]->mGraphFill           = 0;
//    solverData->mLinearSolverDataList[0]->mIlutFill            = 6.0;
//    solverData->mLinearSolverDataList[0]->mPolyOrder           = 3;
//    solverData->mLinearSolverDataList[0]->mAlphaThresh         = 0.00000;
//    solverData->mLinearSolverDataList[0]->mRhoThresh           = 0.0;
//    solverData->mLinearSolverDataList[0]->mConvergence         = LINSOL_CONV_R0;
//    solverData->mLinearSolverDataList[0]->mDiagnostics         = LINSOL_DIAGNOSTICS_NONE;
//    solverData->mLinearSolverDataList[0]->mEpsilon             = 1.0e-16;
//    solverData->mLinearSolverDataList[0]->mOverlap             = 1;
//    solverData->mLinearSolverDataList[0]->mReorder             = 1;

    solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AMESOS;
    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_UMFPACK;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_KLU;
//    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_SUPERLU;
#endif
//    solverData->mLinearSolverDataList[0]->mExport              = LINSOL_EXPORT_MATLAB;


/// Model Data
    modelData->SetModelId( 1 );
    modelData->mMeshId                        = 1;
    modelData->mSolverId                      = 1;
    modelData->mModelFileData->mModelFileType = MODEL_FEMDOC;
//    modelData->mStrucModelData->mRedDimType   = FULL3D;



    modelData->mFluidModelData->mJacCoupling.FullyCoupled = false;
    modelData->mFluidModelData->mEquation.EqnType = FLUID_EQN_INS;
    modelData->mFluidModelData->mEquation.TauType = IFLUID_STAB_TAUTT;
    modelData->mFluidModelData->mEquation.IncludeAleTerms = false;
    modelData->mFluidModelData->mEquation.Nondimensional = false;
    modelData->mFluidModelData->mFreestream.Length = 1.0;
    modelData->mFluidModelData->mFreestream.Velocity = 1.0;
    modelData->mFluidModelData->mReference.Area = 0.1;
    modelData->mFluidModelData->mFreestream.Viscosity = 1./.1;
    modelData->mFluidModelData->mFreestream.Density = 1.;
    modelData->mFluidModelData->mFreestream.Velocity = 1.;

/// Solution writing format
    solutionData->mSolnFormatType = SF_EXODUS;
    const char* OutputMeshBase   = "pipebend.e-s";
    solutionData->mSaveIC         = 0;
    solutionData->mOptSaveFreq    = 1;
    solutionData->mSaveOpt        = true;
    solutionData->mIsXFem         = (sXfemBlockElemType!=NON_XFEM_ELEMENT);
    strcpy(solutionData->mSolnFileBase,OutputMeshBase);
    strcpy(solutionData->mSolnFileName,OutputMeshBase);
#ifdef PARALLEL
    std::sprintf(solutionData->mSolnFileBase, "%s%s", nemFolder, OutputMeshBase );
    std::sprintf(solutionData->mSolnFileName, "%s%s", nemFolder, OutputMeshBase );
#else
    std::sprintf(solutionData->mSolnFileBase, "%s", OutputMeshBase );
    std::sprintf(solutionData->mSolnFileName, "%s", OutputMeshBase );
#endif
/// Post Processing Output Variables
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_RHOVX               , ELEM_VAR, 1,sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_RHOVY               , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_RHO               , ELEM_VAR, 1, sBSetIds );

    UInt one = 1;
    UInt three = 3;
    
    
    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_MP_FLUID_ALPHA     ,           1, sBSetIds );
    
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_MASS_UTILIZATION, ELEM_VAR,    1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_FLUID_TP        , SIDESET_VAR, 1, &one     , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_FLUID_TP        , SIDESET_VAR, 1, &three   , 0, NULL);

    return pd;
}



// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//
// Optimizer dynamically linked functions
//
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
extern "C"
    void ComputeObjectivesAndConstraints(const std::vector <Real> & aCriteriaVals, Real & aObjectiveVal, std::vector <Real> & aConstraintVals)
{
    /// ----------------------------------------------------------------------------
    /// Objective is strain energy with a perimeter penalty
    /// ----------------------------------------------------------------------------
    aObjectiveVal = (aCriteriaVals[0]-aCriteriaVals[1]);

    /// ----------------------------------------------------------------------------
    /// Constraint is a limiting volume of phase 1
    /// ----------------------------------------------------------------------------
    aConstraintVals[0] = aCriteriaVals[2]/.25-1.;
}

extern "C"
    void Compute_dOptScalarsDQs(const std::vector< Real > & aCriteriaVals, std::vector< Real > & aGrad)
{
    UInt numCriteria   = 3;

    aGrad[0*numCriteria + 0] = 1.;
    aGrad[0*numCriteria + 1] = -1.;
    aGrad[0*numCriteria + 2] = 0.;

    aGrad[1*numCriteria + 0] = 0.;
    aGrad[1*numCriteria + 1] = 0;
    aGrad[1*numCriteria + 2] = 1./.25;
}



///----------------------------------------------------------------------------------------
/// function for defining the paramters of the optimization
///----------------------------------------------------------------------------------------
extern "C"
    void SetOptimizationParameters(Optimizer *aOptimizer,
                                   Problem   *aProblem)
{
/*
    // Finite difference checker and dR_ds perturbation value
    aOptimizer->SetAdvEpsilon( 1.0e-5 );

    // Set whether or not to use "wiggling" (true = don't wiggle)
    aOptimizer->SetMapFlag(false);

    // Set the number of constraints
    aOptimizer->SetNumConstraints(1);

    // Set the optimizer and algorithmic parameters
    OptimizerData* optSolData = aOptimizer->GetOptimizerData();
    optSolData->mType     = OPT_SOLVER_SWEEP;
    optSolData->mMaxIts   = 75;
    optSolData->mNormDrop = 1e-6;

    optSolData->mGcmmaData.AsympAdapt0 = 0.5;
    optSolData->mGcmmaData.AsympAdapt  = 0.7;
    optSolData->mGcmmaData.AsympAdaptC = 1./0.7;
    optSolData->mGcmmaData.StepSize    = 0.05;
    optSolData->mGcmmaData.Penalty     = 1.0;



 static const UInt mNumVar			= 100;					// Number of variables to be analyzed
// static UInt mNumSteps[mNumVar]		= {    1,    1,    1};	// Number of steps between variable lower limit and initial variable
// static UInt mVarId[mNumVar]			= {   1,   2,   3};	// Variable ID's
// static Real mVarLow[mNumVar]		= { 0.3, 0.3, 0.3};	// Lower LsVals for each variable ID respectively
// static Real mVarUp[mNumVar]			= {  0.3,  0.3,  0.3};	// Upper LsVals for each variable ID respectively
// static bool mFDCheck[mNumVar]		= {    1,    1,    1};	// Finite difference check switch
 static char mSweepFileName[]		= "sweepResults.txt";
// static Real eps[mNumVar] = {1.0e-5,1.0e-5,1.0e-5};
// static bool fuck[mNumVar] = {true,true,true};
     optSolData->mSweepData.mNumVar			= mNumVar;
//     optSolData->mSweepData.mNumSteps		= mNumSteps;
//     optSolData->mSweepData.mVarId			= mVarId;
//     optSolData->mSweepData.mVarLow			= mVarLow;
//     optSolData->mSweepData.mVarUp			= mVarUp;
//     optSolData->mSweepData.mFDCheck		= mFDCheck;
//     optSolData->mSweepData.mFDEpsilon		= eps;
     optSolData->mSweepData.mSweepFileName	= mSweepFileName;
//     optSolData->mSweepData.mGradFlag            = fuck;


//UInt Ids[] = {72,200,50};
    std::vector<UInt> NumSteps;
    std::vector<UInt> VarId;
    std::vector<Real> VarLow;
    std::vector<Real> VarUp;
    std::vector<UInt> GradFlag;
    std::vector<UInt> FDCheck;
    std::vector<Real> FDEpsilon;
    for (UInt id = 0; id < mNumVar; ++id)
    {
        NumSteps.push_back(0);
//        VarId.push_back(Ids[id]);
        VarId.push_back(id);
        VarLow.push_back(0.3);
        VarUp.push_back(0.3);
        GradFlag.push_back(1);
        FDCheck.push_back(1);
        FDEpsilon.push_back(1.0e-5 );

    }

    optSolData->mSweepData.SetNumStepsVector(NumSteps);
    optSolData->mSweepData.SetVarIdVector(VarId);
    optSolData->mSweepData.SetVarLowVector(VarLow);
    optSolData->mSweepData.SetVarUpVector(VarUp);
    optSolData->mSweepData.SetGradFlagVector(GradFlag);
    optSolData->mSweepData.SetFDCheckVector(FDCheck);
    optSolData->mSweepData.SetFDEpsilonVector(FDEpsilon);

    optSolData->mDumpAdvs = ADV_LAST;
    optSolData->mRestartIteration = 0;
    aOptimizer->mWriteOptScalarValuesToFile = true;

*/


    // Finite difference checker and dR_ds perturbation value
    aOptimizer->SetAdvEpsilon( 1.0e-5 );

    // Set whether or not to use "wiggling" (true = don't wiggle)
    aOptimizer->SetMapFlag(true);

    // Set the number of constraints
    aOptimizer->SetNumConstraints(1);

    // Set the optimizer and algorithmic parameters
    OptimizerData* optSolData = aOptimizer->GetOptimizerData();
    optSolData->mType     = OPT_SOLVER_OLD_FEM_GCMMA;
    optSolData->mMaxIts   = 100;
    optSolData->mNormDrop = 1e-6;

    optSolData->mGcmmaData.AsympAdapt0 = 0.5;
    optSolData->mGcmmaData.AsympAdapt  = 0.7;
    optSolData->mGcmmaData.AsympAdaptC = 1./0.7;
    optSolData->mGcmmaData.StepSize    = 0.05;
    optSolData->mGcmmaData.Penalty     = 1000.;

    optSolData->mDumpAdvs = ADV_LAST;
    optSolData->mRestartIteration = 0;
    aOptimizer->mWriteOptScalarValuesToFile = true;


}

///----------------------------------------------------------------------------------------
/// function for defining the optimization criteria
///----------------------------------------------------------------------------------------
extern "C"
    void SetOptimizationCriteria(Optimizer *aOptimizer,
                                 Problem   *aProblem)
{
    OptimizationCriteria* crit = NULL;
    Function* CritEvalFunc = NULL;
    Real weights = 1.;
    UInt entries[] = {1,0};
    UInt indices = sNumTimeIts; // Last time step

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_FLUID_TP,true,false);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_FLUID_TP,true,false);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_MASS_UTILIZATION,false,true);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

}


///----------------------------------------------------------------------------------------
/// function for defining the Physical design variables and Abstract design variables
///----------------------------------------------------------------------------------------
extern "C"
    void SetOptimizationPdvsAndAdvs(Optimizer *aOptimizer,
                                    Problem   *aProblem)
{
   /// Set up design variables
   // FIXME: right now only set up for one (first) mesh
    Mesh* mesh = aProblem->GetMesh(0);
    UInt modelId = 1;

    // for setting the design variables
    UInt  blockId = 0, cellSetIndex = 0, myNumCellsInSet = 0;
    CellSet*      cset = NULL;
    UInt* cellIndexList = NULL;
        Model*     model = aProblem->GetModelById(modelId);
	model->BuildNodeToElementTable();
	
    // Vars for creating PDVs
    DesignDepObject desDepObj = MAT_PROP;
    UInt propId      = 0;
    UInt propTypeKey = 0;

    // Get block 1 info
    blockId = 1;
    cellSetIndex    = mesh->GetCellSetIndexFromBaseCellSetId(blockId);
    cset            = mesh->GetCellSet(cellSetIndex);
    cellIndexList   = cset->GetCellSetCellIndices();
    myNumCellsInSet = cset->GetNumCellsInSet();


    // Create physical variables for the element's DENSITY in block 1
    desDepObj   = MAT_PROP;
    propId      = 30;
    int propTypeAlpha    = MP_FLUID_ICMP_ALPHA;
    int propTypeSolidDens    = MP_FLUID_ICMP_SOLID_DENSITY;
    
    // Set the number of design variables per processor
    aOptimizer->SetNumDesVars(myNumCellsInSet);
    
 
    PhysicalDesignVar* pdv = NULL;
    ///////////////////////////////////////////////////////////////////////////
    // Assign Pdv's 
    ///////////////////////////////////////////////////////////////////////////
    for (UInt ie = 0; ie < myNumCellsInSet; ++ie)
    {
        propTypeKey  = cellIndexList[ie];

            UInt paramIndex = 0;
        pdv = aOptimizer->CreatePhysDesVar(aProblem,modelId,desDepObj,propId,propTypeAlpha,propTypeKey,paramIndex);
/// Create the evaluation function for this physical design variable
            Real vals[] = {1.,0.,25000.,0.01};
            UInt numDstncs = 1;
            Function* PDVEvalFunc = new FunctionCnvxInterpolation(vals);
            pdv->CreateEvaluationFunction(PDVEvalFunc,&propTypeKey);
	    
	    
	    pdv = aOptimizer->CreatePhysDesVar(aProblem,modelId,desDepObj,propId,propTypeSolidDens,propTypeKey,paramIndex);
/// Create the evaluation function for this physical design variable
Real one = 1.;
            PDVEvalFunc = new FunctionScalarProduct(&numDstncs,&one);
            pdv->CreateEvaluationFunction(PDVEvalFunc,&propTypeKey);
    }
}



///----------------------------------------------------------------------------------------
/// Initialize the Vector of abstract design variables (Adv) mAbsDesVarVec
///----------------------------------------------------------------------------------------
extern "C"
    void InitializeAbsDesVarVector(Optimizer* aOptimizer, UInt aNumMyAbsDesVars)
{
    Problem* problem = aOptimizer->GetProblem();
    Mesh* mesh = problem->GetMesh(0);
    Vector* absDesVarVec = aOptimizer->GetAbsDesVarVec();
    Element* ele = NULL;
    Real* glbCoords = (Real*) alloca(sizeof(Real)*2);

    // for setting the design variables
    UInt  blockId = 0, cellSetIndex = 0, myNumCellsInSet = 0;
    CellSet*      cset = NULL;
    UInt* cellIndexList = NULL;



    blockId = 1;
    cellSetIndex    = mesh->GetCellSetIndexFromBaseCellSetId(blockId);
    cset            = mesh->GetCellSet(cellSetIndex);
    cellIndexList   = cset->GetCellSetCellIndices();
    myNumCellsInSet = cset->GetNumCellsInSet();

    // FIXME right now only setup for 1 model
    UInt modelId = 1;
    Model* model = problem->GetModelById(modelId);

    Real sVal = 0.0e0;


    
    Real* Ex = (Real*) alloca(sizeof(Real)*4);
    Real* Ey = (Real*) alloca(sizeof(Real)*4);
    Real* Ez = NULL;
    
    UInt numFixedElems = 0;
    UInt* fixedElems = (UInt*) alloca(sizeof(UInt)*myNumCellsInSet);
    Real* lowerLimForFixedElems = (Real*) alloca(sizeof(Real)*myNumCellsInSet);
    Real* upperLimForFixedElems = (Real*) alloca(sizeof(Real)*myNumCellsInSet);
   
    UInt sIndex = 0;

    for (UInt ie = 0; ie < myNumCellsInSet; ++ie)
    {
        ele = model->GetElement(cellIndexList[ie]);
	
	ele->GetElementNodeCoords(Ex,Ey,Ez);
	
        glbCoords[0] = Ex[0] + 0.5*(Ex[1] - Ex[0]);
        glbCoords[1] = Ey[0] + 0.5*(Ey[3] - Ey[0]);
	
	sIndex = ie;
	if ( ( glbCoords[0] < 0.125e0 ) &
	   ( ( glbCoords[1] < 2.0e0 ) & ( glbCoords[1] > 1.0e0) ) ) 
	{
	  fixedElems[numFixedElems] = sIndex;
	  lowerLimForFixedElems[numFixedElems] =  0.99e0;
	  upperLimForFixedElems[numFixedElems] =  0.999e0;
	  ++numFixedElems;
	  sVal = 0.99e0;
	}
	else if ( ( glbCoords[1] < -2.375e0 ) &
	        ( ( glbCoords[0] < 4.5e0 ) & ( glbCoords[0] > 3.5e0) ) ) 
	{
	  fixedElems[numFixedElems] = sIndex;
	  lowerLimForFixedElems[numFixedElems] =  0.99e0;
	  upperLimForFixedElems[numFixedElems] =  0.999e0;
	  ++numFixedElems;
	  sVal = 0.99e0;
	}
        else 
	{
          sVal = 0.9e0;
        }
        absDesVarVec->ReplaceGlobalValues(1,&sIndex,&sVal);
    }

    // set the bounds for mAbsDesVarVec 
    const Real AdvLowerBound = 1.0e-2;
    const Real AdvUpperBound =  1.0-1.0e-2;

    aOptimizer->FillAbsDesVarVecBounds(AdvLowerBound,AdvUpperBound);

    aOptimizer->SetUpperAndLowerBounds(numFixedElems,
                                       fixedElems, 
                                       lowerLimForFixedElems,
                                       numFixedElems,
                                       fixedElems, 
                                       upperLimForFixedElems);


}




