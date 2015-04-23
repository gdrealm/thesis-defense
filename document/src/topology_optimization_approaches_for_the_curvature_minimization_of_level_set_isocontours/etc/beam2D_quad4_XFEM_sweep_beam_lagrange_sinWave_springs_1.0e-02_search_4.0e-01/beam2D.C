#include <fstream>
#include "DynLinkTools/StandardIncludes.h"

/// Element Information Creation Block
static const UInt sBSetIds[] = {1};
static const UInt sSSetIds[] = {1,2,3,4,5};
static const UInt sNSetIds[] = {1,2,3,4,5,6};
static const UInt sAnchNSSetIndx = 3;  // Inverter anchor nodeset&sideset index
static const UInt sLoadNSSetIndx = 4;  // Load surface nodeset&sideset index

static const Real sCenterX = 1.500000;
static const Real sCenterY = 1.000000;
static const Real sCenterZ = 0.000000;

static const Real sRadiusSteps = 100;
static const Real sRadiusLower = 0.25;
static const Real sRadiusUpper = 0.75;

static const enum XFemElemType  sXfemBlockElemType  = STRUC_XFEM_ELEMENT;
static const UInt sNumDofsPerNodeMax                = 2;

/// Time Integration options
static const Real sTotalTime             = 200000.;
static const Real sTimeStep              = 200000.;
static const UInt sNumIts                = sTotalTime/sTimeStep;

/// -------------------------------------------------------------------------------------
///  SECTION FOR OPTIMIZATION INPUTS
/// -------------------------------------------------------------------------------------
static const Real sSmthng_rad = 0.1;

static const Real sPerPen    = 0.0e+00;
static const Real sCurPen    = 1.0e+03;
static const Real sDispPen   = 1.0e+02;
static const Real sMassRatio = 0.50;
static const UInt sNumOptIts = 0;

static const Real sCurvatureSearchRadius = 0.4;
static const enum GlbSolVarType GSOL_CURVATURE_TYPE = GSOL_INTERFACE_CURVATURE_BEAM_LAGRANGE;

static const bool sSaveIC = 0;

static const bool sOptFDCheck = true;
static const Real sOptFDEpsilon = 1.0e-03;
static const Real sOptAdvEpsilon = 1.0e-04;

static const UInt sDesignDomainBlock = 1;

// User defined function declaration
Real mySweepFunc(const Real* aIndep, const Real* aParameters);
void mySweepFuncDeriv(Real* aDerivVals, const Real* aIndep, const Real* aParameters, const enum DerivativeOrder& aDerivOrder, const Real* aDCoefs);

// ----------------------------------------------------------------------------

Real mySweepFunc(const Real* aIndep, const Real* aParameters)
{
	Real x = aParameters[0];
	Real y = aParameters[1];

	x = x; y = y;

	return y - aIndep[0] * std::sin( MATH_PI * x );
}

void mySweepFuncDeriv(Real* aDerivVals, const Real* aIndep, const Real* aParameters, const enum DerivativeOrder& aDerivOrder, const Real* aDCoefs)
{
	return;
}

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

    UInt numFunctions   = 0;

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
    UInt numE[3] = {0,0,0};
    UInt *NodeIndices=NULL, *NodeIds=NULL;
    numFunctions = parseModelHelper.GetAllNodesInfo(msh,NodeIndices,NodeIds);
    PropFunction** PropFunctionPtrs = (PropFunction**) alloca(sizeof(PropFunction*)*numFunctions);
    Real funcVal = 1.0;
    const enum FunctionVarType funcVarType = FUNC_VAR_DESVAR;
    for ( UInt in = 0; in < numFunctions; ++in )
        PropFunctionPtrs[in] = model->CreatePropFunction(NodeIds[in],FUNCTION_CONSTANT,&funcVarType, numE,&funcVal);

/// Node Properties
   /// Node Properties
	std::vector<UInt> NpIds;
	NpIds.push_back(10), NpIds.push_back(11);
      model->CreateNodeProperty ( NpIds[0],NP_STRUC );
      
    NodeProperty* nodeProp = model->CreateNodeProperty ( NpIds[1],NP_GENERIC );
    nodeProp->InitializeProperty(NP_LS,0, numFunctions, NodeIndices, numFunctions, PropFunctionPtrs );


/// Element Properties block
   std::vector<UInt> BlkEpIds;
    BlkEpIds.push_back(20);
    
    ElementProperty* elemProp = model->CreateElementProperty ( BlkEpIds[0],EP_STRUC );
    elemProp->SetDirective( (UInt) ED_STRUC_DAMP_TYPE,DIRECTIVE_FALSE );
    elemProp->InitializeProperty( EP_STRUC_ALPHAD, 0.0, 1, NULL, 1, NULL );
    elemProp->InitializeProperty( EP_STRUC_BETAD , 0.0, 1, NULL, 1, NULL );
    elemProp->InitializeProperty( EP_THICKNESS   , 1.0, 1, NULL, 1, NULL );

	elemProp->InitializeProperty( EP_INTERFACE_PARAM, 1.0e+04, 1, NULL, 1, NULL );
    elemProp->SetDirective( (UInt) ED_ADD_SPRINGS,DIRECTIVE_TRUE );
	elemProp->InitializeProperty( EP_SPRINGS_PARAM, 1.0e-06, 1, NULL, 1, NULL );

//   std::vector<UInt> SSetEpIds;
//    SSetEpIds.push_back(21);
//    elemProp = model->CreateElementProperty ( SSetEpIds[0],EP_STRUCSURF );
//    elemProp->InitializeProperty( EP_STRUC_YTRAC  ,-2.0e-2 , 1, NULL, 1, NULL );
//	elemProp->InitializeProperty( EP_INTERFACE_PARAM, 1.0e+04, 1, NULL, 1, NULL );

/// Materials Block
    std::vector<UInt> BlkMatIds;
    BlkMatIds.push_back(30);
    Material* material = model->CreateMaterial ( BlkMatIds[0], MATERIAL_STRUC_ELASTIC_ISOTROPIC);
    material->InitializeProperty( MP_THERM_RHO, 1.0, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_STRUC_EX, 1.0, 1, NULL, 1, NULL );
    material->InitializeProperty( MP_STRUC_PRXY, 0.2, 1, NULL, 1, NULL );


    /// Nodes
    if (sXfemBlockElemType == NON_XFEM_ELEMENT)
        parseModelHelper.CreateAllNodes(model,NpIds);
    else
        parseModelHelper.CreateAllNodes(model,NpIds);

/// Elements Block
    enum ElementType BlockElemTypes[] = {STRUC_LIN_QUAD4,GNRC_DIFF};
    parseModelHelper.CreateElementsByBaseSetId(model,sBSetIds[0],NON_COUPLED_ELEMENT,NON_IC_ELEMENT,sXfemBlockElemType,BlockElemTypes,BlkEpIds,BlkMatIds);

/// SideSet Elements
//    enum ElementType SSetElemType = STRUC_SURF_BAR2;
//    parseModelHelper.CreateElementsBySideSetId(model,sSSetIds[sLoadNSSetIndx],NON_COUPLED_ELEMENT,NON_IC_ELEMENT,NON_XFEM_ELEMENT,&SSetElemType,SSetEpIds,BlkMatIds);


/// Finalize Model
    parseModelHelper.FinalizeModel(model,STATIC_XFEM_INTERFACE);


/// MPC Block
    model->SetMpcDofs();


/// Dirichlet Conditions
    UInt phase1Ind = 0;
    UInt phase2Ind = 1;
    UInt phaseEnrichInd = -1;
    
    parseModelHelper.ApplyNodeConditionToSSetByIdXFem( model, sSSetIds[sAnchNSSetIndx], DIRICHLET, REGULAR_DOF, UX, phase1Ind, phaseEnrichInd, 0.0);
    parseModelHelper.ApplyNodeConditionToSSetByIdXFem( model, sSSetIds[sAnchNSSetIndx], DIRICHLET, REGULAR_DOF, UY, phase1Ind, phaseEnrichInd, 0.0);
    
    parseModelHelper.ApplyNodeConditionToAllXFem( model, DIRICHLET, REGULAR_DOF, UX, phase2Ind, phaseEnrichInd, 0.0);
    parseModelHelper.ApplyNodeConditionToAllXFem( model, DIRICHLET, REGULAR_DOF, UY, phase2Ind, phaseEnrichInd, 0.0);
    
	parseModelHelper.ApplyNodeConditionToSSetByIdXFem( model, sSSetIds[sLoadNSSetIndx], NEUMANN, REGULAR_DOF, UY, phase1Ind, phaseEnrichInd, -2.0e-2 );

    parseModelHelper.ApplyNodeConditionToAllXFem( model, IDISP, REGULAR_DOF, UX, phase1Ind, phaseEnrichInd, 0.0);
    parseModelHelper.ApplyNodeConditionToAllXFem( model, IDISP, REGULAR_DOF, UY, phase1Ind, phaseEnrichInd, 0.0);
    parseModelHelper.ApplyNodeConditionToAllXFem( model, IDISP, REGULAR_DOF, UX, phase2Ind, phaseEnrichInd, 0.0);
    parseModelHelper.ApplyNodeConditionToAllXFem( model, IDISP, REGULAR_DOF, UY, phase2Ind, phaseEnrichInd, 0.0);
    
    /*
    model->CreateXFemDirichletCondition(UX , PhaseIds[0], sSSetIds[sAnchNSSetIndx], SIDESET_VAR, 0.0);
    model->CreateXFemDirichletCondition(UY , PhaseIds[0], sSSetIds[sAnchNSSetIndx], SIDESET_VAR, 0.0);
    model->CreateXFemDirichletCondition(UX , PhaseIds[1], 0                       , GLOBAL_VAR , 0.0);
    model->CreateXFemDirichletCondition(UY , PhaseIds[1], 0                       , GLOBAL_VAR , 0.0);


/// Initial Conditions Block
    model->CreateXFemInitialCondition(UX   , PhaseIds[0], 0.0);
    model->CreateXFemInitialCondition(UX   , PhaseIds[1], 0.0);
    model->CreateXFemInitialCondition(UY   , PhaseIds[0], 0.0);
    model->CreateXFemInitialCondition(UY   , PhaseIds[1], 0.0);*/
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
    modelData->mStrucModelData        = new StrucModelData();
    SolutionData* solutionData        = new SolutionData();
    modelData->mSolutionData.push_back( solutionData );


///Set Problem and Mesh
    meshData->SetMeshId ( 1 );
    meshData->mMeshFileData->mMeshFileType = MESH_EXODUS;
    const char* InputMeshName             = "beam2D.g";
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
    solverData->mTimeSolverData->mBdfOrder       = 2;
    solverData->mTimeSolverData->mTime           = 0.0;
    solverData->mTimeSolverData->mTimeIts        = sNumIts;
    solverData->mTimeSolverData->mTimeItList[0]  = solverData->mTimeSolverData->mTimeIts;
    solverData->mTimeSolverData->mDt             = sTimeStep;
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
    solverData->mTimeSolverData->mNaNCheck = false;
    
        solverData->mTimeSolverData->mReinitByStandaloneFile = true;
    const char* tempPath = "tempDir/";
    solverData->mTimeSolverData->mSolVecDir = new char[200];
    strcpy(solverData->mTimeSolverData->mSolVecDir,tempPath);


/// Newton solver and linear solver options
    solverData->mNonlinearSolverDataList[0]->mType             = NONLINEAR_SOLVER_FEMDOC_NEWTON;
    solverData->mNonlinearSolverDataList[0]->mMaxItsList[0]    = 10;
    solverData->mNonlinearSolverDataList[0]->mRelaxation       = 1.0;
    solverData->mNonlinearSolverDataList[0]->mRlxList[0]       = solverData->mNonlinearSolverDataList[0]->mRelaxation;
    solverData->mNonlinearSolverDataList[0]->mTotResNormDrop   = 1.0e-3;
    solverData->mNonlinearSolverDataList[0]->mTotResNorm       = 1.0e-20;

/// Linear solver options
	solverData->mLinearSolverDataList[0]->mSparseMatrixType    = SPARSE_MATRIX_TRILINOS_EPETRA_FECRS;

//	solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AMESOS;
//#ifdef PARALLEL
//	solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_MUMPS;
//#else
//	solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_UMFPACK;
//#endif
//	solverData->mLinearSolverDataList[0]->mExport          = LINSOL_EXPORT_MATLAB;

	solverData->mLinearSolverDataList[0]->mType             = LINSOL_TRILINOS_AZTEC;
	solverData->mLinearSolverDataList[0]->mMethod           = LINSOL_AZTEC_GMRES;
	solverData->mLinearSolverDataList[0]->mPreconditioner   = LINSOL_PREC_AZTEC_DD_ILUT;
	solverData->mLinearSolverDataList[0]->mKrylovSpace      = 500;
	solverData->mLinearSolverDataList[0]->mMaxIts           = 3000;
	solverData->mLinearSolverDataList[0]->mDropTol          = 1.0e-09;
	solverData->mLinearSolverDataList[0]->mGraphFill        = 1;
	solverData->mLinearSolverDataList[0]->mIlutFill         = 7.0;
	solverData->mLinearSolverDataList[0]->mPolyOrder        = 3;
	solverData->mLinearSolverDataList[0]->mAlphaThresh      = 0.0;
	solverData->mLinearSolverDataList[0]->mRhoThresh        = 0.0;
	solverData->mLinearSolverDataList[0]->mEpsilon          = 1.0e-09;
	solverData->mLinearSolverDataList[0]->mOverlap          = 1;
	solverData->mLinearSolverDataList[0]->mReorder          = 1;

	solverData->mLinearSolverDataList[0]->mHardBreak        = false;
	solverData->mLinearSolverDataList[0]->mOutput           = LINSOL_OUTPUT_ALL;



/// Model Data
    modelData->SetModelId( 1 );
    modelData->mMeshId                        = 1;
    modelData->mSolverId                      = 1;
    modelData->mModelFileData->mModelFileType = MODEL_FEMDOC;
    modelData->mStrucModelData->mRedDimType   = PLANE_STRAIN;
//    strcpy(modelData->mModelFileData->mIniLevelSetFile,"./IniLSStrucBeam.so");

	modelData->mXFemCurvatureSearchRadius = sCurvatureSearchRadius;
    modelData->mdScalar_ds_byAdv = false;

/// Solution writing format
    solutionData->mSolnFormatType = SF_EXODUS;
    const char* OutputMeshBase   = "beam2D.e-s";
    solutionData->mSaveIC         = sSaveIC;
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
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_UX               , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_UY               , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_UX_XFEM          , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_UY_XFEM          , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_NP_LEVELSET      , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SXX        , ELEM_VAR, 1, sBSetIds  );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SYY        , ELEM_VAR, 1, sBSetIds  );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SZZ        , ELEM_VAR, 1, sBSetIds  );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SXY        , ELEM_VAR, 1, sBSetIds  );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SHY        , ELEM_VAR, 1, sBSetIds  );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SVM        , ELEM_VAR, 1, sBSetIds  );

	parseInputHelper.CreateNodalOutputVariables( solutionData,                          SOL_OBJ_SENS, ELEM_VAR, 1, sBSetIds );
	parseInputHelper.CreateNodalOutputVariables( solutionData,                      SOL_CONST_1_SENS, ELEM_VAR, 1, sBSetIds );
	parseInputHelper.CreateNodalOutputVariables( solutionData, SOL_INTERFACE_CURVATURE_BEAM_LAGRANGE, ELEM_VAR, 1, sBSetIds );

//    for (UInt ie = 0; ie <(sNumDofsPerNodeMax*MAX_NUM_ENRICHMENTS); ++ie)
//    for (UInt ie = 0; ie <(sNumDofsPerNodeMax*MAX_NUM_ENRICHMENTS_2D); ++ie)  // NOTE:: Only use if code was compiled with the hack MAX_NUM_ENRICHMENTS=MAX_NUM_ENRICHMENTS_2D
//        parseInputHelper.CreateNodalOutputVariables   (solutionData, (SolVarType) (SOL_XFEM000+ie), ELEM_VAR, 1, sBSetIds);

    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_EP_MAIN_PHASE    ,           1, sBSetIds );
    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_STRUC_STEN       ,           1, sBSetIds );

    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_VOLUME_PHASE_1, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_VOLUME_PHASE_2, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,  GSOL_STRAIN_ENERGY, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,      GSOL_PERIMETER, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_CURVATURE_TYPE, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,             GSOL_UX, ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,             GSOL_UY, ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,          GSOL_UX_SQ, ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData,          GSOL_UY_SQ, ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));

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
	void ComputeObjectivesAndConstraints(
			const std::vector <Real> & aCriteriaVals,
			Real & aObjectiveVal,
			std::vector <Real> & aConstraintVals )
{
//	aObjectiveVal = sPerPen*aCriteriaVals[0] + sCurPen*aCriteriaVals[1] + sDispPen * aCriteriaVals[4];
//
//	Real totalVol = (aCriteriaVals[2]+aCriteriaVals[3]);
//	aConstraintVals[0] = aCriteriaVals[2]/totalVol/sMassRatio - 1.0;
}

extern "C"
	void Compute_dOptScalarsDQs(
			const std::vector< Real > & aCriteriaVals,
			std::vector< Real > & aGrad )
{
//	UInt numCriteria   = 5;
//
//	aGrad[0*numCriteria + 0] = sPerPen;
//	aGrad[0*numCriteria + 1] = sCurPen;
//	aGrad[0*numCriteria + 2] = 0.0;
//	aGrad[0*numCriteria + 3] = 0.0;
//	aGrad[0*numCriteria + 4] = sDispPen;
//
//	Real totalVol = (aCriteriaVals[2]+aCriteriaVals[3]);
//
//	aGrad[1*numCriteria + 0] = 0.0;
//	aGrad[1*numCriteria + 1] = 0.0;
//	aGrad[1*numCriteria + 2] = (1./totalVol-aCriteriaVals[2]/totalVol/totalVol)/sMassRatio;
//	aGrad[1*numCriteria + 3] = -aCriteriaVals[2]/totalVol/totalVol/sMassRatio;
//	aGrad[1*numCriteria + 4] = 0.0;
}


///----------------------------------------------------------------------------------------
/// function for defining the paramters of the optimization
///----------------------------------------------------------------------------------------
extern "C"
	void SetOptimizationParameters(
			Optimizer * aOptimizer,
			Problem   * aProblem )
{
//	// Finite difference checker and dR_ds perturbation value
//	aOptimizer->SetAdvEpsilon( sOptAdvEpsilon );
//
//	// Set whether or not to use "wiggling" (true = don't wiggle)
//	aOptimizer->SetMapFlag( true );
//
//	// Set the number of constraints
//	aOptimizer->SetNumConstraints(1);
//
//	// Set the optimizer and algorithmic parameters
//	OptimizerData* optSolData = aOptimizer->GetOptimizerData();
//	optSolData->mType     = OPT_SOLVER_OLD_FEM_GCMMA;
//	optSolData->mMaxIts   = sNumOptIts;
//	optSolData->mNormDrop = 1e-6;
//
//	optSolData->mGcmmaData.AsympAdapt0 = 0.5;
//	optSolData->mGcmmaData.AsympAdapt  = 0.7;
//	optSolData->mGcmmaData.AsympAdaptC = 1./0.7;
//	optSolData->mGcmmaData.StepSize    = 0.01;
//	optSolData->mGcmmaData.Penalty     = 100.0;
//
//	optSolData->mDumpAdvs = ADV_EVERY;
//	optSolData->mRestartIteration = 0;
//	aOptimizer->mWriteOptScalarValuesToFile = true;
//
//	/// Finite difference check.
//
//	std::vector<UInt> VarId;
//	std::vector<UInt> FDCheck;
//	std::vector<Real> FDEpsilon;
//
//	// Mega hack
//	std::vector< UInt > varIDs;
//	varIDs.push_back( 547 );
//
//	UInt numVar = varIDs.size();
//
//	printf( "Number of variables: %i\n", numVar );
//
//	for (UInt in = 0; in < numVar; in++)
//	{
//		printf( "Performing sensitivities on Node Id: %i\n", varIDs[in] );
//		VarId.push_back(varIDs[in]);
//		FDCheck.push_back( sOptFDCheck );
//		FDEpsilon.push_back( sOptFDEpsilon );
//	}
//
//	optSolData->mSweepData.mNumVar = numVar;
//	optSolData->mSweepData.SetVarIdVector( VarId );
//	optSolData->mSweepData.SetFDCheckVector( FDCheck );
//	optSolData->mSweepData.SetFDEpsilonVector( FDEpsilon );
//
//	optSolData->mFDCheck = sOptFDCheck;
//	optSolData->mFDEps = sOptFDEpsilon;
//	optSolData->mFDSens = false;

	// Set the optimizer and algorithmic parameters
	OptimizerData* optSolData = aOptimizer->GetOptimizerData();
	optSolData->mType     = OPT_SOLVER_SWEEP;

	static const UInt numVariables		  = 1;					// Number of variables to be analyzed
	static char mSweepFileName[]		  = "sweepResults.txt";
	optSolData->mSweepData.mNumVar		  = numVariables;
	optSolData->mSweepData.mSweepFileName = mSweepFileName;

	std::vector<UInt> NumSteps  ( 1,  sRadiusSteps  );	// This will equal to 500 different meshes, since it starts counting from zero
	std::vector<UInt> VarId     ( 1,  0    );
	std::vector<Real> VarLow    ( 1,  sRadiusLower );
	std::vector<Real> VarUp     ( 1,  sRadiusUpper );
	std::vector<UInt> GradFlag  ( 1,  0    );
	std::vector<UInt> FDCheck   ( 1,  0    );
	std::vector<Real> FDEpsilon ( 1,  0    );

	optSolData->mSweepData.SetNumStepsVector(NumSteps);
	optSolData->mSweepData.SetVarIdVector(VarId);
	optSolData->mSweepData.SetVarLowVector(VarLow);
	optSolData->mSweepData.SetVarUpVector(VarUp);
	optSolData->mSweepData.SetGradFlagVector(GradFlag);
	optSolData->mSweepData.SetFDCheckVector(FDCheck);
	optSolData->mSweepData.SetFDEpsilonVector(FDEpsilon);

	optSolData->mDumpAdvs = ADV_EVERY;
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
//	OptimizationCriteria* crit = NULL;
//	Function* CritEvalFunc = NULL;
//	Real weights = 1.;
//	UInt entries[] = {1,0};
//	UInt indices = sNumIts; // Last time step
//
//	// First bool is dependent on states. Second bool dependent on Adv's
//	crit = aOptimizer->CreateCriteria(GSOL_PERIMETER,false,true);
//	CritEvalFunc = new FunctionScalarProduct(entries,&weights);
//	crit->CreateEvaluationFunction(CritEvalFunc,&indices);
//
//	crit = aOptimizer->CreateCriteria(GSOL_CURVATURE_TYPE,false,true);
//	CritEvalFunc = new FunctionScalarProduct(entries,&weights);
//	crit->CreateEvaluationFunction(CritEvalFunc,&indices);
//
//	// First bool is dependent on states. Second bool dependent on Adv's
//	crit = aOptimizer->CreateCriteria(GSOL_VOLUME_PHASE_1,false,true);
//	CritEvalFunc = new FunctionScalarProduct(entries,&weights);
//	crit->CreateEvaluationFunction(CritEvalFunc,&indices);
//
//	// First bool is dependent on states. Second bool dependent on Adv's
//	crit = aOptimizer->CreateCriteria(GSOL_VOLUME_PHASE_2,false,true);
//	CritEvalFunc = new FunctionScalarProduct(entries,&weights);
//	crit->CreateEvaluationFunction(CritEvalFunc,&indices);
//
//	// First bool is dependent on states. Second bool dependent on Adv's
//	crit = aOptimizer->CreateCriteria(GSOL_STRAIN_ENERGY,true,true);
//	CritEvalFunc = new FunctionScalarProduct(entries,&weights);
//	crit->CreateEvaluationFunction(CritEvalFunc,&indices);
}

///----------------------------------------------------------------------------------------
/// function for defining the Physical design variables and Abstract design variables
///----------------------------------------------------------------------------------------
extern "C"
	void SetOptimizationPdvsAndAdvs(
			Optimizer * aOptimizer,
			Problem   * aProblem )
{
//	UInt modelId = 1;
//	Model* model = aProblem->GetModelById( 1 );
//	Mesh* mesh = model->GetMesh();
//
//	model->BuildNodeToElementTable();
//
//	// Vars for creating PDVs
//	DesignDepObject desDepObj = NOD_PROP;
//	UInt propId = 11;
//	int propType = NP_LS;
//
//	UInt myNumNodesInSet = model->GetNumNodes();
//
//	// ------------------------------------------------------------------------
//	std::vector < UInt > &DesVarKeys = aOptimizer->GetMyDesVarKeys();
//
//	DesVarKeys.assign( myNumNodesInSet, 0 );
//
//	for ( UInt in = 0; in < myNumNodesInSet; ++in )
//		DesVarKeys[in] = model->GetNode( in )->GetIdNum();
//
//	// Set the number of design variables per processor
//	aOptimizer->FinalizeDesVarCreation();
//
//	// ------------------------------------------------------------------------
//	// Create PDVs
//	std::vector < UInt > globalNodesIds( myNumNodesInSet, 0 );
//	for ( UInt in = 0; in < myNumNodesInSet; ++in )
//		globalNodesIds[in] = model->GetNode( in )->GetIdNum();
//
//	UInt *ElemIndices = NULL, *ElemIds = NULL;
//	UInt numLoadElems = 0;
//
//	// Initialize a ParseModelHelper
//	ParseModelHelper parseModelHelper = ParseModelHelper();
//
//	// Kill filter in inflow and outflow.
//	std::vector< UInt > noFilterNodeIds;
//
//	numLoadElems = parseModelHelper.GetSideSetNeighborElemsInfoByID( mesh, sSSetIds[sLoadNSSetIndx], ElemIndices, ElemIds );
//	for ( UInt ie = 0; ie < numLoadElems; ++ie ) {
//		std::vector < UInt > NodeIdsVec( FP_MAX_NODES, 0 );
//		Element* elem = model->GetElement( ElemIndices[ie] );
//		UInt numNodes = elem->GetNodeIds( &NodeIdsVec[0] );
//
//		for ( UInt in = 0; in < numNodes; ++in ) {
//			noFilterNodeIds.push_back( NodeIdsVec[in] );
//		}
//	}
//
//	FemdocAlgorithms::SortUniqueResizeStdVector( noFilterNodeIds );
//
//	// Read file
//	std::string line;
//	std::ifstream myfile ( "disttable" );
//
//	if ( myfile.is_open() )
//	{
//		std::fprintf( stdout, "\n ... Reading PDVs file\n" );
//		while ( myfile.good() )
//		{
//			getline( myfile, line );
//
//			//          char* lineChar = new char[line.size() + 1];
//			//          lineChar[line.size()] = 0;
//			//          memcpy( lineChar, line.c_str(), line.size() );
//
//			const char* lineChar = line.c_str();
//
//			Real test1 = 0.0;
//			Real test2 = 0.0;
//
//			std::sscanf( lineChar, "%lf %lf",  &test1, &test2 );
//			UInt glbNodeId = test1;
//			UInt numNeigb = test2;
//
//			assert( numNeigb > 0 );
//
//			if ( std::find( &globalNodesIds[0], &globalNodesIds[0] + myNumNodesInSet, glbNodeId ) - &globalNodesIds[0] != myNumNodesInSet )
//			{
//				Node* node = model->GetNodeById( glbNodeId );
//
//				std::vector < UInt > neigbIdsVec ( ( UInt ) numNeigb , 0.0 );
//				std::vector < Real > neigbWgtVec ( ( UInt ) numNeigb , 0.0 );
//
//				if ( node )
//				{
//					// If found in no filter IDs.
//					if ( std::find( &noFilterNodeIds[0], &noFilterNodeIds[0] + ( UInt ) noFilterNodeIds.size(), glbNodeId ) - &noFilterNodeIds[0] != ( UInt ) noFilterNodeIds.size() ) {
//						UInt numEntries[] = {1,0};
//
//						neigbIdsVec[0] = glbNodeId - 1;
//						neigbWgtVec[0] = 1.0;
//
//						UInt paramIndex = 0;
//						UInt propTypeKey = node->GetIndexNum();
//						PhysicalDesignVar* pdv = aOptimizer->CreatePhysDesVar( aProblem, modelId, desDepObj, propId, propType, propTypeKey, paramIndex );
//
//						Function* PDVEvalFunc = new FunctionScalarProduct( numEntries, &neigbWgtVec[0] );
//						pdv->CreateEvaluationFunction( PDVEvalFunc, &neigbIdsVec[0] );
//					}
//					else {
//						Real neighborIds = 0;
//						Real neighborWgt = 0;
//						UInt position = 28;
//						for ( UInt ii = 0; ii < numNeigb; ii++ )
//						{
//							std::sscanf( lineChar + position,"%lf",&neighborIds );
//							neigbIdsVec[ii] = neighborIds - 1;
//							position += 14;
//						}
//						for ( UInt ii = 0; ii < numNeigb; ii++ )
//						{
//							std::sscanf( lineChar + position,"%lf",&neighborWgt );
//							neigbWgtVec[ii] = neighborWgt;
//							position += 14;
//						}
//
//						UInt numEntries[] = {numNeigb, 0};
//						UInt paramIndex = 0;
//						UInt propTypeKey = node->GetIndexNum();
//						PhysicalDesignVar* pdv = aOptimizer->CreatePhysDesVar( aProblem, modelId, desDepObj, propId, propType, propTypeKey, paramIndex );
//
//						Function* PDVEvalFunc = new FunctionScalarProduct( numEntries, &neigbWgtVec[0] );
//						pdv->CreateEvaluationFunction( PDVEvalFunc, &neigbIdsVec[0] );
//					}
//				}
//			}
//		}
//	}
//	else
//	{
//		std::fprintf( stdout, " Unable to open file!!!\n" );
//	}

	// Set number of abstract design variables (radius)
	UInt modelId = 1;
	UInt numAbsDesVar   = 1;
	UInt AbsDesVarInd[] = {0};
	// ------------------------------------------------------------------------
	// Set number of abstract design variables
    aOptimizer->SetNumDesVars(numAbsDesVar);

    // ------------------------------------------------------------------------
    // Create PDV's
    DesignDepObject desDepObj = NOD_PROP;
    UInt propId   = 11;
    Int  propType = NP_LS;
    PhysicalDesignVar* pdv = NULL;

    // Create optimizer sweep function
	UInt numParameters[] = { 3, 1 };      // number of coefficients and independent variables
	Real ParameterVals[] = { 0.0, 0.0, 0.0 };  // radius and center of circle

	Real CntrXGlbCoords = sCenterX;
	Real CntrYGlbCoords = sCenterY;
	Real CntrZGlbCoords = sCenterZ;

	Model* model = aProblem->GetModelById(modelId);
	UInt numNodes = model->GetNumNodes();

	model->BuildNodeToElementTable();

	Node* node = NULL;
	Real NodeXGlbCoords = 0.0;
	Real NodeYGlbCoords = 0.0;
	Real NodeZGlbCoords = 0.0;

	for ( UInt in = 0; in < numNodes; ++in )
	{
		node = model->GetNode(in);
		NodeXGlbCoords = node->GetGlobalXCoord();
		NodeYGlbCoords = node->GetGlobalYCoord();
		NodeZGlbCoords = node->GetGlobalZCoord();

		ParameterVals[0] = NodeXGlbCoords-CntrXGlbCoords;
		ParameterVals[1] = NodeYGlbCoords-CntrYGlbCoords;
		ParameterVals[2] = NodeZGlbCoords-CntrZGlbCoords;

		// Parameters for moving circle will be [distance, xcoords of node, ycoords of node, zcoords of node]

		pdv = aOptimizer->CreatePhysDesVar(aProblem,modelId,desDepObj,propId,propType,node->GetIndexNum());
		Function* SweepFunc = new FunctionUserDefined(numParameters,ParameterVals,mySweepFunc,mySweepFuncDeriv);
		pdv->CreateEvaluationFunction(SweepFunc,AbsDesVarInd);
	}
}


///----------------------------------------------------------------------------------------
/// Initialize the Vector of abstract design variables (Adv) mAbsDesVarVec
///----------------------------------------------------------------------------------------
extern "C"
    void InitializeAbsDesVarVector(Optimizer* aOptimizer, UInt aNumMyAbsDesVars)
{
//	// Get abstract design variables vector
//	Vector * absDesVarVec = aOptimizer->GetAbsDesVarVec();
//
//	// Get mesh
//	Problem* problem = aOptimizer->GetProblem();
//
//	// Get the model
//	Model* model = problem->GetModelById( 1 );
//	Mesh* mesh = model->GetMesh();
//
//	// Initialize a ParseModelHelper
//	ParseModelHelper parseModelHelper = ParseModelHelper();
//
//	// ------------------------------------------------------------------------
//	// Generate swiss cheese
//	UInt numNodes = model->GetNumNodes();
//
//	Real InclusionRadius = 0.25;
////	Real InclusionRadius = 0.50;
//	InclusionRadius= InclusionRadius;
//
//	Real exponent = 2.0;
////	exponent = 100.0;
//	exponent = exponent;
//
//	bool restartFromFile = false;
//	if ( restartFromFile )
//	{
//		// Get all nodes information
//		ParseModelHelper parseModelHelper = ParseModelHelper();
//		UInt *NodeIndices = NULL, *NodeIds = NULL;
//		UInt  numAllNodes = parseModelHelper.GetAllNodesInfo(mesh,NodeIndices,NodeIds);
//		UInt  totNumNodes = model->GetGlbNumNodes();
//
//		// allocate memory space for storing variables read from file
//		std::vector<double> Advs(totNumNodes,0.0);
//
//		// read otimization variables from restart file
//		char buffer[100];
//		//int iteration = 199;
//		std::sprintf(buffer,"./AbsDesVariablesXFEM.dat");
//		FILE* restartFile = std::fopen(buffer,"rb");
//		if ( restartFile )
//		{
//			std::fread(&Advs[0],sizeof(double),totNumNodes,restartFile);
//			fclose(restartFile);
//
////			// ------------------------------------------------------------------------
////			// Scale ADVs
////			Real maxElement = *std::max_element( &Advs[0], &Advs[0] + numNodes );
////			Real minElement = *std::min_element( &Advs[0], &Advs[0] + numNodes );
////
////			Real glbMaxElement = maxElement;
////			Real glbMinElement = minElement;
////		#ifdef PARALLEL
////			MPI_Allreduce( &maxElement, &glbMaxElement, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
////			MPI_Allreduce( &minElement, &glbMinElement, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
////		#endif
////
////			for ( UInt nInd = 0; nInd < numNodes; ++ nInd ) {
////				Advs[nInd] = ( ( 1 - ( -1 ) ) * ( Advs[nInd] - glbMinElement ) ) / ( glbMaxElement - glbMinElement ) + ( -1 );
////			}
//
//			for (UInt inod=0, nodInd=0; inod<numAllNodes; ++inod)
//			{
//				nodInd = NodeIds[inod]-1;
//				if ( model->GetNode(inod)->IsMaster() )
//				{
//					absDesVarVec->ReplaceGlobalValues(1,&nodInd,&Advs[nodInd]);
//				}
//			}
//		}
//	}
//	else {
//		std::vector < Real > Advs( numNodes, 0.0 );
//		std::vector < Real > DesignOptimizationDomain;
//		std::vector < UInt > NumberInclusions( MAX_ELEMENT_DIMS, 0 );
//
//		Real InclusionsXYZ[3] = {4,3,1};
//
//		for ( UInt idim = 0; idim < MAX_ELEMENT_DIMS; ++idim )
//			NumberInclusions[idim] = InclusionsXYZ[idim];
//
//		std::vector < bool > OffSet ( MAX_ELEMENT_DIMS, true );
//
//		mesh->GenerateCirclesOrSquares( Advs, OffSet, DesignOptimizationDomain, NumberInclusions, InclusionRadius );
//
//		// Change
////		for ( UInt in = 0; in < numNodes; ++in ) {
////			if ( model->GetNode( in )->GetGlobalCoords()[1] < 0.95 ) {
////				Advs[in] = +1.0;
////			}
////			else {
////				Advs[in] = -1.0;
////			}
////		}
////
////		for ( UInt in = 0; in < numNodes; ++in ) {
////			Advs[in] = InclusionRadius - std::pow( std::pow( model->GetNode( in )->GetGlobalCoords()[0] - 1.50, exponent ) + std::pow( model->GetNode( in )->GetGlobalCoords()[1] - 1.00, exponent ), 1.0 / exponent );
////		}
////
////		for ( UInt in = 0; in < numNodes; ++in ) {
////			Advs[in] = std::max( Advs[in], InclusionRadius - std::sqrt( std::pow( model->GetNode( in )->GetGlobalCoords()[0] - 2.5, 2.0 ) + std::pow( model->GetNode( in )->GetGlobalCoords()[1], 2.0 ) ) );
////		}
////
////		for ( UInt in = 0; in < numNodes; ++in ) {
////			Advs[in] = std::max( Advs[in], InclusionRadius - std::sqrt( std::pow( model->GetNode( in )->GetGlobalCoords()[0] - 5.0, 2.0 ) + std::pow( model->GetNode( in )->GetGlobalCoords()[1], 2.0 ) ) );
////		}
//
//		// ------------------------------------------------------------------------
//		// Assign ADVs
//		for ( UInt nInd = 0; nInd < numNodes; ++ nInd ) {
//			Node* node = model->GetNode( nInd );
//			if ( node->IsMaster() ) {
//				UInt myAdvInd = node->GetIdNum() - 1;       // This is the warning above
//				if ( std::abs( Advs[nInd] - 0.0 ) < 1.0e-14 ) {
//					printf( " absDesVarVec: %5.8e\n", Advs[nInd] );
//				}
//
//				absDesVarVec->ReplaceGlobalValues( 1, &myAdvInd, NULL, &Advs[nInd] );
//			}
//		}
//	}
//
//	// ------------------------------------------------------------------------
//	// Get max and min value in vector
//	Real maxVal = 0.0, minVal = 0.0;
//	absDesVarVec->MaxValue( &maxVal );
//	absDesVarVec->MinValue( &minVal );
//	Real maxLsVal = std::max( std::abs( maxVal ), std::abs( minVal ) );
//
//	Real glbMaxLsVal = maxLsVal;
//#ifdef PARALLEL
//	MPI_Allreduce( &maxLsVal, &glbMaxLsVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
//#endif
//	Real AdvLowerBound = -1.2*glbMaxLsVal;
//	Real AdvUpperBound =  1.2*glbMaxLsVal;
//
//	for ( UInt nInd = 0; nInd < numNodes; ++ nInd )
//	{
//		Node* node = model->GetNode( nInd );
//		if ( node->IsMaster() )
//		{
//			UInt AdvIndex = node->GetIdNum() - 1;
//			aOptimizer->SetUpperAndLowerBounds( 1, &AdvIndex, &AdvLowerBound, 1, &AdvIndex, &AdvUpperBound );
//		}
//	}
//
//	// ------------------------------------------------------------------------
//	/// Find and fix design variables at load points
//	// Get sideset nodes information
//	Real fixVal = -1.0;
//	Real fixValMin = fixVal;
//	Real fixValMax = fixVal + 1.0e-04;
//
//	fixValMin = fixValMin; fixValMax = fixValMax;
//
//	UInt *ElemIndices = NULL, *ElemIds = NULL;
//	UInt numLoadElems = 0;
//
//	// Constrain abstract design variables for the inflow
//	for ( UInt in = 0; in < numNodes; ++in ) {
//		Node * node = model->GetNode( in );
//		Real * glbCoords = node->GetGlobalCoords();
//
//		UInt nodeIndex = node->GetIdNum() - 1;
//		if ( glbCoords[0] > 2.75 && glbCoords[1] > 0.75 && glbCoords[1] < 1.25) {
//			absDesVarVec->ReplaceGlobalValues ( 1, &nodeIndex, &fixVal );
//		}
//	}
//
//	numLoadElems = parseModelHelper.GetSideSetNeighborElemsInfoByID( mesh, sSSetIds[sLoadNSSetIndx], ElemIndices, ElemIds );
//	for ( UInt ie = 0; ie < numLoadElems; ++ie ) {
//
//		std::vector < UInt > NodeIdsVec( FP_MAX_NODES, 0 );
//		Element* elem = model->GetElement( ElemIndices[ie] );
//		UInt numNodes = elem->GetNodeIds( &NodeIdsVec[0] );
//
//		for ( UInt in = 0; in < numNodes; ++in ) {
//
//			Node* node = model->GetNodeById( NodeIdsVec[in] );
//
//			if ( node->IsMaster() ) {
//				UInt nodeIndex = NodeIdsVec[in] - 1;
//				absDesVarVec->ReplaceGlobalValues ( 1, &nodeIndex, &fixVal );
//				aOptimizer->SetUpperAndLowerBounds( 1, &nodeIndex, &fixValMin, 1, &nodeIndex, &fixValMax );
//			}
//
//		}
//
//	}

    // Get abstract design variables vector
    Vector* absDesVarVec = aOptimizer->GetAbsDesVarVec();

    // set the bounds for mAbsDesVarVec
    UInt indices[]= {0};
    Real lowerBounds[]= {-1.0};
    Real upperBounds[]= { 1.0};
    aOptimizer->SetUpperAndLowerBounds(1,indices,lowerBounds,1,indices,upperBounds);

    Real initVals[]= {0.5};

    absDesVarVec->ReplaceGlobalValues(1,indices,initVals);
}
