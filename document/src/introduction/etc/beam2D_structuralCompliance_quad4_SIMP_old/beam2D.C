#include "DynLinkTools/davesMacros.h"
#include "DynLinkTools/StandardIncludes.h"

/// Element Information Creation Block
static const UInt sBSetIds[] = {1};
static const UInt sSSetIds[] = {1,2,3,4,5};
static const UInt sNSetIds[] = {1,2,3,4,5,6};
static const UInt sAnchNSSetIndx = 3;  // Anchor nodeset&sideset index
static const UInt sLoadNSSetIndx = 4;  // Load surface nodeset&sideset index

static const enum XFemElemType  sXfemBlockElemType  = NON_XFEM_ELEMENT;
static const UInt sNumDofsPerNodeMax                = 2;

/// Time Integration options
static const Real sTimeStep              = 200000.;
static const UInt sNumIts                = 1;
static const Real sTotalTime             = sTimeStep*sNumIts;

/// -------------------------------------------------------------------------------------
///  SECTION FOR OPTIMIZATION INPUTS
/// -------------------------------------------------------------------------------------
static const Real sPerPen    = 1.0e-2;
static const Real sMassRatio = 0.5;
static const Real sSimpExp   = 3.0;

static const bool sFDCheck = false;
static const bool sFDSens = false;
static const Real sFDEps = 1.0e-5;

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

    UInt modelId        = aModelDataStruct->GetModelId();
    UInt meshId         = aModelDataStruct->GetMeshId();
    UInt numSolnWriters = aModelDataStruct->mSolutionData.size();

/// MODEL --------------------------------------------------------
    Model* model = problem->CreateModel ( aModelDataStruct,modelId,meshId,0,0,0,0,
                                          0,0,0,numSolnWriters);
    Mesh* msh = model->GetMesh();


/// Functions
    // Physical design variables function creation
    UInt numE[3] = {0,0,0};
    UInt *IndicesList=NULL, *IdsList=NULL;
    UInt numPropFunctions = parseModelHelper.GetBaseSetElemsInfoByID(msh,sBSetIds[0],IndicesList,IdsList);
    std::vector<UInt>       ElemIndicesVec(IndicesList,IndicesList+numPropFunctions);
    std::vector<PropFunction*>  EmodFunctionPtrsVec(numPropFunctions,NULL), RhoFunctionPtrsVec(numPropFunctions,NULL);
    Real funcVal            = 1.0;
    enum FunctionVarType funcVarType = FUNC_VAR_DESVAR;
    for ( UInt in = 0; in < numPropFunctions; ++in )
    {
        EmodFunctionPtrsVec[in] = model->CreatePropFunction(                 in+1,FUNCTION_CONSTANT,&funcVarType, numE,&funcVal);
        RhoFunctionPtrsVec[in]  = model->CreatePropFunction(numPropFunctions+in+1,FUNCTION_CONSTANT,&funcVarType, numE,&funcVal);
    }

    // Loading function creation
    const UInt numPts1    = 2;
    funcVarType           = FUNC_VAR_TIME;
    Real InterpFuncParams[] = {0.0, 0.0, sTotalTime, 1.0};
    UInt bcFuncId = 1;
    model->CreateBCFunction(bcFuncId,FUNCTION_INTERPOLATION_1D,&funcVarType,&numPts1,InterpFuncParams);


/// Node Properties
	std::vector<UInt> NpIds;
	NpIds.push_back(10);
      model->CreateNodeProperty ( NpIds[0],NP_STRUC );



/// Element Properties block

    std::vector<UInt> BlkEpIds;
    BlkEpIds.push_back(20);
    
    ElementProperty* elemProp = model->CreateElementProperty ( BlkEpIds[0],EP_STRUC );
    elemProp->InitializeProperty( EP_THICKNESS   , 1.0, 1, NULL, 1, NULL );
    std::vector<UInt> SSetEpIds;
    SSetEpIds.push_back(21);

    elemProp = model->CreateElementProperty ( SSetEpIds[0],EP_STRUCSURF );
    elemProp->InitializeProperty( EP_STRUC_YTRAC  ,-2.0e-2 , 1, NULL, 1, NULL );

/// Materials Block
    std::vector<UInt> BlkMatIds;
    BlkMatIds.push_back(30);
    Material* material = model->CreateMaterial ( BlkMatIds[0], MATERIAL_STRUC_ELASTIC_ISOTROPIC);
//    Material* material = model->CreateMaterial ( BlkMatIds[0], MATERIAL_STRUC_HYPERELASTIC_KIRCHHOFF);
    material->InitializeProperty( MP_STRUC_RHO , 1.0, numPropFunctions, &ElemIndicesVec[0], numPropFunctions, &RhoFunctionPtrsVec[0]);
    material->InitializeProperty( MP_STRUC_EX  , 1.0, numPropFunctions, &ElemIndicesVec[0], numPropFunctions, &EmodFunctionPtrsVec[0]);
    material->InitializeProperty( MP_STRUC_PRXY, 0.2, 1, NULL, 1, NULL );

/// Nodes
    parseModelHelper.CreateAllNodes(model,NpIds);


/// Block Elements
//    enum ElementType BlockElemType = STRUC_NLTL_QUAD4;
    enum ElementType BlockElemType = STRUC_LIN_QUAD4;
    parseModelHelper.CreateElementsByBaseSetId(model,sBSetIds[0],NON_COUPLED_ELEMENT,NON_IC_ELEMENT,NON_XFEM_ELEMENT,&BlockElemType,BlkEpIds,BlkMatIds);

/// SideSet Elements
    enum ElementType SSetElemType = STRUC_SURF_BAR2;
//    enum ElementType SSetElemType = STRUC_NLTLS_QUAD4;
    parseModelHelper.CreateElementsBySideSetId(model,sSSetIds[sLoadNSSetIndx],NON_COUPLED_ELEMENT,NON_IC_ELEMENT,NON_XFEM_ELEMENT,&SSetElemType,SSetEpIds,BlkMatIds);


/// Finalize Model
    parseModelHelper.FinalizeModel(model);


///MPC Block
    model->SetMpcDofs();


/// Dirichlet Conditions
    parseModelHelper.ApplyNodeConditionToSSetById( model, sSSetIds[sAnchNSSetIndx], DIRICHLET, REGULAR_DOF, UX, 0.0, 0);
    parseModelHelper.ApplyNodeConditionToSSetById( model, sSSetIds[sAnchNSSetIndx], DIRICHLET, REGULAR_DOF, UY, 0.0, 0);


/// Neumann Conditions
//    parseModelHelper.ApplyNodeConditionToNSetById( model, sNSetIds[sLoadNSSetIndx], NEUMANN, REGULAR_DOF, UY, 1.0, bcFuncId1);


/// Initial Conditions Block
    parseModelHelper.ApplyNodeConditionToAll( model, IDISP, REGULAR_DOF, UX, 0.0, 0);
    parseModelHelper.ApplyNodeConditionToAll( model, IDISP, REGULAR_DOF, UY, 0.0, 0);
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
//    const char* InputMeshName             = "beam2D_30x20.g";
    const char* InputMeshName             = "beam2D_60x40.g";
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
    solverData->mNonlinearSolverDataList[0]->mMaxItsList[0]    = 6;
    solverData->mNonlinearSolverDataList[0]->mRelaxation       = 1.0;
    solverData->mNonlinearSolverDataList[0]->mRlxList[0]       = solverData->mNonlinearSolverDataList[0]->mRelaxation;
    solverData->mNonlinearSolverDataList[0]->mTotResNormDrop   = 1.0e-4;
    solverData->mNonlinearSolverDataList[0]->mTotResNorm       = 1.0e-20;


/// Linear solver options
    solverData->mLinearSolverDataList[0]->mType                = LINSOL_TRILINOS_AMESOS;
    solverData->mLinearSolverDataList[0]->mSparseMatrixType    = SPARSE_MATRIX_TRILINOS_EPETRA_FECRS;
    solverData->mLinearSolverDataList[0]->mMethod              = LINSOL_AMESOS_UMFPACK;
    solverData->mLinearSolverDataList[0]->mPreconditioner      = LINSOL_PREC_AZTEC_DD_ILUT;
//    solverData->mLinearSolverDataList[0]->mExport              = LINSOL_EXPORT_MATLAB;


/// Model Data
    modelData->SetModelId( 1 );
    modelData->mMeshId                        = 1;
    modelData->mSolverId                      = 1;
    modelData->mModelFileData->mModelFileType = MODEL_FEMDOC;
    modelData->mStrucModelData->mRedDimType   = PLANE_STRAIN;


/// Solution writing format
    solutionData->mSolnFormatType = SF_EXODUS;
    const char* OutputMeshBase   = "beam2D.e-s";
    solutionData->mSaveIC         = 0;
    solutionData->mOptSaveFreq    = 1;
    solutionData->mSaveOpt        = true;
    solutionData->mIsXFem         = (sXfemBlockElemType!=NON_XFEM_ELEMENT);
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
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SXX        , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SYY        , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SXY        , ELEM_VAR, 1, sBSetIds );
    parseInputHelper.CreateNodalOutputVariables   (solutionData, SOL_STRUC_SVM        , ELEM_VAR, 1, sBSetIds );

    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_MP_STRUC_EX      ,           1, sBSetIds );
    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_MP_STRUC_RHO     ,           1, sBSetIds );
    parseInputHelper.CreateBlockSetOutputVariables(solutionData, SOL_STRUC_STEN       ,           1, sBSetIds );

    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_STRAIN_ENERGY   , ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_VOLUME          , ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_MASS_UTILIZATION, ELEM_VAR, 1, sBSetIds , 0, NULL);
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_UX              , ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_UY              , ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_UX_SQ           , ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));
    parseInputHelper.CreateGlobalOutputVariables  (solutionData, GSOL_UY_SQ           , ELEM_VAR, 1, sBSetIds , 1, (sNSetIds+sLoadNSSetIndx));

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
    aObjectiveVal = 0.0*aCriteriaVals[0] + 100.0*std::pow((aCriteriaVals[2] + aCriteriaVals[3]),0.5);

    /// ----------------------------------------------------------------------------
    /// Constraint is a limiting volume of phase 1
    /// ----------------------------------------------------------------------------
    aConstraintVals[0] = aCriteriaVals[1]/sMassRatio - 1.;

}

extern "C"
    void Compute_dOptScalarsDQs(const std::vector< Real > & aCriteriaVals, std::vector< Real > & aGrad)
{
    UInt numCriteria   = 4;
//    UInt numOptScalars = 2;

    aGrad[0*numCriteria + 0] = 0.;
    aGrad[0*numCriteria + 1] = 0.;
    aGrad[0*numCriteria + 2] = 100.0*0.5/std::sqrt((aCriteriaVals[2] + aCriteriaVals[3]));
    aGrad[0*numCriteria + 3] = 100.0*0.5/std::sqrt((aCriteriaVals[2] + aCriteriaVals[3]));

    aGrad[1*numCriteria + 0] = 0.;
    aGrad[1*numCriteria + 1] = 1./sMassRatio;
    aGrad[1*numCriteria + 2] = 0.;
    aGrad[1*numCriteria + 3] = 0.;

}

///----------------------------------------------------------------------------------------
/// function for defining the paramters of the optimization
///----------------------------------------------------------------------------------------
extern "C"
    void SetOptimizationParameters(Optimizer *aOptimizer,
                                   Problem   *aProblem)
{
    // Finite difference checker and dR_ds perturbation value
    aOptimizer->SetAdvEpsilon( 1.0e-5);

    // Set whether or not to use "wiggling" (true = don't wiggle)
    aOptimizer->SetMapFlag(true);

    // Set the number of constraints
    aOptimizer->SetNumConstraints(1);

    // Set the optimizer and algorithmic parameters
    OptimizerData* optSolData = aOptimizer->GetOptimizerData();
    optSolData->mType     = OPT_SOLVER_OLD_FEM_GCMMA;
    optSolData->mMaxIts   = 200;
    optSolData->mNormDrop = 1e-6;

    optSolData->mGcmmaData.AsympAdapt0 = 0.5;
    optSolData->mGcmmaData.AsympAdapt  = 0.7;
    optSolData->mGcmmaData.AsympAdaptC = 1./0.7;
    optSolData->mGcmmaData.StepSize    = 0.05;
    optSolData->mGcmmaData.Penalty     = 100.0;

    optSolData->mDumpAdvs = ADV_LAST;
    optSolData->mRestartIteration = 0; 
    aOptimizer->mWriteOptScalarValuesToFile = true;
    
    optSolData->mFDCheck = sFDCheck;    
    optSolData->mFDEps = sFDEps;
    optSolData->mFDSens = sFDSens; 

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
    UInt indices = sNumIts; // Last time step

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_STRAIN_ENERGY,true,true);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_MASS_UTILIZATION,false,true);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_UX_SQ,true,false);
    CritEvalFunc = new FunctionScalarProduct(entries,&weights);
    crit->CreateEvaluationFunction(CritEvalFunc,&indices);

    // First bool is dependent on states. Second bool dependent on Adv's
    crit = aOptimizer->CreateCriteria(GSOL_UY_SQ,true,false);
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
// ----------------------------------------------
// Set number of abstract design variables
// ----------------------------------------------
// Do some craziness to get the number of abstract design variables
    Mesh* mesh   = aProblem->GetMesh(0);
    UInt modelId = 1;
    Model* model = aProblem->GetModelById(modelId);

    // Get block 1 info
    ParseModelHelper parseModelHelper = ParseModelHelper();
    UInt *ElemIndices  = NULL, *ElemIds = NULL;
    UInt  numElemInSet = parseModelHelper.GetBaseSetElemsInfoByID(mesh,sBSetIds[0],ElemIndices,ElemIds);
    ParseModelHelper parseModelHelper1 = ParseModelHelper();
    UInt *NodeIndices  = NULL, *NodeIds = NULL;
    UInt  numNodeInSet = parseModelHelper1.GetBaseSetMasterNodesInfoByID(mesh,sBSetIds[0],NodeIndices,NodeIds), totNumNodesInSet=numNodeInSet;

    // Set the number of design variables per processor
#ifdef PARALLEL
    MPI_Allreduce(&numNodeInSet,&totNumNodesInSet,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
#endif
///----------------------------------------------
    aOptimizer->SetNumDesVars(totNumNodesInSet);
///----------------------------------------------




    // Create physical variables for the element's DENSITY in block 1
    UInt propId                = 30;
    DesignDepObject desDepObj  = MAT_PROP;
    int propTypeEmod           = MP_STRUC_EX;
    int propTypeRho            = MP_STRUC_RHO;
    PhysicalDesignVar* pdv     = NULL;
    ///////////////////////////////////////////////////////////////////////////
    // Assign Pdv's
    ///////////////////////////////////////////////////////////////////////////
    std::vector<UInt> NodeIndicesVec(FP_MAX_NODES,0), NodeIdsVec(FP_MAX_NODES,0);
    Function* EmodPDVEvalFunc=NULL;
    Function* RhoPDVEvalFunc=NULL;

    for (UInt ie=0, paramIndex=0, numNodes=0, inod=0, eleInd=0; ie < numElemInSet; ++ie)
    {
        eleInd      = ElemIndices[ie];
        paramIndex  = 0;
        pdv         = aOptimizer->CreatePhysDesVar(aProblem,modelId,desDepObj,propId,propTypeEmod,eleInd,paramIndex);

        /// Create the evaluation function for this physical design variable
        numNodes    = model->GetElement(eleInd)->GetNodeIds(&NodeIdsVec[0]);

        Real* EModFuncParamsVec = (Real*) alloca(sizeof(Real)*numNodes);
        Real* RhoFuncParamsVec = (Real*) alloca(sizeof(Real)*numNodes);
        UInt EModNumParams[] = {numNodes,1};
        UInt RhoNumParams[] = {numNodes,0};

        for (inod=0; inod<numNodes; ++inod)
        {
            NodeIndicesVec[inod] = NodeIdsVec[inod] - 1;
            EModFuncParamsVec[inod] = 1.0/numNodes;
            RhoFuncParamsVec[inod] = 1.0/numNodes;
        }
        EModFuncParamsVec[numNodes] =		sSimpExp;

        EmodPDVEvalFunc = new FunctionScalarProduct(EModNumParams,EModFuncParamsVec);
        pdv->CreateEvaluationFunction(EmodPDVEvalFunc,&NodeIndicesVec[0]);

        pdv      = aOptimizer->CreatePhysDesVar(aProblem,modelId,desDepObj,propId,propTypeRho,eleInd,paramIndex);
        RhoPDVEvalFunc = new FunctionScalarProduct(RhoNumParams,RhoFuncParamsVec);
        pdv->CreateEvaluationFunction(RhoPDVEvalFunc,&NodeIndicesVec[0]);
    }

    ///////////////////////////////////////////////////////////////////////////
    // For writing the adjoints
    ///////////////////////////////////////////////////////////////////////////
//    const enum SolVarType solNodalVars[] = {SOL_UX,SOL_UY,SOL_UZ,SOL_STRUC_SXX,SOL_STRUC_SYY,SOL_STRUC_SZZ,SOL_STRUC_SXY,SOL_STRUC_SXZ,SOL_STRUC_SYZ,SOL_STRUC_SVM,SOL_BATT_CONC1}; // End with Battery Conc1 because c++ sucks
//    int BlockIdsToPostProc[] = {1,-1}; // End with -1 because c++ sucks
//    const bool dumpNodalSensFlag= true;
//    const enum CriteriaTypes critTypes[] = {STRAIN_ENERGY_FINAL_STEP,DENSITY_VOLUME,XDISPLACEMENT_AVERAGE,YDISPLACEMENT_AVERAGE,ZDISPLACEMENT_AVERAGE,OTHER_CRIT}; // End with OTHER_CRIT because c++ sucks
//
//    WriteNodalAdjointsAndSensitivities(aProblem,solNodalVars,critTypes,BlockIdsToPostProc,sNumDofsPerNodeMax,dumpNodalSensFlag);
}



///----------------------------------------------------------------------------------------
/// Initialize the Vector of abstract design variables (Adv) mAbsDesVarVec
///----------------------------------------------------------------------------------------
extern "C"
    void InitializeAbsDesVarVector(Optimizer* aOptimizer, UInt aNumMyAbsDesVars)
{
    // Get abstract design variables vector
    Vector* absDesVarVec = aOptimizer->GetAbsDesVarVec();

    // get mesh
    Problem* problem = aOptimizer->GetProblem();
    Mesh*    mesh = problem->GetMesh(0);

    // Get the model
    Model* model = problem->GetModelById(1);

    // set the bounds for mAbsDesVarVec
    const Real AdvLowerBound = 1.0e-3;
    const Real AdvUpperBound =  1.0;

    // Read designvariables from file
    bool restartFromFile=false;
    if (restartFromFile)
    {
        // Get all nodes information
        ParseModelHelper parseModelHelper = ParseModelHelper();
        UInt *NodeIndices  = NULL, *NodeIds = NULL;
        UInt  numAllNodes = parseModelHelper.GetAllNodesInfo(mesh,NodeIndices,NodeIds);
//        UInt  numAllNodes = parseModelHelper.GetAllGlobalNodesInfo(model,mesh,NodeIndices,NodeIds);
        UInt  totNumNodes = model->GetGlbNumNodes();

        // allocate memory space for storing variables read from file
        std::vector<double> rstrtAbsDesVarsVec(totNumNodes,0.0);

        // read otimization variables from restart file
        char buffer[100];
        int iteration = 199;
        std::sprintf(buffer,"./Results_20x10/VolumePct_30/Pll_BICGSTAB/Run01_SIMP3/AbsDesVariables%04i.dat",iteration);
        FILE* restartFile = std::fopen(buffer,"rb");
        if (restartFile)
        {
            std::fprintf(stdout,"\nJust before read\n");
            std::fread(&rstrtAbsDesVarsVec[0],sizeof(double),totNumNodes,restartFile);
            std::fprintf(stdout,"\nJust after read\n");
            fclose(restartFile);
            for (UInt inod=0, nodInd=0; inod<numAllNodes; ++inod)
            {
                nodInd = NodeIds[inod]-1;
                absDesVarVec->ReplaceGlobalValues(1,&nodInd,&rstrtAbsDesVarsVec[nodInd]);
            }
        }

        absDesVarVec->GlobalAssemble();
    }
    else
    {
        absDesVarVec->PutScalar(sMassRatio);
    }
    aOptimizer->FillAbsDesVarVecBounds(AdvLowerBound,AdvUpperBound);

    /// Find and fix design variables at load points
    // Get nodeset nodes information
    ParseModelHelper parseModelHelper = ParseModelHelper();
    UInt *ElemIndices  = NULL, *ElemIds = NULL;
    UInt  numLoadElems = parseModelHelper.GetSideSetNeighborElemsInfoByID(mesh, sSSetIds[sLoadNSSetIndx],ElemIndices,ElemIds);
    std::vector<UInt> NodeIdsVec(FP_MAX_NODES,0);
    Real fixVal=AdvUpperBound, fixValLo=AdvUpperBound-1.0e-5;
    for (UInt ie=0, inod=0, nodInd=0, numNodes=0; ie<numLoadElems; ++ie)
    {
        numNodes    = model->GetElement(ElemIndices[ie])->GetNodeIds(&NodeIdsVec[0]);
        for (inod=0; inod<numNodes; ++inod)
        {
            nodInd = NodeIdsVec[inod]-1;
            absDesVarVec->ReplaceGlobalValues(1,&nodInd,&fixVal);
            aOptimizer->SetUpperAndLowerBounds(1,&nodInd, &fixValLo,1,&nodInd,&fixVal);
        }
    }
}





