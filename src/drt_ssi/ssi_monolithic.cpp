/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"

#include "ssi_coupling.H"
#include "ssi_manifold_flux_evaluator.H"
#include "ssi_monolithic_assemble_strategy.H"
#include "ssi_monolithic_contact_strategy.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_monolithic_dbc_handler.H"
#include "ssi_monolithic_evaluate_OffDiag.H"
#include "ssi_monolithic_meshtying_strategy.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_contact/contact_nitsche_strategy_ssi.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_elch.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_Time.h>


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SSIMono::SSIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams),
      contact_strategy_nitsche_(Teuchos::null),
      dbc_handler_(Teuchos::null),
      dtele_(0.0),
      dtsolve_(0.0),
      equilibration_method_{Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE")},
      manifoldscatraflux_(Teuchos::null),
      map_structure_(Teuchos::null),
      maps_scatra_(Teuchos::null),
      maps_sub_problems_(Teuchos::null),
      maps_systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      scatrastructureOffDiagcoupling_(Teuchos::null),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      ssi_maps_(Teuchos::null),
      ssi_matrices_(Teuchos::null),
      ssi_vectors_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_contact_(Teuchos::null),
      strategy_convcheck_(Teuchos::null),
      strategy_equilibration_(Teuchos::null),
      strategy_meshtying_(Teuchos::null),
      timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyContactToSubProblems()
{
  // uncomplete matrices; we need to do this here since in contact simulations the dofs that
  // interact with each other can change and thus the graph of the matrix can also change.
  ssi_matrices_->ScaTraMatrix()->UnComplete();
  ssi_matrices_->ScaTraStructureMatrix()->UnComplete();
  ssi_matrices_->StructureScaTraMatrix()->UnComplete();

  // add contributions
  strategy_contact_->ApplyContactToScatraResidual(ssi_vectors_->ScatraResidual());
  strategy_contact_->ApplyContactToScatraScatra(ssi_matrices_->ScaTraMatrix());
  strategy_contact_->ApplyContactToScatraStructure(ssi_matrices_->ScaTraStructureMatrix());
  strategy_contact_->ApplyContactToStructureScatra(ssi_matrices_->StructureScaTraMatrix());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyDBCToSystem()
{
  // apply Dirichlet boundary conditions to global system matrix
  dbc_handler_->ApplyDBCToSystemMatrix(ssi_matrices_->SystemMatrix());

  // apply Dirichlet boundary conditions to global RHS
  dbc_handler_->ApplyDBCToRHS(ssi_vectors_->Residual());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool SSI::SSIMono::IsUncompleteOfMatricesNecessaryForMeshTying() const
{
  // check for first iteration in calculation of initial time derivative
  if (IterationCount() == 0 and Step() == 0 and !DoCalculateInitialPotentialField()) return true;

  if (IterationCount() == 1)
  {
    // check for first iteration in calculation of initial potential field
    if (Step() == 0 and DoCalculateInitialPotentialField()) return true;

    // check for first iteration in restart simulations
    if (IsRestart())
    {
      auto* problem = DRT::Problem::Instance();
      // restart based on time step
      if (Step() == problem->Restart() + 1) return true;
      // restart based on time
      if (Time() == problem->RestartTime() + Dt()) return true;
    }
  }

  return false;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyMeshtyingToSubProblems()
{
  if (SSIInterfaceMeshtying())
  {
    // check if matrices are filled because they have to be for the below methods
    if (!ssi_matrices_->StructureScaTraMatrix()->Filled())
      ssi_matrices_->CompleteStructureScaTraMatrix();
    if (!ssi_matrices_->ScaTraStructureMatrix()->Filled())
      ssi_matrices_->CompleteScaTraStructureMatrix();

    if (IsScaTraManifold())
    {
      if (!ssi_matrices_->ScaTraManifoldStructureMatrix()->Filled())
        ssi_matrices_->CompleteScaTraManifoldStructureMatrix();

      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          ssi_matrices_->ScaTraManifoldStructureMatrix(),
          IsUncompleteOfMatricesNecessaryForMeshTying());

      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          manifoldscatraflux_->MatrixManifoldStructure(),
          IsUncompleteOfMatricesNecessaryForMeshTying());

      strategy_meshtying_->ApplyMeshtyingToScatraStructure(
          manifoldscatraflux_->MatrixScaTraStructure(), true);
    }

    strategy_meshtying_->ApplyMeshtyingToScatraStructure(
        ssi_matrices_->ScaTraStructureMatrix(), IsUncompleteOfMatricesNecessaryForMeshTying());

    strategy_meshtying_->ApplyMeshtyingToStructureMatrix(
        *ssi_matrices_->StructureMatrix(), StructureField()->SystemMatrix());

    strategy_meshtying_->ApplyMeshtyingToStructureScatra(
        ssi_matrices_->StructureScaTraMatrix(), IsUncompleteOfMatricesNecessaryForMeshTying());

    ssi_vectors_->StructureResidual()->Update(
        1.0, strategy_meshtying_->ApplyMeshtyingToStructureRHS(StructureField()->RHS()), 1.0);
  }
  // copy the structure residual and matrix if we do not have a mesh tying problem
  else
  {
    ssi_vectors_->StructureResidual()->Update(1.0, *(StructureField()->RHS()), 1.0);
    ssi_matrices_->StructureMatrix()->Add(*StructureField()->SystemMatrix(), false, 1.0, 1.0);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatAndRHS()
{
  AssembleMatScaTra();

  AssembleMatStructure();

  if (IsScaTraManifold()) AssembleMatScaTraManifold();

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(ssi_vectors_->Residual(), ssi_vectors_->ScatraResidual(),
      ssi_vectors_->StructureResidual(),
      IsScaTraManifold() ? ScaTraManifold()->Residual() : Teuchos::null,
      IsScaTraManifold() ? manifoldscatraflux_->RHSManifold() : Teuchos::null,
      IsScaTraManifold() ? manifoldscatraflux_->RHSScaTra() : Teuchos::null);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatScaTra()
{
  // assemble scatra-scatra block into system matrix
  strategy_assemble_->AssembleScatraScatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraMatrix());

  // assemble scatra-structure block into system matrix
  strategy_assemble_->AssembleScatraStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureMatrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatScaTraManifold()
{
  // assemble scatra manifold - scatra manifold block into system matrix
  strategy_assemble_->AssembleScatramanifoldScatramanifold(
      ssi_matrices_->SystemMatrix(), ScaTraManifold()->SystemMatrixOperator());

  // assemble scatra manifold-structure block into system matrix
  strategy_assemble_->AssembleScatramanifoldStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraManifoldStructureMatrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // manifold side
  strategy_assemble_->AssembleScatramanifoldScatramanifold(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixManifold());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixScaTra());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatramanifold(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixScaTraManifold());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // manifold side
  strategy_assemble_->AssembleScatramanifoldScatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixManifoldScatra());

  strategy_assemble_->AssembleScatramanifoldStructure(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixManifoldStructure());

  strategy_assemble_->AssembleScatraStructure(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixScaTraStructure());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatStructure()
{  // assemble structure-scatra block into system matrix
  strategy_assemble_->AssembleStructureScatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureScaTraMatrix());

  // assemble structure-structure block into system matrix
  strategy_assemble_->AssembleStructureStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureMatrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateSubproblems()
{
  // clear all matrices and residuals from previous Newton iteration
  ssi_matrices_->ClearMatrices();
  ssi_vectors_->ClearResiduals();

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // build system matrix and residual for structure field
  StructureField()->Evaluate();

  // build system matrix and residual for scalar transport field
  EvaluateScaTra();

  // build system matrix and residual for scalar transport field on manifold
  if (IsScaTraManifold()) EvaluateScaTraManifold();

  // build all off diagonal matrices
  EvaluateOffDiagContributions();

  // apply mesh tying to sub problems
  ApplyMeshtyingToSubProblems();

  // apply contact contributions to sub problems
  if (SSIInterfaceContact()) ApplyContactToSubProblems();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateOffDiagContributions()
{
  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      ssi_matrices_->ScaTraStructureMatrix());

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (SSIInterfaceMeshtying())
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssi_matrices_->ScaTraStructureMatrix());

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      ssi_matrices_->StructureScaTraMatrix());

  if (IsScaTraManifold())
  {
    // evaluate off-diagonal manifold-structure block of global system matrix
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraManifoldStructureDomain(
        ssi_matrices_->ScaTraManifoldStructureMatrix());
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::BuildNullSpaces() const
{
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      // equip smoother for scatra matrix blocks with null space
      ScaTraField()->BuildBlockNullSpaces(
          solver_, GetBlockPositions(Subproblem::scalar_transport)->at(0));
      if (IsScaTraManifold())
      {
        ScaTraManifold()->BuildBlockNullSpaces(
            solver_, GetBlockPositions(Subproblem::manifold)->at(0));
      }
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      std::ostringstream scatrablockstr;
      scatrablockstr << GetBlockPositions(Subproblem::scalar_transport)->at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_->Params().sublist("Inverse" + scatrablockstr.str());
      blocksmootherparamsscatra.sublist("Aztec Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsscatra);

      if (IsScaTraManifold())
      {
        std::ostringstream scatramanifoldblockstr;
        scatramanifoldblockstr << GetBlockPositions(Subproblem::manifold)->at(0) + 1;
        Teuchos::ParameterList& blocksmootherparamsscatramanifold =
            solver_->Params().sublist("Inverse" + scatramanifoldblockstr.str());
        blocksmootherparamsscatramanifold.sublist("Aztec Parameters");
        blocksmootherparamsscatramanifold.sublist("MueLu Parameters");

        // equip smoother for scatra matrix block with null space associated with all degrees of
        // freedom on scatra discretization
        ScaTraManifold()->Discretization()->ComputeNullSpaceIfNecessary(
            blocksmootherparamsscatramanifold);
      }

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with structural field as string
  std::stringstream iblockstr;
  iblockstr << GetBlockPositions(Subproblem::structure)->at(0) + 1;

  // equip smoother for structural matrix block with empty parameter sub lists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
}  // SSI::SSIMono::BuildNullSpaces

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::CompleteSubproblemMatrices()
{
  ssi_matrices_->ScaTraMatrix()->Complete();
  ssi_matrices_->CompleteScaTraStructureMatrix();
  ssi_matrices_->CompleteStructureScaTraMatrix();
  ssi_matrices_->StructureMatrix()->Complete();

  if (IsScaTraManifold())
  {
    ssi_matrices_->CompleteScaTraManifoldStructureMatrix();
    manifoldscatraflux_->CompleteMatrixManifoldStructure();
    manifoldscatraflux_->CompleteMatrixScaTraStructure();
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSIMono::DofRowMap() const
{
  return MapsSubProblems()->FullMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::SetupContactStrategy()
{
  // get the contact solution strategy
  auto contact_solution_type = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
      DRT::Problem::Instance()->ContactDynamicParams(), "STRATEGY");

  if (contact_solution_type == INPAR::CONTACT::solution_nitsche)
  {
    if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
            DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") !=
        INPAR::STR::int_standard)
      dserror("ssi contact only with new structural time integration");

    // get the contact model evaluator and store a pointer to the strategy
    auto& model_evaluator_contact = dynamic_cast<STR::MODELEVALUATOR::Contact&>(
        StructureField()->ModelEvaluator(INPAR::STR::model_contact));
    contact_strategy_nitsche_ = Teuchos::rcp_dynamic_cast<CONTACT::CoNitscheStrategySsi>(
        model_evaluator_contact.StrategyPtr(), true);
  }
  else
    dserror("Only Nitsche contact implemented for SSI problems at the moment!");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // check input parameters for scalar transport field
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

  // initialize strategy for Newton-Raphson convergence check
  switch (
      Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(globaltimeparams, "SCATRATIMINTTYPE"))
  {
    case INPAR::SSI::ScaTraTimIntType::elch:
    {
      if (IsScaTraManifold())
      {
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyElchScaTraManifold(globaltimeparams));
      }
      else
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyElch(globaltimeparams));
      break;
    }

    case INPAR::SSI::ScaTraTimIntType::standard:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyStd(globaltimeparams));
      break;
    }

    default:
    {
      dserror("Type of scalar transport time integrator currently not supported!");
      break;
    }
  }

  // call base class routine
  SSIBase::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Output()
{
  // output scalar transport field
  ScaTraField()->Output();
  if (IsScaTraManifold())
  {
    // domain output
    ScaTraManifold()->Output();
    // coupling output
    if (manifoldscatraflux_->DoOutput()) manifoldscatraflux_->Output();
  }

  // output structure field
  StructureField()->Output();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ReadRestart(int restart)
{
  // call base class
  SSIBase::ReadRestart(restart);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    SetupContactStrategy();
    SetSSIContactStates(ScaTraField()->Phinp());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ReadRestartfromTime(double restarttime)
{
  // call base class
  SSIBase::ReadRestartfromTime(restarttime);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    SetupContactStrategy();
    SetSSIContactStates(ScaTraField()->Phinp());
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareTimeLoop()
{
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  ScaTraField()->Output();
  if (IsScaTraManifold()) ScaTraManifold()->Output();

  // calculate initial potential field if needed
  if (DoCalculateInitialPotentialField()) CalcInitialPotentialField();

  // calculate initial time derivatives
  CalcInitialTimeDerivative();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // prepare time step for scalar transport field
  ScaTraField()->PrepareTimeStep();
  if (IsScaTraManifold()) ScaTraManifold()->PrepareTimeStep();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (PrepareTimeStep() of Scatra) and pass to other fields
  if (ScaTraField()->TimeStepAdapted()) SetDtFromScaTraToSSI();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(ScaTraField()->Phinp());
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // prepare time step for structural field
  StructureField()->PrepareTimeStep();

  // print time step information to screen
  ScaTraField()->PrintTimeStepInfo();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Setup()
{
  // call base class routine
  SSIBase::Setup();

  // safety checks
  if (ScaTraField()->NumScal() != 1)
  {
    dserror(
        "Since the ssi_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only "
        "one transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }
  const auto ssi_params = DRT::Problem::Instance()->SSIControlParams();

  const bool calc_initial_pot_elch =
      DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->ELCHControlParams(), "INITPOTCALC");
  const bool calc_initial_pot_ssi =
      DRT::INPUT::IntegralValue<bool>(ssi_params.sublist("ELCH"), "INITPOTCALC");

  if (ScaTraField()->EquilibrationMethod() != LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");
  }
  if (equilibration_method_.global != LINALG::EquilibrationMethod::local and
      (equilibration_method_.structure != LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != LINALG::EquilibrationMethod::none))
    dserror("Either global equilibration or local equilibration");

  if (matrixtype_ == LINALG::MatrixType::sparse and
      (equilibration_method_.structure != LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != LINALG::EquilibrationMethod::none))
    dserror("Block based equilibration only for block matrices");

  if (!DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ScalarTransportDynamicParams(), "SKIPINITDER"))
  {
    dserror(
        "Initial derivatives are already calculated in monolithic SSI. Enable 'SKIPINITDER' in the "
        "input file.");
  }

  if (calc_initial_pot_elch)
    dserror("Initial potential is calculated by SSI. Disable in Elch section.");
  if (calc_initial_pot_ssi and Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(ssi_params,
                                   "SCATRATIMINTTYPE") != INPAR::SSI::ScaTraTimIntType::elch)
    dserror("Calculation of initial potential only in case of Elch");

  if (!ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for monolithic scalar-structure interaction!");

  if (SSIInterfaceMeshtying() and
      MeshtyingStrategyS2I()->CouplingType() != INPAR::S2I::coupling_matching_nodes)
  {
    dserror(
        "Monolithic scalar-structure interaction only implemented for scatra-scatra "
        "interface coupling with matching interface nodes!");
  }

  if (SSIInterfaceContact() and !IsRestart()) SetupContactStrategy();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::SetupSystem()
{
  // setup the ssi maps object
  ssi_maps_ = Teuchos::rcp(new SSI::UTILS::SSIMaps(*this));

  // merge slave and master side block maps for interface matrix for scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);

  if (SSIInterfaceMeshtying())
  {
    // check whether slave-side degrees of freedom are Dirichlet-free
    std::vector<Teuchos::RCP<const Epetra_Map>> maps(2, Teuchos::null);
    maps[0] = InterfaceCouplingAdapterStructure()->SlaveDofMap();
    maps[1] = StructureField()->GetDBCMapExtractor()->CondMap();
    if (LINALG::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      dserror("Must not apply Dirichlet conditions to slave-side structural displacements!");

    interface_map_scatra = LINALG::MultiMapExtractor::MergeMaps(
        {MeshtyingStrategyS2I()->CouplingAdapter()->MasterDofMap(),
            MeshtyingStrategyS2I()->CouplingAdapter()->SlaveDofMap()});
  }

  // initialize global map extractor
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(
      IsScaTraManifold() ? 3 : 2, Teuchos::null);
  Teuchos::RCP<const Epetra_Map> merged_map;

  partial_maps[UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)] =
      Teuchos::rcp(new Epetra_Map(*ScaTraField()->DofRowMap()));
  partial_maps[UTILS::SSIMaps::GetProblemPosition(Subproblem::structure)] =
      Teuchos::rcp(new Epetra_Map(*StructureField()->DofRowMap()));
  if (IsScaTraManifold())
  {
    partial_maps[UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)] =
        Teuchos::rcp(new Epetra_Map(*ScaTraManifold()->DofRowMap()));
    auto temp_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
    merged_map = LINALG::MergeMap(temp_map, partial_maps[2], false);
  }
  else
    merged_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);

  maps_sub_problems_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_sub_problems_->CheckForValidMapExtractor();

  // initialize map extractors associated with blocks of global system matrix
  switch (ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      maps_systemmatrix_ = MapsSubProblems();
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix;

      // store an RCP to the block maps of the scatra field
      maps_scatra_ = Teuchos::rcpFromRef(ScaTraField()->BlockMaps());
      maps_scatra_->CheckForValidMapExtractor();

      if (IsScaTraManifold())
      {
        auto maps_scatra_manifold = Teuchos::rcpFromRef(ScaTraManifold()->BlockMaps());
        maps_scatra_manifold->CheckForValidMapExtractor();
        maps_systemmatrix.resize(GetBlockPositions(Subproblem::scalar_transport)->size() +
                                 GetBlockPositions(Subproblem::structure)->size() +
                                 GetBlockPositions(Subproblem::manifold)->size());

        for (int imap = 0; imap < static_cast<int>(GetBlockPositions(Subproblem::manifold)->size());
             ++imap)
        {
          maps_systemmatrix[GetBlockPositions(Subproblem::manifold)->at(imap)] =
              maps_scatra_manifold->Map(imap);
        }
      }
      else
      {
        // extract maps underlying main-diagonal matrix blocks associated with scalar transport
        // field
        maps_systemmatrix.resize(GetBlockPositions(Subproblem::scalar_transport)->size() +
                                 GetBlockPositions(Subproblem::structure)->size());
      }

      for (int imap = 0;
           imap < static_cast<int>(GetBlockPositions(Subproblem::scalar_transport)->size()); ++imap)
      {
        maps_systemmatrix[GetBlockPositions(Subproblem::scalar_transport)->at(imap)] =
            maps_scatra_->Map(imap);
      }

      // extract map underlying single main-diagonal matrix block associated with structural
      // field
      maps_systemmatrix[GetBlockPositions(Subproblem::structure)->at(0)] =
          StructureField()->DofRowMap();

      // initialize map extractor associated with blocks of global system matrix
      maps_systemmatrix_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*DofRowMap(), maps_systemmatrix));

      // initialize map extractor associated with all degrees of freedom inside structural field
      map_structure_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*StructureField()->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, StructureField()->DofRowMap())));

      // safety check
      map_structure_->CheckForValidMapExtractor();

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  maps_systemmatrix_->CheckForValidMapExtractor();

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_field:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // feed AMGnxn block preconditioner with null space information for each block of global
      // block system matrix
      BuildNullSpaces();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // safety check
      if (ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize sub blocks and system matrix
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(*this));

  // initialize residual and increment vectors
  ssi_vectors_ = Teuchos::rcp(new SSI::UTILS::SSIVectors(*this));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(*this, matrixtype_, ScaTraField()->MatrixType());

  if (IsScaTraManifold())
  {
    // initialize object, that performs evaluations of OD coupling
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraManifoldStructureOffDiagCoupling(
        MapStructure(),
        MapsSubProblems()->Map(UTILS::SSIMaps::GetProblemPosition(Subproblem::structure)),
        InterfaceCouplingAdapterStructure(), InterfaceCouplingAdapterStructure3DomainIntersection(),
        interface_map_scatra, MeshtyingStrategyS2I(), ScaTraBaseAlgorithm(),
        ScaTraManifoldBaseAlgorithm(), StructureField(), Meshtying3DomainIntersection()));

    // initialize object, that performs evaluations of scatra - scatra on manifold coupling
    manifoldscatraflux_ = Teuchos::rcp(new SSI::ScaTraManifoldScaTraFluxEvaluator(*this));
  }
  else
  {
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
        MapStructure(),
        MapsSubProblems()->Map(UTILS::SSIMaps::GetProblemPosition(Subproblem::structure)),
        InterfaceCouplingAdapterStructure(), InterfaceCouplingAdapterStructure3DomainIntersection(),
        interface_map_scatra, MeshtyingStrategyS2I(), ScaTraBaseAlgorithm(), StructureField(),
        Meshtying3DomainIntersection()));
  }
  // instantiate appropriate equilibration class
  strategy_equilibration_ = LINALG::BuildEquilibration(
      matrixtype_, GetBlockEquilibration(), MapsSubProblems()->FullMap());

  // instantiate appropriate contact class
  strategy_contact_ = SSI::BuildContactStrategy(*this, ScaTraField()->MatrixType());

  // instantiate appropriate mesh tying class
  strategy_meshtying_ =
      SSI::BuildMeshtyingStrategy(*this, matrixtype_, ScaTraField()->MatrixType());

  // instantiate Dirichlet boundary condition handler class
  dbc_handler_ = SSI::BuildDBCHandler(Teuchos::rcp(this, false), matrixtype_);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetupModelEvaluator() const
{
  // construct and register structural model evaluator if necessary

  const bool do_output_stress =
      DRT::INPUT::IntegralValue<INPAR::STR::StressType>(
          DRT::Problem::Instance()->IOParams(), "STRUCT_STRESS") != INPAR::STR::stress_none;
  const bool smooth_output_interface_stress = DRT::INPUT::IntegralValue<bool>(
      DRT::Problem::Instance()->SSIControlParams().sublist("MONOLITHIC"),
      "SMOOTH_OUTPUT_INTERFACE_STRESS");

  if (Meshtying3DomainIntersection() and smooth_output_interface_stress)
    dserror("Smoothing of interface stresses not implemented for triple meshtying.");

  if (smooth_output_interface_stress and !do_output_stress)
    dserror("Smoothing of interface stresses only when stress output is written.");

  if (do_output_stress and SSIInterfaceMeshtying())
  {
    StructureBaseAlgorithm()->RegisterModelEvaluator("Monolithic Coupling Model",
        Teuchos::rcp(new STR::MODELEVALUATOR::MonolithicSSI(
            Teuchos::rcp(this, false), smooth_output_interface_stress)));
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // call base class
  SSIBase::SetScatraSolution(phi);

  // set state for contact evaluation
  SetSSIContactStates(phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetSSIContactStates(Teuchos::RCP<const Epetra_Vector> phi) const
{
  if (contact_strategy_nitsche_ != Teuchos::null)
    contact_strategy_nitsche_->SetState(MORTAR::state_scalar, *phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SolveLinearSystem()
{
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), ssi_vectors_->Residual(), *MapsSystemMatrix());

  // solve global system of equations
  // Dirichlet boundary conditions have already been applied to global system of equations
  solver_->Solve(ssi_matrices_->SystemMatrix()->EpetraOperator(), ssi_vectors_->Increment(),
      ssi_vectors_->Residual(), true, IterationCount() == 1);

  strategy_equilibration_->UnequilibrateIncrement(ssi_vectors_->Increment());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::NewtonLoop()
{
  // reset counter for Newton-Raphson iteration
  ResetIterationCount();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    IncrementIterationCount();

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // set solution from last Newton step to all fields
    DistributeSolutionAllFields();

    // evaluate sub problems and get all matrices and right-hand-sides
    EvaluateSubproblems();

    // complete the sub problem matrices
    CompleteSubproblemMatrices();

    // assemble global system of equations
    AssembleMatAndRHS();

    // apply the Dirichlet boundary conditions to global system
    ApplyDBCToSystem();

    // determine time needed for evaluating elements and assembling global system of
    // equations, and take maximum over all processors via communication
    double mydtele = timer_->WallTime() - time;
    Comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // clear the global increment vector
    ssi_vectors_->ClearIncrement();

    // store time before solving global system of equations
    time = timer_->WallTime();

    SolveLinearSystem();

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime() - time;
    Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (DRT::INPUT::IntegralValue<bool>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->OutputLinSolverStats(*solver_, dtsolve_, Step(), IterationCount(),
          ssi_vectors_->Residual()->Map().NumGlobalElements());

    // update states for next Newton iteration
    UpdateIterScaTra();
    UpdateIterStructure();

  }  // Newton-Raphson iteration
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Timeloop()
{
  if (Step() == 0) PrepareTimeLoop();

  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    const double time = timer_->WallTime();

    // evaluate time step
    NewtonLoop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->WallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if
    // applicable
    if (DRT::INPUT::IntegralValue<int>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      ScaTraField()->OutputNonlinSolverStats(IterationCount(), dtnonlinsolve, Step(), Comm());

    PrepareOutput();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::Update()
{
  // update scalar transport field
  ScaTraField()->Update();
  if (IsScaTraManifold()) ScaTraManifold()->Update();

  // update structure field
  StructureField()->Update();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterScaTra()
{
  // update scalar transport field
  ScaTraField()->UpdateIter(MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)));
  ScaTraField()->ComputeIntermediateValues();

  if (IsScaTraManifold())
  {
    ScaTraManifold()->UpdateIter(MapsSubProblems()->ExtractVector(
        ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)));
    ScaTraManifold()->ComputeIntermediateValues();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterStructure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure = MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::structure));

  // consider structural meshtying. Copy master increments and displacements to slave side.
  if (SSIInterfaceMeshtying())
  {
    MapsCoupStruct()->InsertVector(
        InterfaceCouplingAdapterStructure()->MasterToSlave(
            MapsCoupStruct()->ExtractVector(StructureField()->Dispnp(), 2)),
        1, StructureField()->WriteAccessDispnp());
    StructureField()->SetState(StructureField()->WriteAccessDispnp());
    MapsCoupStruct()->InsertVector(InterfaceCouplingAdapterStructure()->MasterToSlave(
                                       MapsCoupStruct()->ExtractVector(increment_structure, 2)),
        1, increment_structure);

    if (Meshtying3DomainIntersection())
    {
      MapsCoupStruct3DomainIntersection()->InsertVector(
          InterfaceCouplingAdapterStructure3DomainIntersection()->MasterToSlave(
              MapsCoupStruct3DomainIntersection()->ExtractVector(StructureField()->Dispnp(), 2)),
          1, StructureField()->WriteAccessDispnp());
      StructureField()->SetState(StructureField()->WriteAccessDispnp());
      MapsCoupStruct3DomainIntersection()->InsertVector(
          InterfaceCouplingAdapterStructure3DomainIntersection()->MasterToSlave(
              MapsCoupStruct3DomainIntersection()->ExtractVector(increment_structure, 2)),
          1, increment_structure);
    }
  }

  // update displacement of structure field
  StructureField()->UpdateStateIncrementally(increment_structure);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<int>> SSI::SSIMono::GetBlockPositions(Subproblem subproblem) const
{
  if (matrixtype_ == LINALG::MatrixType::sparse) dserror("Sparse matrices have just one block");

  Teuchos::RCP<std::vector<int>> block_position = Teuchos::rcp(new std::vector<int>(0));

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      if (ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(1);
      else
        block_position->emplace_back(ScaTraField()->BlockMaps().NumMaps());
      break;
    }
    case Subproblem::scalar_transport:
    {
      if (ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(0);
      else
      {
        for (int i = 0; i < ScaTraField()->BlockMaps().NumMaps(); ++i)
          block_position->emplace_back(i);
      }
      break;
    }
    case Subproblem::manifold:
    {
      if (ScaTraManifold()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(2);
      else
      {
        for (int i = 0; i < static_cast<int>(ScaTraManifold()->BlockMaps().NumMaps()); ++i)
          block_position->emplace_back(ScaTraField()->BlockMaps().NumMaps() + 1 + i);
      }
      break;
    }
    default:
    {
      dserror("Unknown type of subproblem");
      break;
    }
  }

  return block_position;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<LINALG::EquilibrationMethod>> SSI::SSIMono::GetBlockEquilibration()
{
  Teuchos::RCP<std::vector<LINALG::EquilibrationMethod>> equilibration_method_vector;
  switch (matrixtype_)
  {
    case LINALG::MatrixType::sparse:
    {
      equilibration_method_vector = Teuchos::rcp(
          new std::vector<LINALG::EquilibrationMethod>(1, equilibration_method_.global));
      break;
    }
    case LINALG::MatrixType::block_field:
    {
      if (equilibration_method_.global != LINALG::EquilibrationMethod::local)
      {
        equilibration_method_vector = Teuchos::rcp(
            new std::vector<LINALG::EquilibrationMethod>(1, equilibration_method_.global));
      }
      else if (equilibration_method_.structure == LINALG::EquilibrationMethod::none and
               equilibration_method_.scatra == LINALG::EquilibrationMethod::none)
      {
        equilibration_method_vector = Teuchos::rcp(
            new std::vector<LINALG::EquilibrationMethod>(1, LINALG::EquilibrationMethod::none));
      }
      else
      {
        auto block_positions_scatra = GetBlockPositions(Subproblem::scalar_transport);
        auto block_position_structure = GetBlockPositions(Subproblem::structure);
        auto block_positions_scatra_manifold =
            IsScaTraManifold() ? GetBlockPositions(Subproblem::manifold) : Teuchos::null;

        equilibration_method_vector = Teuchos::rcp(new std::vector<LINALG::EquilibrationMethod>(
            block_positions_scatra->size() + block_position_structure->size() +
                (IsScaTraManifold() ? block_positions_scatra_manifold->size() : 0),
            LINALG::EquilibrationMethod::none));

        for (const int block_position_scatra : *block_positions_scatra)
          equilibration_method_vector->at(block_position_scatra) = equilibration_method_.scatra;

        equilibration_method_vector->at(block_position_structure->at(0)) =
            equilibration_method_.structure;

        if (IsScaTraManifold())
        {
          for (const int block_position_scatra_manifold : *block_positions_scatra_manifold)
          {
            equilibration_method_vector->at(block_position_scatra_manifold) =
                equilibration_method_.scatra;
          }
        }
      }

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with system matrix field!");
      break;
    }
  }
  return equilibration_method_vector;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateScaTra()
{
  // evaluate the scatra field
  ScaTraField()->PrepareLinearSolve();

  // copy the matrix to the corresponding ssi matrix and complete it such that additional
  // contributions like contact contributions can be added before assembly
  ssi_matrices_->ScaTraMatrix()->Add(*ScaTraField()->SystemMatrixOperator(), false, 1.0, 1.0);

  // copy the residual to the corresponding ssi vector to enable application of contact
  // contributions before assembly
  ssi_vectors_->ScatraResidual()->Update(1.0, *ScaTraField()->Residual(), 1.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateScaTraManifold()
{
  // evaluate single problem
  ScaTraManifold()->PrepareLinearSolve();

  // evaluate coupling fluxes
  manifoldscatraflux_->Evaluate();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareOutput()
{
  StructureField()->PrepareOutput();

  // prepare output of coupling sctra manifold - scatra
  if (IsScaTraManifold() and manifoldscatraflux_->DoOutput())
    manifoldscatraflux_->EvaluateScaTraManifoldInflow();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::DistributeSolutionAllFields(const bool restore_velocity)
{
  // needed to communicate to NOX state
  if (restore_velocity)
  {
    auto vel_temp = *StructureField()->Velnp();
    StructureField()->SetState(StructureField()->WriteAccessDispnp());
    StructureField()->WriteAccessVelnp()->Update(1.0, vel_temp, 0.0);
  }
  else
    StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // distribute states to other fields
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  SetScatraSolution(ScaTraField()->Phinp());
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::CalcInitialPotentialField()
{
  const auto equpot = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(
      DRT::Problem::Instance()->ELCHControlParams(), "EQUPOT");
  if (equpot != INPAR::ELCH::equpot_divi and equpot != INPAR::ELCH::equpot_enc_pde and
      equpot != INPAR::ELCH::equpot_enc_pde_elim)
  {
    dserror(
        "Initial potential field cannot be computed for chosen closing equation for electric "
        "potential!");
  }

  // store initial velocity to restore them afterwards
  auto init_velocity = *StructureField()->Velnp();

  // cast scatra time integrators to elch to call elch specific methods
  auto scatra_elch = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraField());
  auto manifold_elch = IsScaTraManifold()
                           ? Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraManifold())
                           : Teuchos::null;
  if (scatra_elch == Teuchos::null or (IsScaTraManifold() and manifold_elch == Teuchos::null))
    dserror("Cast to Elch time integrator faild. Scatra is not an Elch problem");

  // prepare specific time integrators
  scatra_elch->PreCalcInitialPotentialField();
  if (IsScaTraManifold()) manifold_elch->PreCalcInitialPotentialField();

  auto scatra_elch_splitter = ScaTraField()->Splitter();
  auto manifold_elch_splitter = IsScaTraManifold() ? ScaTraManifold()->Splitter() : Teuchos::null;

  ResetIterationCount();

  while (true)
  {
    IncrementIterationCount();

    // prepare full SSI system
    DistributeSolutionAllFields(true);
    EvaluateSubproblems();

    // complete the sub problem matrices
    CompleteSubproblemMatrices();

    AssembleMatAndRHS();
    ApplyDBCToSystem();

    // apply artificial Dirichlet boundary conditions to system of equations (on concentration
    // dofs and on structure dofs)
    Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
    if (IsScaTraManifold())
    {
      auto conc_map =
          LINALG::MergeMap(scatra_elch_splitter->OtherMap(), manifold_elch_splitter->OtherMap());
      pseudo_dbc_map = LINALG::MergeMap(conc_map, StructureField()->DofRowMap());
    }
    else
    {
      pseudo_dbc_map =
          LINALG::MergeMap(scatra_elch_splitter->OtherMap(), StructureField()->DofRowMap());
    }

    auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

    auto rhs = ssi_vectors_->Residual();
    ApplyDirichlettoSystem(
        ssi_matrices_->SystemMatrix(), rhs, Teuchos::null, dbc_zeros, *pseudo_dbc_map);
    ssi_vectors_->Residual()->Update(1.0, *rhs, 0.0);

    if (strategy_convcheck_->ExitNewtonRaphsonInitPotCalc(*this)) break;

    // solve for potential increments
    ssi_vectors_->ClearIncrement();
    SolveLinearSystem();

    // update potential dofs in scatra and manifold fields
    ScaTraField()->UpdateIter(MapsSubProblems()->ExtractVector(ssi_vectors_->Increment(),
        UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)));
    if (IsScaTraManifold())
      ScaTraManifold()->UpdateIter(MapsSubProblems()->ExtractVector(
          ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)));

    // copy initial state vector
    ScaTraField()->Phin()->Update(1.0, *ScaTraField()->Phinp(), 0.0);
    if (IsScaTraManifold()) ScaTraManifold()->Phin()->Update(1.0, *ScaTraManifold()->Phinp(), 0.0);

    // update state vectors for intermediate time steps (only for generalized alpha)
    ScaTraField()->ComputeIntermediateValues();
    if (IsScaTraManifold()) ScaTraManifold()->ComputeIntermediateValues();
  }

  scatra_elch->PostCalcInitialPotentialField();
  if (IsScaTraManifold()) manifold_elch->PostCalcInitialPotentialField();

  StructureField()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::CalcInitialTimeDerivative()
{
  // store initial velocity to restore them afterwards
  auto init_velocity = *StructureField()->Velnp();

  const bool is_elch = IsElchScaTraTimIntType();

  // prepare specific time integrators
  ScaTraField()->PreCalcInitialTimeDerivative();
  if (IsScaTraManifold()) ScaTraManifold()->PreCalcInitialTimeDerivative();

  auto scatra_elch_splitter = is_elch ? ScaTraField()->Splitter() : Teuchos::null;
  auto manifold_elch_splitter =
      (is_elch and IsScaTraManifold()) ? ScaTraManifold()->Splitter() : Teuchos::null;

  // initial screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << "Calculating initial time derivative of state variables on discretization "
              << ScaTraField()->Discretization()->Name();
    if (IsScaTraManifold())
      std::cout << " and discretization " << ScaTraManifold()->Discretization()->Name();
    std::cout << std::endl;
  }

  // evaluate Dirichlet and Neumann boundary conditions
  ScaTraField()->ApplyBCToSystem();
  if (IsScaTraManifold()) ScaTraManifold()->ApplyBCToSystem();

  // clear history values (this is the first step)
  ScaTraField()->Hist()->PutScalar(0.0);
  if (IsScaTraManifold()) ScaTraManifold()->Hist()->PutScalar(0.0);

  // In a first step, we assemble the standard global system of equations (we need the residual)
  DistributeSolutionAllFields(true);
  EvaluateSubproblems();

  // complete the sub problem matrices
  CompleteSubproblemMatrices();

  AssembleMatAndRHS();
  ApplyDBCToSystem();

  // prepare mass matrices of sub problems and global system
  auto massmatrix_scatra =
      ScaTraField()->MatrixType() == LINALG::MatrixType::sparse
          ? Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(
                UTILS::SSIMatrices::SetupSparseMatrix(ScaTraField()->DofRowMap()))
          : Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(UTILS::SSIMatrices::SetupBlockMatrix(
                Teuchos::rcpFromRef(ScaTraField()->BlockMaps()),
                Teuchos::rcpFromRef(ScaTraField()->BlockMaps())));

  auto massmatrix_manifold =
      IsScaTraManifold()
          ? (ScaTraManifold()->MatrixType() == LINALG::MatrixType::sparse
                    ? Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(
                          UTILS::SSIMatrices::SetupSparseMatrix(ScaTraManifold()->DofRowMap()))
                    : Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(
                          UTILS::SSIMatrices::SetupBlockMatrix(
                              Teuchos::rcpFromRef(ScaTraManifold()->BlockMaps()),
                              Teuchos::rcpFromRef(ScaTraManifold()->BlockMaps()))))
          : Teuchos::null;

  auto massmatrix_system =
      MatrixType() == LINALG::MatrixType::sparse
          ? Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(
                UTILS::SSIMatrices::SetupSparseMatrix(DofRowMap()))
          : Teuchos::rcp_dynamic_cast<LINALG::SparseOperator>(
                UTILS::SSIMatrices::SetupBlockMatrix(MapsSystemMatrix(), MapsSystemMatrix()));

  // fill ones on main diag of structure block (not solved)
  auto ones_struct = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(), true));
  ones_struct->PutScalar(1.0);
  MatrixType() == LINALG::MatrixType::sparse
      ? LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            *LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_system), *ones_struct)
      : LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system)
                ->Matrix(GetBlockPositions(Subproblem::structure)->at(0),
                    GetBlockPositions(Subproblem::structure)->at(0)),
            *ones_struct);

  // extract residuals of scatra and manifold from global residual
  auto rhs_scatra = Teuchos::rcp(new Epetra_Vector(*ScaTraField()->DofRowMap(), true));
  auto rhs_manifold = IsScaTraManifold()
                          ? Teuchos::rcp(new Epetra_Vector(*ScaTraManifold()->DofRowMap(), true))
                          : Teuchos::null;

  rhs_scatra->Update(1.0,
      *MapsSubProblems()->ExtractVector(ssi_vectors_->Residual(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)),
      0.0);
  if (IsScaTraManifold())
  {
    rhs_manifold->Update(1.0,
        *MapsSubProblems()->ExtractVector(
            ssi_vectors_->Residual(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)),
        0.0);
  }

  // In a second step, we need to modify the assembled system of equations, since we want to solve
  // M phidt^0 = f^n - K\phi^n - C(u_n)\phi^n
  // In particular, we need to replace the global system matrix by a global mass matrix,
  // and we need to remove all transient contributions associated with time discretization from the
  // global residual vector.

  // Evaluate mass matrix and modify residual
  ScaTraField()->EvaluateInitialTimeDerivative(massmatrix_scatra, rhs_scatra);
  if (IsScaTraManifold())
    ScaTraManifold()->EvaluateInitialTimeDerivative(massmatrix_manifold, rhs_manifold);

  // assemble global mass matrix from
  switch (MatrixType())
  {
    case LINALG::MatrixType::block_field:
    {
      switch (ScaTraField()->MatrixType())
      {
        case LINALG::MatrixType::block_condition:
        case LINALG::MatrixType::block_condition_dof:
        {
          auto massmatrix_system_block =
              LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system);

          auto massmatrix_scatra_block =
              LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_scatra);

          auto positions_scatra = GetBlockPositions(Subproblem::scalar_transport);

          for (int i = 0; i < static_cast<int>(positions_scatra->size()); ++i)
          {
            const int position_scatra = positions_scatra->at(i);
            massmatrix_system_block->Matrix(position_scatra, position_scatra)
                .Add(massmatrix_scatra_block->Matrix(i, i), false, 1.0, 1.0);
          }
          if (IsScaTraManifold())
          {
            auto positions_manifold = GetBlockPositions(Subproblem::manifold);

            auto massmatrix_manifold_block =
                LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_manifold);

            for (int i = 0; i < static_cast<int>(positions_manifold->size()); ++i)
            {
              const int position_manifold = positions_manifold->at(i);
              massmatrix_system_block->Matrix(position_manifold, position_manifold)
                  .Add(massmatrix_manifold_block->Matrix(i, i), false, 1.0, 1.0);
            }
          }

          break;
        }

        case LINALG::MatrixType::sparse:
        {
          auto massmatrix_system_block =
              LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system);

          const int position_scatra = GetBlockPositions(Subproblem::scalar_transport)->at(0);

          massmatrix_system_block->Matrix(position_scatra, position_scatra)
              .Add(*LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_scatra), false, 1.0, 1.0);

          if (IsScaTraManifold())
          {
            const int position_manifold = GetBlockPositions(Subproblem::manifold)->at(0);

            massmatrix_system_block->Matrix(position_manifold, position_manifold)
                .Add(*LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false, 1.0,
                    1.0);
          }
          break;
        }

        default:
        {
          dserror("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }
      massmatrix_system->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      auto massmatrix_system_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_system);
      massmatrix_system_sparse->Add(
          *LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_scatra), false, 1.0, 1.0);

      if (IsScaTraManifold())
        massmatrix_system_sparse->Add(
            *LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false, 1.0, 1.0);

      massmatrix_system->Complete(*DofRowMap(), *DofRowMap());
      break;
    }
    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // reconstruct global residual from partial residuals
  auto rhs_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*DofRowMap(), true));
  MapsSubProblems()->InsertVector(
      rhs_scatra, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport), rhs_system);
  if (IsScaTraManifold())
    MapsSubProblems()->InsertVector(
        rhs_manifold, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold), rhs_system);

  // apply artificial Dirichlet boundary conditions to system of equations to non-transported
  // scalars and structure
  Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
  if (IsScaTraManifold() and is_elch)
  {
    auto conc_map =
        LINALG::MergeMap(scatra_elch_splitter->CondMap(), manifold_elch_splitter->CondMap());
    pseudo_dbc_map = LINALG::MergeMap(conc_map, StructureField()->DofRowMap());
  }
  else if (is_elch)
  {
    pseudo_dbc_map =
        LINALG::MergeMap(scatra_elch_splitter->CondMap(), StructureField()->DofRowMap());
  }
  else
    pseudo_dbc_map = Teuchos::rcp(new Epetra_Map(*StructureField()->DofRowMap()));

  auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

  // temporal derivative of transported scalars
  auto phidtnp_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*DofRowMap(), true));
  LINALG::ApplyDirichlettoSystem(
      massmatrix_system, phidtnp_system, rhs_system, dbc_zeros, *(pseudo_dbc_map));

  // solve global system of equations for initial time derivative of state variables
  solver_->Solve(massmatrix_system->EpetraOperator(), phidtnp_system, rhs_system, true, true);

  // copy solution to sub problmes
  auto phidtnp_scatra = MapsSubProblems()->ExtractVector(
      phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport));
  ScaTraField()->Phidtnp()->Update(1.0, *phidtnp_scatra, 0.0);
  ScaTraField()->Phidtn()->Update(1.0, *phidtnp_scatra, 0.0);

  if (IsScaTraManifold())
  {
    auto phidtnp_manifold = MapsSubProblems()->ExtractVector(
        phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold));
    ScaTraManifold()->Phidtnp()->Update(1.0, *phidtnp_manifold, 0.0);
    ScaTraManifold()->Phidtn()->Update(1.0, *phidtnp_manifold, 0.0);
  }

  // reset solver
  solver_->Reset();

  ScaTraField()->PostCalcInitialTimeDerivative();
  if (IsScaTraManifold()) ScaTraManifold()->PostCalcInitialTimeDerivative();

  StructureField()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}