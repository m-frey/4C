/*----------------------------------------------------------------------------*/
/*! \file

 \brief ALE field adapter

 \level 1

 */
/*----------------------------------------------------------------------------*/

#include "4C_adapter_ale.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_ale_fpsi.hpp"
#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_ale_fsi_msht.hpp"
#include "4C_adapter_ale_xffsi.hpp"
#include "4C_ale.hpp"
#include "4C_fem_condition_periodic.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleBaseAlgorithm::AleBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  setup_ale(prbdyn, actdis);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleBaseAlgorithm::setup_ale(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ALE::AleBaseAlgorithm::setup_ale");
  Teuchos::TimeMonitor monitor(*t);

  // what's the current problem type?
  const Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();

  // ---------------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // ---------------------------------------------------------------------------
  if (!actdis->filled()) actdis->fill_complete();

  // ---------------------------------------------------------------------------
  // connect degrees of freedom for coupled nodes
  // ---------------------------------------------------------------------------
  Core::Conditions::PeriodicBoundaryConditions pbc(actdis);
  pbc.update_dofs_for_periodic_boundary_conditions();

  // ---------------------------------------------------------------------------
  // context for output and restart
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->writer();

  // Output for these problems are not necessary because we write
  // restart data at each time step for visualization
  output->write_mesh(0, 0.0);

  // ---------------------------------------------------------------------------
  // set some pointers and variables
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> adyn =
      Teuchos::RCP(new Teuchos::ParameterList(Global::Problem::instance()->ale_dynamic_params()));

  // ---------------------------------------------------------------------------
  // create a linear solver
  // ---------------------------------------------------------------------------
  // get the linear solver number
  const int linsolvernumber = adyn->get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "No linear solver defined for ALE problems. Please set "
        "LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::RCP(
      new Core::LinAlg::Solver(Global::Problem::instance()->solver_params(linsolvernumber),
          actdis->get_comm(), Global::Problem::instance()->solver_params_callback(),
          Teuchos::getIntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY")));
  actdis->compute_null_space_if_necessary(solver->params());

  // ---------------------------------------------------------------------------
  // overwrite certain parameters when ALE is part of a multi-field problem
  // ---------------------------------------------------------------------------
  adyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  adyn->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  adyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  adyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  adyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  if (probtype == Core::ProblemType::fpsi)
  {
    // FPSI input parameters
    const Teuchos::ParameterList& fpsidyn = Global::Problem::instance()->fpsi_dynamic_params();
    auto coupling = Teuchos::getIntegralValue<FpsiCouplingType>(fpsidyn, "COUPALGO");
    if (coupling == partitioned)
    {
      FOUR_C_THROW("partitioned fpsi solution scheme has not been implemented yet.");
    }
  }

  // create the ALE time integrator
  auto aletype = Teuchos::getIntegralValue<Inpar::ALE::AleDynamic>(*adyn, "ALE_TYPE");
  Teuchos::RCP<ALE::Ale> ale = Teuchos::null;
  switch (aletype)
  {
    // catch all nonlinear cases
    case Inpar::ALE::solid:
    case Inpar::ALE::laplace_spatial:
    case Inpar::ALE::springs_spatial:
    {
      ale = Teuchos::RCP(new ALE::Ale(actdis, solver, adyn, output));

      break;
    }
    // catch and linear cases
    case Inpar::ALE::solid_linear:
    case Inpar::ALE::laplace_material:
    case Inpar::ALE::springs_material:
    {
      ale = Teuchos::RCP(new ALE::AleLinear(actdis, solver, adyn, output));

      break;
    }
    default:
    {
      FOUR_C_THROW("Decide, whether ALE_TYPE = '%s' is linear or nonlinear.",
          adyn->get<std::string>("ALE_TYPE").c_str());
      break;
    }
  }

  /* determine problem type and then wrap the ALE time integrator into a
   * problem-specific wrapper */
  switch (probtype)
  {
    case Core::ProblemType::ale:
    {
      ale_ = ale;
      break;
    }
    case Core::ProblemType::fsi:
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::biofilm_fsi:
    {
      const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
      auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit_saddlepoint)
      {
        ale_ = Teuchos::RCP(new Adapter::AleFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_sliding_monolithicfluidsplit or
               coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        ale_ = Teuchos::RCP(new Adapter::AleFsiMshtWrapper(ale));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit or
               coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
               coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nonox or
               coupling == fsi_iter_fluidfluid_monolithicfluidsplit_nonox)
      {
        ale_ = Teuchos::RCP(new Adapter::AleXFFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_stagg_AITKEN_rel_force or
               coupling == fsi_iter_stagg_AITKEN_rel_param or
               coupling == fsi_iter_stagg_CHEB_rel_param or coupling == fsi_iter_stagg_MFNK_FD or
               coupling == fsi_iter_stagg_MFNK_FSI or coupling == fsi_iter_stagg_MPE or
               coupling == fsi_iter_stagg_NLCG or coupling == fsi_iter_stagg_Newton_FD or
               coupling == fsi_iter_stagg_Newton_I or coupling == fsi_iter_stagg_RRE or
               coupling == fsi_iter_stagg_fixed_rel_param or
               coupling == fsi_iter_stagg_steep_desc or coupling == fsi_iter_stagg_steep_desc_force)
      {
        ale_ = Teuchos::RCP(new Adapter::AleFluidWrapper(ale));
      }
      else
      {
        FOUR_C_THROW(
            "No ALE adapter available yet for your chosen FSI coupling "
            "algorithm!");
      }
      break;
    }
    case Core::ProblemType::fsi_redmodels:
    {
      const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
      auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_sliding_monolithicfluidsplit or
          coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        ale_ = Teuchos::RCP(new Adapter::AleFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_stagg_AITKEN_rel_force or
               coupling == fsi_iter_stagg_AITKEN_rel_param or
               coupling == fsi_iter_stagg_CHEB_rel_param or coupling == fsi_iter_stagg_MFNK_FD or
               coupling == fsi_iter_stagg_MFNK_FSI or coupling == fsi_iter_stagg_MPE or
               coupling == fsi_iter_stagg_NLCG or coupling == fsi_iter_stagg_Newton_FD or
               coupling == fsi_iter_stagg_Newton_I or coupling == fsi_iter_stagg_RRE or
               coupling == fsi_iter_stagg_fixed_rel_param or
               coupling == fsi_iter_stagg_steep_desc or coupling == fsi_iter_stagg_steep_desc_force)
      {
        ale_ = Teuchos::RCP(new Adapter::AleFluidWrapper(ale));
      }
      else
      {
        FOUR_C_THROW(
            "No ALE adapter available yet for your chosen FSI coupling "
            "algorithm!");
      }
      break;
    }
    case Core::ProblemType::fpsi:
    case Core::ProblemType::fps3i:
    case Core::ProblemType::fsi_xfem:
    case Core::ProblemType::fpsi_xfem:
    {
      ale_ = Teuchos::RCP(new Adapter::AleFpsiWrapper(ale));
      break;
    }
    case Core::ProblemType::fluid_ale:
    case Core::ProblemType::elch:
    case Core::ProblemType::fluid_xfem:
    {
      ale_ = Teuchos::RCP(new Adapter::AleFluidWrapper(ale));
      break;
    }
    default:
      FOUR_C_THROW("ALE type not implemented yet!!");
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
