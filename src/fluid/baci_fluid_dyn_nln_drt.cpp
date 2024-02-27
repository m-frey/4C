/*-----------------------------------------------------------*/
/*! \file

\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_dyn_nln_drt.hpp"

#include "baci_adapter_fld_base_algorithm.hpp"
#include "baci_fluid_turbulence_turbulent_flow_algorithm.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_condition_utils.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha (two versions)
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void dyn_fluid_drt(const int restart)
{
  // create a communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("fluid")->Comm();

  // access to some parameter lists
  // const Teuchos::ParameterList& probtype = GLOBAL::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();

  // prepares a turbulent flow simulation with generation of turbulent inflow during the
  // actual simulation
  // this is done in two steps
  // 1. computation of inflow until it reaches a fully turbulent state
  // 2. computation of the main problem after restart
  // Remark: we restart the simulation to save procs!
  if ((CORE::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
          true) and
      (restart < fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
  {
    if (comm.MyPID() == 0)
    {
      std::cout << "#-----------------------------------------------#" << std::endl;
      std::cout << "#      ENTER TURBULENT INFLOW COMPUTATION       #" << std::endl;
      std::cout << "#-----------------------------------------------#" << std::endl;
    }

    // create instance of fluid turbulent flow algorithm
    Teuchos::RCP<FLD::TurbulentFlowAlgorithm> turbfluidalgo =
        Teuchos::rcp(new FLD::TurbulentFlowAlgorithm(comm, fdyn));

    // read the restart information, set vectors and variables
    if (restart) turbfluidalgo->ReadRestart(restart);

    // run simulation for a separate part of the complete domain to get turbulent flow in it
    // after restart a turbulent inflow profile is computed in the separate inflow section and
    // transferred as dirichlet boundary condition to the problem domain of interest
    // this finally allows to get high quality turbulent inflow conditions during simulation of the
    // actual flow
    turbfluidalgo->TimeLoop();

    // perform result tests if required
    GLOBAL::Problem::Instance()->AddFieldTest(turbfluidalgo->DoResultCheck());
    GLOBAL::Problem::Instance()->TestAll(comm);
  }
  // solve a simple fluid problem
  else
  {
    // create instance of fluid basis algorithm
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo =
        Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

    // read the restart information, set vectors and variables
    if (restart) fluidalgo->FluidField()->ReadRestart(restart);

    // run the simulation
    //    fluidalgo->FluidField()->TimeLoop();
    fluidalgo->FluidField()->Integrate();

    // perform result tests if required
    GLOBAL::Problem::Instance()->AddFieldTest(fluidalgo->FluidField()->CreateFieldTest());
    GLOBAL::Problem::Instance()->TestAll(comm);
  }

  // have fun with your results!
  return;

}  // end of dyn_fluid_drt()

BACI_NAMESPACE_CLOSE
