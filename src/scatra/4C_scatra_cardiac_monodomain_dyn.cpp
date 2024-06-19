/*----------------------------------------------------------------------*/
/*! \file

\brief entry point for cardiac monodomain scalar transport problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_cardiac_monodomain_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_binstrategy.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_resulttest.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, incl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditions)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *        o scalar transport in velocity field given by Navier-Stokes with natural convection
 *          (two-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_cardiac_monodomain_dyn(int restart)
{
  // pointer to problem
  Global::Problem* problem = Global::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

  //  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << Global::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }


  // print CardiacMonodomain-Logo to screen
  if (comm.MyPID() == 0) printheartlogo();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  // access the fluid discretization
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra dof
  fluiddis->fill_complete();
  scatradis->fill_complete();

  // set velocity field
  const Inpar::ScaTra::VelocityField veltype =
      Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(scatradyn, "VELOCITYFIELD");
  switch (veltype)
  {
    case Inpar::ScaTra::velocity_zero:      // zero  (see case 1)
    case Inpar::ScaTra::velocity_function:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes() == 0)
        FOUR_C_THROW("No elements in the ---TRANSPORT ELEMENTS section");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
      {
        FOUR_C_THROW(
            "no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in "
            "SCALAR TRANSPORT DYNAMIC to a valid number!");
      }

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatraonly =
          Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
              scatradyn, scatradyn, Global::Problem::Instance()->SolverParams(linsolvernumber)));

      // add proxy of velocity related degrees of freedom to scatra discretization
      Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
          Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
              Global::Problem::Instance()->NDim() + 1, 0, 0, true));
      if (scatradis->AddDofSet(dofsetaux) != 1)
        FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      scatraonly->ScaTraField()->set_number_of_dof_set_velocity(1);

      // allow TRANSPORT conditions, too
      // NOTE: we can not use the conditions given by 'conditions_to_copy =
      // clonestrategy.conditions_to_copy()' since we than may have some scatra condition twice. So
      // we only copy the Dirichlet and Neumann conditions:
      const std::map<std::string, std::string> conditions_to_copy = {
          {"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
          {"TransportLineNeumann", "LineNeumann"}, {"TransportSurfaceNeumann", "SurfaceNeumann"},
          {"TransportVolumeNeumann", "VolumeNeumann"}};

      Core::FE::DiscretizationCreatorBase creator;
      creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

      // finalize discretization
      scatradis->fill_complete();



      // We have do use the binningstrategy and extended ghosting when we use p-adaptivity. This
      // guarantees that the elements at the boarder between the processor calculate correctly since
      // one face is shared with the neighboring element (which is owned by an other processor =
      // ghosted element) which again is sharing other faces with elements on other processors
      // (extended ghosted element)
      if (Core::UTILS::IntegralValue<bool>(scatradyn, "PADAPTIVITY"))
      {
        // redistribute discr. with help of binning strategy
        if (scatradis->Comm().NumProc() > 1)
        {
          // create vector of discr.
          std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
          dis.push_back(scatradis);

          // binning strategy for parallel redistribution
          Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy;

          std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
          std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

          // binning strategy is created and parallel redistribution is performed
          Teuchos::ParameterList binning_params =
              Global::Problem::Instance()->binning_strategy_params();
          Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
              "spatial_approximation_type",
              Global::Problem::Instance()->spatial_approximation_type(), binning_params);
          binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(binning_params,
              Global::Problem::Instance()->OutputControlFile(), scatradis->Comm(),
              scatradis->Comm().MyPID(), dis));
          binningstrategy
              ->do_weighted_partitioning_of_bins_and_extend_ghosting_of_discret_to_one_bin_layer(
                  dis, stdelecolmap, stdnodecolmap);
        }
      }

      // now we can call init() on the base algo.
      // time integrator is constructed and initialized inside
      scatraonly->init();

      // NOTE : At this point we may redistribute and/or
      //        ghost our discretizations at will.

      // now we must call setup()
      scatraonly->setup();

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->read_restart(restart);

      // set initial velocity field
      // note: The order read_restart() before set_velocity_field() is important here!!
      // for time-dependent velocity fields, set_velocity_field() is additionally called in each
      // prepare_time_step()-call
      (scatraonly->ScaTraField())->set_velocity_field();

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField())->TimeLoop();

      // perform the result test if required
      Global::Problem::Instance()->AddFieldTest(scatraonly->create_sca_tra_field_test());
      Global::Problem::Instance()->TestAll(comm);

      break;
    }
    case Inpar::ScaTra::velocity_Navier_Stokes:  // Navier_Stokes
    {
      FOUR_C_THROW(
          "Navier Stokes case not implemented for cardiac monodomain scalar transport problem");

      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      const Inpar::ScaTra::FieldCoupling fieldcoupling =
          Core::UTILS::IntegralValue<Inpar::ScaTra::FieldCoupling>(
              Global::Problem::Instance()->scalar_transport_dynamic_params(), "FIELDCOUPLING");

      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes() == 0)
      {
        if (fieldcoupling != Inpar::ScaTra::coupling_match)
          FOUR_C_THROW(
              "If you want matching fluid and scatra meshes, do clone you fluid mesh and use "
              "FIELDCOUPLING match!");

        fluiddis->fill_complete();
        scatradis->fill_complete();

        // fill scatra discretization by cloning fluid discretization
        Core::FE::CloneDiscretization<ScaTra::ScatraFluidCloneStrategy>(
            fluiddis, scatradis, Global::Problem::Instance()->CloningMaterialMap());

        // set implementation type of cloned scatra elements
        for (int i = 0; i < scatradis->NumMyColElements(); ++i)
        {
          Discret::ELEMENTS::Transport* element =
              dynamic_cast<Discret::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if (element == nullptr)
            FOUR_C_THROW("Invalid element type!");
          else
            element->SetImplType(Inpar::ScaTra::impltype_std);
        }

        // add proxy of fluid transport degrees of freedom to scatra discretization
        if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
          FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      }
      else
      {
        if (fieldcoupling != Inpar::ScaTra::coupling_volmortar)
          FOUR_C_THROW(
              "If you want non-matching fluid and scatra meshes, you need to use FIELDCOUPLING "
              "volmortar!");

        // allow TRANSPORT conditions, too
        ScaTra::ScatraFluidCloneStrategy clonestrategy;
        const auto conditions_to_copy = clonestrategy.conditions_to_copy();
        Core::FE::DiscretizationCreatorBase creator;
        creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

        // first call fill_complete for single discretizations.
        // This way the physical dofs are numbered successively
        fluiddis->fill_complete();
        scatradis->fill_complete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
        const int ndofperelement_scatra = 0;
        const int ndofpernode_fluid = fluiddis->NumDof(0, fluiddis->lRowNode(0));
        const int ndofperelement_fluid = 0;
        Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (fluiddis->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in fluid field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");

        // call assign_degrees_of_freedom also for auxiliary dofsets
        // note: the order of fill_complete() calls determines the gid numbering!
        // 1. fluid dofs
        // 2. scatra dofs
        // 3. fluid auxiliary dofs
        // 4. scatra auxiliary dofs
        fluiddis->fill_complete(true, false, false);
        scatradis->fill_complete(true, false, false);

        // NOTE: we have do use the binningstrategy here since we build our fluid and scatra
        // problems by inheritance, i.e. by calling the constructor of the corresponding class. But
        // since we have to use the binning-strategy before creating the single field we have to do
        // it here :-( We would prefer to to it like the SSI since than we could extended ghosting
        // TODO (thon): make this if-case obsolete and allow for redistribution within
        // volmortar->setup() by removing inheitance-building of fields
        {
          // redistribute discr. with help of binning strategy
          if (fluiddis->Comm().NumProc() > 1)
          {
            // create vector of discr.
            std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
            dis.push_back(fluiddis);
            dis.push_back(scatradis);

            // binning strategy for parallel redistribution
            Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy;

            std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
            std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

            /// binning strategy is created and parallel redistribution is performed
            Teuchos::ParameterList binning_params =
                Global::Problem::Instance()->binning_strategy_params();
            Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
                "spatial_approximation_type",
                Global::Problem::Instance()->spatial_approximation_type(), binning_params);
            binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(binning_params,
                Global::Problem::Instance()->OutputControlFile(), scatradis->Comm(),
                scatradis->Comm().MyPID(), dis));
            binningstrategy
                ->do_weighted_partitioning_of_bins_and_extend_ghosting_of_discret_to_one_bin_layer(
                    dis, stdelecolmap, stdnodecolmap);
          }
        }
      }

      // support for turbulent flow statistics
      const Teuchos::ParameterList& fdyn = (Global::Problem::Instance()->FluidDynamicParams());

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        FOUR_C_THROW(
            "no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in "
            "SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create a scalar transport algorithm instance
      Teuchos::RCP<ScaTra::ScaTraAlgorithm> algo = Teuchos::rcp(new ScaTra::ScaTraAlgorithm(comm,
          scatradyn, fdyn, "scatra", Global::Problem::Instance()->SolverParams(linsolvernumber)));

      // init algo (init fluid time integrator and scatra time integrator inside)
      algo->init();

      // setup algo (setup fluid time integrator and scatra time integrator inside)
      algo->setup();

      // read restart information
      // in case a inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        if ((Core::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
                true) and
            (restart == fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
          algo->ReadInflowRestart(restart);
        else
          algo->read_restart(restart);
      }
      else if (Core::UTILS::IntegralValue<int>(
                   fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") == true)
        FOUR_C_THROW(
            "Turbulent inflow generation for passive scalar transport should be performed as fluid "
            "problem!");

      // solve the whole scalar transport problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      Global::Problem::Instance()->AddFieldTest(algo->fluid_field()->CreateFieldTest());
      Global::Problem::Instance()->AddFieldTest(algo->create_sca_tra_field_test());
      Global::Problem::Instance()->TestAll(comm);

      break;
    }  // case 2
    default:
    {
      FOUR_C_THROW("unknown velocity field type for transport of passive scalar");
      break;
    }
  }

  return;

}  // end of scatra_cardiac_monodomain_dyn()


/*----------------------------------------------------------------------*/
// print ELCH-Module logo
/*----------------------------------------------------------------------*/
void printheartlogo()
{
  // more at http://www.ascii-art.de
  std::cout << "                                                         " << std::endl;
  std::cout << "               |  \\ \\ | |/ /                           " << std::endl;
  std::cout << "               |  |\\ `' ' /                             " << std::endl;
  std::cout << "               |  ;'aorta \\      / , pulmonary          " << std::endl;
  std::cout << "               | ;    _,   |    / / ,  arteries          " << std::endl;
  std::cout << "      superior | |   (  `-.;_,-' '-' ,                   " << std::endl;
  std::cout << "     vena cava | `,   `-._       _,-'_                   " << std::endl;
  std::cout
      << "               |,-`.    `.)    ,<_,-'_, pulmonary                     ______ _____    "
      << std::endl;
  std::cout
      << "              ,'    `.   /   ,'  `;-' _,  veins                      |  ____|  __ \\   "
      << std::endl;
  std::cout
      << "             ;        `./   /`,    \\-'                               | |__  | |__) |  "
      << std::endl;
  std::cout
      << "             | right   /   |  ;\\   |\\                                |  __| |  ___/   "
      << std::endl;
  std::cout
      << "             | atrium ;_,._|_,  `, ' \\                               | |____| |       "
      << std::endl;
  std::cout
      << "             |        \\    \\ `       `,                              |______|_|       "
      << std::endl;
  std::cout
      << "             `      __ `    \\   left  ;,                                              "
      << std::endl;
  std::cout
      << "              \\   ,'  `      \\,  ventricle                                            "
      << std::endl;
  std::cout << "               \\_(            ;,      ;;                " << std::endl;
  std::cout << "               |  \\           `;,     ;;                " << std::endl;
  std::cout << "      inferior |  |`.          `;;,   ;'                 " << std::endl;
  std::cout << "     vena cava |  |  `-.        ;;;;,;'                  " << std::endl;
  std::cout << "               |  |    |`-.._  ,;;;;;'                   " << std::endl;
  std::cout << "               |  |    |   | ``';;;'                     " << std::endl;
  std::cout << "                       aorta                             " << std::endl;
  std::cout << "                                                         " << std::endl;
}

FOUR_C_NAMESPACE_CLOSE
