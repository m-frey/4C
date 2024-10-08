/*----------------------------------------------------------------------*/
/*! \file

\brief testing of FSI calculation results

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_resulttest.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_fluidfluidmonolithic_fluidsplit_nonox.hpp"
#include "4C_fsi_fluidfluidmonolithic_structuresplit_nonox.hpp"
#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_fsi_monolithicstructuresplit.hpp"
#include "4C_fsi_mortarmonolithic_fluidsplit.hpp"
#include "4C_fsi_mortarmonolithic_fluidsplit_sp.hpp"
#include "4C_fsi_mortarmonolithic_structuresplit.hpp"
#include "4C_fsi_slidingmonolithic_fluidsplit.hpp"
#include "4C_fsi_slidingmonolithic_structuresplit.hpp"
#include "4C_io_linedefinition.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <string>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FSIResultTest::FSIResultTest(
    Teuchos::RCP<FSI::Monolithic>& fsi, const Teuchos::ParameterList& fsidyn)
    : Core::UTILS::ResultTest("FSI"), fsi_(fsi)
{
  const auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null) FOUR_C_THROW("Cast to FSI::MonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->fluid_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_monolithicstructuresplit:
    case fsi_iter_fluidfluid_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null) FOUR_C_THROW("Cast to FSI::MonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->structure_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::MortarMonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->fluid_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicfluidsplit_saddlepoint:
    {
      const Teuchos::RCP<FSI::MortarMonolithicFluidSplitSaddlePoint>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicFluidSplitSaddlePoint>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::MortarMonolithicFluidSplitSaddlePoint failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->fluid_field()->discretization();
      auto copy = Teuchos::RCP(new Core::LinAlg::Vector<double>(*fsiobject->lag_mult_));
      copy->ReplaceMap(*fsiobject->fluid_field()->interface()->fsi_cond_map());
      fsilambda_ = copy;

      break;
    }
    case fsi_iter_mortar_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::MortarMonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->structure_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_sliding_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::SlidingMonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::SlidingMonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::SlidingMonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->fluid_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_sliding_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::SlidingMonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::SlidingMonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::SlidingMonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->structure_field()->discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    default:
    {
      slavedisc_ = Teuchos::null;
      fsilambda_ = Teuchos::null;

      std::cout << "\nNo FSI test routines implemented for this coupling algorithm." << std::endl;

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FSIResultTest::FSIResultTest(
    Teuchos::RCP<FSI::MonolithicNoNOX> fsi, const Teuchos::ParameterList& fsidyn)
    : Core::UTILS::ResultTest("FSI")
{
  const auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_fluidfluid_monolithicstructuresplit_nonox:
    {
      // Lagrange multipliers live on the slave field
      slavedisc_ = fsi->structure_field()->discretization();

      const Teuchos::RCP<FSI::FluidFluidMonolithicStructureSplitNoNOX>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::FluidFluidMonolithicStructureSplitNoNOX>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::FluidFluidMonolithicStructureSplitNoNOX failed.");
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_fluidfluid_monolithicfluidsplit_nonox:
    {
      // Lagrange multiplier lives on the slave field (fluid in this case!)
      slavedisc_ = fsi->fluid_field()->discretization();

      const Teuchos::RCP<FSI::FluidFluidMonolithicFluidSplitNoNOX>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::FluidFluidMonolithicFluidSplitNoNOX>(fsi);

      if (fsiobject == Teuchos::null)
        FOUR_C_THROW("Cast to FSI::FluidFluidMonolithicFluidSplitNoNOX failed.");

      fsilambda_ = fsiobject->lambda_;

      break;
    }
    default:
    {
      slavedisc_ = Teuchos::null;
      fsilambda_ = Teuchos::null;

      std::cout << "\nNo FSI test routines implemented for this coupling algorithm." << std::endl;

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  int node = container.get<int>("NODE");
  node -= 1;

  int havenode(slavedisc_->have_global_node(node));
  int isnodeofanybody(0);
  slavedisc_->get_comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW(
        "Node %d does not belong to discretization %s", node + 1, slavedisc_->name().c_str());
  }
  else
  {
    if (slavedisc_->have_global_node(node))
    {
      const Core::Nodes::Node* actnode = slavedisc_->g_node(node);

      // Strange! It seems we might actually have a global node around
      // even if it does not belong to us. But here we are just
      // interested in our nodes!
      if (actnode->owner() != slavedisc_->get_comm().MyPID()) return;

      std::string quantity = container.get<std::string>("QUANTITY");
      bool unknownquantity = true;  // make sure the result value std::string can be handled
      double result = 0.0;          // will hold the actual result of run

      // test Lagrange multipliers
      if (fsilambda_ != Teuchos::null)
      {
        const Epetra_BlockMap& fsilambdamap = fsilambda_->Map();
        if (quantity == "lambdax")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->dof(0, actnode, 0))];
        }
        else if (quantity == "lambday")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->dof(0, actnode, 1))];
        }
        else if (quantity == "lambdaz")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->dof(0, actnode, 2))];
        }
      }

      // catch quantity strings, which are not handled by fsi result test
      if (unknownquantity)
        FOUR_C_THROW("Quantity '%s' not supported in fsi testing", quantity.c_str());

      // compare values
      const int err = compare_values(result, "NODE", container);
      nerr += err;
      test_count++;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::test_element(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  FOUR_C_THROW("FSI ELEMENT test not implemented, yet.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::test_special(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  std::string quantity = container.get<std::string>("QUANTITY");
  bool unknownquantity = true;  // make sure the result value std::string can be handled
  double result = 0.0;          // will hold the actual result of run

  // test for time step size
  if (quantity == "dt")
  {
    unknownquantity = false;
    result = fsi_->dt();
  }

  // test for number of repetitions of time step in case of time step size adaptivity
  if (quantity == "adasteps")
  {
    unknownquantity = false;
    result = fsi_->get_num_adapt_steps();
  }

  // test for simulation time in case of time step size adaptivity
  if (quantity == "time")
  {
    unknownquantity = false;
    result = fsi_->time();
  }

  // catch quantity strings, which are not handled by fsi result test
  if (unknownquantity) FOUR_C_THROW("Quantity '%s' not supported in fsi testing", quantity.c_str());

  // compare values
  const int err = compare_values(result, "SPECIAL", container);
  nerr += err;
  test_count++;
}

FOUR_C_NAMESPACE_CLOSE
