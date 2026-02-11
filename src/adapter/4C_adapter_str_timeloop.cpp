// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timeloop.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_tensor.hpp"

#include <boost/math/special_functions/math_fwd.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

std::shared_ptr<Core::LinAlg::SparseMatrix> global_fullstiff = nullptr;
;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimeLoop::integrate()
{
  // Print optional user-provided macro deformation gradient once, if present
  {
    static bool defgrad_printed = false;
    if (!defgrad_printed)
    {
      const auto& constraint = Global::Problem::instance()->constraint_params();
      if (constraint.isParameter("F_MACRO"))
      {
        const auto defgrad = constraint.get<Core::LinAlg::Tensor<double, 3, 3>>("F_MACRO");
        std::cout << "CONSTRAINT/F_MACRO:" << std::endl;
        for (int i = 0; i < 3; ++i)
        {
          std::cout << "  ";
          for (int j = 0; j < 3; ++j)
          {
            std::cout << defgrad(i, j) << (j < 2 ? " " : "");
          }
          std::cout << std::endl;
        }
      }
      defgrad_printed = true;
    }
  }

  // error checking variables
  Inpar::Solid::ConvergenceStatus convergencestatus = Inpar::Solid::conv_success;

  // target time #timen_ and step #stepn_ already set
  // time loop
  while (not_finished() and (convergencestatus == Inpar::Solid::conv_success or
                                convergencestatus == Inpar::Solid::conv_fail_repeat))
  {
    // call the predictor
    pre_predict();
    prepare_time_step();

    // integrate time step, i.e. do corrector steps
    // after this step we hold disn_, etc
    pre_solve();
    convergencestatus = solve();

    // if everything is fine
    if (convergencestatus == Inpar::Solid::conv_success)
    {
      // calculate stresses, strains and energies
      // note: this has to be done before the update since otherwise a potential
      // material history is overwritten
      constexpr bool force_prepare = false;
      prepare_output(force_prepare);

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      // update time and step
      // update everything on the element level
      pre_update();
      update();
      post_update();

      // write output
      output();
      post_output();

      // print info about finished time step
      print_step();
    }
    // todo: remove this as soon as old structure time integration is gone
    else if (Teuchos::getIntegralValue<Inpar::Solid::IntegrationStrategy>(
                 Global::Problem::instance()->structural_dynamic_params(), "INT_STRATEGY") ==
             Inpar::Solid::int_old)
    {
      convergencestatus =
          perform_error_action(convergencestatus);  // something went wrong update error code
                                                    // according to chosen divcont action
    }
  }

  post_time_loop();

  // that's it say what went wrong
  return convergencestatus;
}

FOUR_C_NAMESPACE_CLOSE
