// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timeloop.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_stru_multi_microstatic.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"
#include <4C_solver_nonlin_nox_problem.hpp>

#include <boost/math/special_functions/math_fwd.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimeLoop::integrate()
{
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

      /*  =========================================================================================
       *              STATIC HOMOGENIZAITON
       *  =========================================================================================
       */
      {
        auto problem = Global::Problem::instance();
        auto discret_ = problem->get_dis("structure");
        const Epetra_Map* dofrowmap = discret_->dof_row_map();
        auto stiff_ = std::make_shared<Core::LinAlg::SparseMatrix>(*dofrowmap, 81, true, true);
        auto dt_ = 0.;
        auto timen_ = 1.;

        auto dispn_ = dispn();
        std::shared_ptr<Core::LinAlg::Vector<double>> disn_ =
            std::make_shared<Core::LinAlg::Vector<double>>(*dispn_);
        std::shared_ptr<Core::LinAlg::Vector<double>> disi_ =
            std::make_shared<Core::LinAlg::Vector<double>>(*dispn_);
        std::shared_ptr<Core::LinAlg::Vector<double>> fint_ =
            std::make_shared<Core::LinAlg::Vector<double>>(*freact());
        std::shared_ptr<Core::LinAlg::Vector<double>> fintn_ =
            std::make_shared<Core::LinAlg::Vector<double>>(*freact());
        fintn_->PutScalar(0.0);
        disn_->PutScalar(0.0);

        //---------------------------- compute internal forces and stiffness
        // zero out stiffness
        stiff_->zero();
        // create the parameters for the discretization
        Teuchos::ParameterList p;
        // action for elements
        p.set("action", "calc_struct_nlnstiff");
        // other parameters that might be needed by the elements
        p.set("total time", timen_);
        p.set("delta time", dt_);
        // set vector values needed by elements
        discret_->clear_state();
        // we do not need to scale disi_ here with 1-alphaf (cf. strugenalpha), since
        // everything on the microscale "lives" at the pseudo generalized midpoint
        // -> we solve our quasi-static problem there and only update data to the "end"
        // of the time step after having finished a macroscopic dt
        discret_->set_state("residual displacement", disi_);
        discret_->set_state("displacement", disn_);

        // std::cout << "\ndisi_ (residual displacement:) " << disi_-> << std::endl;
        // std::cout << "\ndisn_ (displacement:) " << *disn_ << std::endl;


        fintn_->PutScalar(0.0);  // initialise internal force vector
        discret_->evaluate(p, stiff_, nullptr, fintn_, nullptr, nullptr);

        // std::cout << "\nFINT (AFTER evaluate NEWTIMINT) = " << *fintn_ << std::endl;
        // std::cout << "\nstiff_ (AFTER evaluate NEWTIMINT)" << std::endl;
        // stiff_->epetra_matrix()->Print(std::cout);
        Core::LinAlg::print_matrix_in_matlab_format("stiff_in_nti", *stiff_->epetra_matrix(), true);
        discret_->clear_state();
      }
      auto dispn_ = dispn();
      // std::cout << "====== dispn is: (NEW TIMEINT) =====\n" << *dispn_;

      auto fr = freact();
      // std::cout << "\n====== freact is (IN NEW TIMINT) ========\n" << *fr;

      // static homogen. as OUTPUT:
      auto MicroStatic_ = Teuchos::rcp(new MultiScale::MicroStatic(0, 1.0, true));
      // MicroStatic_->import_test_freat();
      MicroStatic_->import_freact(freact());

      // Use reaktion force from minimal test just to check everything works as planed:
      // std::cout << "\n Manul def frext \n" << *MicroStatic_->freactn_;


      // overwrite reaction force with result from the
      // MicroStatic_->freactn_ = freact_from_micro;

      Core::LinAlg::Matrix<6, 1> stress(true);
      Core::LinAlg::Matrix<6, 6> cmat(true);

      // ===================================================
      Core::LinAlg::Matrix<3, 3> defgrd(true);

      // defgrd(0, 0) = 0.977655;
      // defgrd(0, 1) = 2.61724e-17;
      // defgrd(0, 2) = 1.49806e-17;
      // defgrd(1, 0) = 2.32727e-17;
      // defgrd(1, 1) = 0.977655;
      // defgrd(1, 2) = 2.63945e-17;
      // defgrd(2, 0) = -1.11022e-16;
      // defgrd(2, 1) = 0;
      // defgrd(2, 2) = 1.09394;

      // from 520 Fe2 results
      // Fs = F + df * beta

      const Teuchos::ParameterList& sdyn_macro =
          Global::Problem::instance()->structural_dynamic_params();

      auto maxtime = sdyn_macro.get<double>("MAXTIME");
      std::cout << "\nmaxtime used for scaling: " << maxtime << " \n";
      double endTime = 1.0;
      double beta = (1. - time_old() / endTime);
      std::cout << "\n beta used for scaling:" << beta << " \n";


      const double F_11 = 0.977655;
      const double F_12 = 1.00491e-17;
      const double F_13 = 2.47629e-18;
      const double F_21 = 4.28386e-17;
      const double F_22 = 0.977655;
      const double F_23 = -2.61468e-17;
      const double F_31 = 0.0;
      const double F_32 = 1.11022e-16;
      const double F_33 = 1.09394;

      defgrd(0, 0) = F_11 + beta * (1. - F_11);
      std::cout << "beta * (1 - F11)" << beta * (1. - F_11) << std::endl;
      defgrd(0, 1) = -1.00491e-17 - beta * (F_12);
      defgrd(0, 2) = 2.47629e-18 - beta * (F_13);
      defgrd(1, 0) = 4.28386e-17 - beta * (F_21);
      defgrd(1, 1) = 0.977655 + beta * (1. - F_22);
      defgrd(1, 2) = -2.61468e-17 - beta * (F_23);
      defgrd(2, 0) = 0.0 - beta * (F_31);
      defgrd(2, 1) = 1.11022e-16 - beta * (F_32);
      defgrd(2, 2) = 1.09394 + beta * (1. - F_33);

      // defgrd(0, 0) = 0.9;  // 0.977655;
      // defgrd(0, 1) = 0.0;  // 4.90714e-17;
      // defgrd(0, 2) = 0.0;  //-4.81204e-18;
      // defgrd(1, 0) = 0.0;  // 7.62351e-17;
      // defgrd(1, 1) = 0.9;  // 0.977655;
      // defgrd(1, 2) = 0.0;  //-4.89088e-17;
      // defgrd(2, 1) = 0.0;  // 1.11022e-16;
      // defgrd(2, 1) = 0.0;  // 0.0;
      // defgrd(2, 2) = 0.9;  // 1.09394;

      // ==== scale def grad:

      // Get defgrad from this section

      auto dt = sdyn_macro.get<double>("TIMESTEP");


      std::cout << "\n time now, it is: time()-dt =  " << time() - dt << "\n";
      std::cout << "\n time now, it is: structure time_old =  " << structure_->time_old() << "\n";
      std::cout << "\n time now, it is: structure time =  " << structure_->time() << "\n";
      // =====================================================
      const bool mod_newton = false;
      bool build_stiff = true;
      MicroStatic_->static_homogenization(&stress, &cmat, &defgrd, mod_newton, build_stiff);
      /*  =========================================================================================
       *              STATIC HOMOGENIZAITON END
       *  =========================================================================================
       */


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
