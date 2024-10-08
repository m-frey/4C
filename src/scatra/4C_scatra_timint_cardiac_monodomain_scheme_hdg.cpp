/*----------------------------------------------------------------------*/
/*! \file
\brief time-integration scheme for HDG with extensions for
       cardiac monodomain problems

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_timint_cardiac_monodomain_scheme_hdg.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc_hdg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomainHDG::TimIntCardiacMonodomainHDG(
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntHDG(actdis, solver, sctratimintparams, extraparams, output)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::setup()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntHDG::setup();
  TimIntCardiacMonodomain::setup();

  // Activation time at time n+1
  activation_time_interpol_.reset(new Core::LinAlg::Vector<double>(*discret_->node_row_map()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::update()
{
  // time update of myocard material
  element_material_time_update();

  // Standard Update
  TimIntHDG::update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::element_material_time_update()
{
  discret_->clear_state(true);

  Teuchos::ParameterList eleparams;
  Core::UTILS::add_enum_class_to_parameter_list<ScaTra::Action>(
      "action", ScaTra::Action::time_update_material, eleparams);

  discret_->set_state("phiaf", phinp_);
  discret_->set_state(nds_intvar_, "intphin", intphin_);
  discret_->set_state(0, "phin", phin_);


  Core::LinAlg::SerialDenseMatrix dummyMat;
  Core::LinAlg::SerialDenseVector dummyVec;
  Core::Elements::LocationArray la(discret_->num_dof_sets());


  for (int iele = 0; iele < discret_->num_my_col_elements(); ++iele)
  {
    Core::Elements::Element *ele = discret_->l_col_element(iele);
    ele->location_vector(*discret_, la, false);

    ele->evaluate(eleparams, *discret_, la, dummyMat, dummyMat, dummyVec, dummyVec, dummyVec);
  }

  discret_->clear_state(true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::collect_runtime_output_data()
{
  // call base class first
  TimIntHDG::collect_runtime_output_data();

  if (nb_max_mat_int_state_vars_)
  {
    material_internal_state_np_->PutScalar(0.0);
    Teuchos::ParameterList params;
    Core::UTILS::add_enum_class_to_parameter_list<ScaTra::Action>(
        "action", ScaTra::Action::get_material_internal_state, params);
    params.set<Teuchos::RCP<Epetra_MultiVector>>(
        "material_internal_state", material_internal_state_np_);
    discret_->evaluate(params);
    material_internal_state_np_ =
        params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");
    if (material_internal_state_np_ == Teuchos::null)
      FOUR_C_THROW("Cannot get state vector material internal state");

    std::vector<std::optional<std::string>> context(
        material_internal_state_np_->NumVectors(), "ionic_currents");
    visualization_writer().append_result_data_vector_with_context(
        *material_internal_state_np_, Core::IO::OutputEntity::element, context);

    for (int k = 0; k < material_internal_state_np_->NumVectors(); ++k)
    {
      std::ostringstream temp;
      temp << k + 1;
      material_internal_state_np_component_ =
          Teuchos::RCP(new Core::LinAlg::Vector<double>(*(*material_internal_state_np_)(k)));

      visualization_writer().append_result_data_vector_with_context(
          *material_internal_state_np_component_, Core::IO::OutputEntity::element,
          {"mat_int_state" + temp.str()});
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::write_restart() const
{
  // step number and time (only after that data output is possible)
  output_->new_step(step_, time_);

  // output restart information associated with mesh tying strategy
  strategy_->write_restart();

  output_->write_vector("intphinp", intphinp_);
  output_->write_vector("phinp_trace", phinp_);
  output_->write_vector("intphin", intphin_);

  // copy values from node to dof vector
  Teuchos::RCP<Core::LinAlg::Vector<double>> dofphi =
      Core::LinAlg::create_vector(*discret_->node_row_map());

  for (int i = 0; i < dofphi->MyLength(); ++i)
  {
    int dofgid = discret_->node_row_map()->GID(i);
    dofphi->ReplaceMyValue(discret_->node_row_map()->LID(dofgid), 0, (*interpolatedPhinp_)[i]);
  }

  output_->write_vector("phinp", dofphi);

  // add info to control file for reading all variables in restart
  output_->write_mesh(step_, time_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::collect_problem_specific_runtime_output_data(
    Teuchos::RCP<Core::LinAlg::Vector<double>> interpolatedPhi)
{
  // Compute and write activation time
  if (activation_time_interpol_ != Teuchos::null)
  {
    for (int k = 0; k < interpolatedPhi->MyLength(); k++)
    {
      if ((*interpolatedPhi)[k] >= activation_threshold_ &&
          (*activation_time_interpol_)[k] <= dta_ * 0.9)
        (*activation_time_interpol_)[k] = time_;
    }
    visualization_writer().append_result_data_vector_with_context(
        *activation_time_interpol_, Core::IO::OutputEntity::node, {"activation_time"});
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::pack_material()
{
  Core::Communication::PackBuffer buffer;

  // loop over elements
  for (int iele = 0; iele < discret_->num_my_col_elements(); ++iele)
  {
    auto *hdgele = dynamic_cast<Discret::ELEMENTS::ScaTraHDG *>(discret_->l_col_element(iele));
    hdgele->pack_material(buffer);
  }

  Teuchos::RCP<std::vector<char>> block = Teuchos::RCP(new std::vector<char>);
  std::swap(*block, buffer());
  data_ = block;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::unpack_material()
{
  // loop over elements
  std::vector<char> data;
  Core::Communication::UnpackBuffer buffer(*data_);
  for (int iele = 0; iele < discret_->num_my_col_elements(); ++iele)
  {
    auto *hdgele = dynamic_cast<Discret::ELEMENTS::ScaTraHDG *>(discret_->l_col_element(iele));
    extract_from_pack(buffer, data);
    Core::Communication::UnpackBuffer buffer_mat(data);
    hdgele->unpack_material(buffer_mat);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::project_material()
{
  discret_->clear_state(true);
  // set action
  Teuchos::ParameterList eleparams;
  Core::UTILS::add_enum_class_to_parameter_list<ScaTra::Action>(
      "action", ScaTra::Action::project_material_field, eleparams);

  Core::LinAlg::SerialDenseMatrix dummyMat;
  Core::LinAlg::SerialDenseVector dummyVec;
  Core::Elements::LocationArray dummy(1);

  for (int iele = 0; iele < discret_->num_my_col_elements(); ++iele)
  {
    Core::Elements::Element *ele = discret_->l_col_element(iele);

    // call routine on elements to project material field
    ele->evaluate(eleparams, *discret_, dummy, dummyMat, dummyMat, dummyVec, dummyVec, dummyVec);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainHDG::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // Call function from base class
  ScaTra::TimIntHDG::read_restart(step, input);

  activation_time_interpol_.reset(new Core::LinAlg::Vector<double>(*discret_->node_row_map()));
}

FOUR_C_NAMESPACE_CLOSE
