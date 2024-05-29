/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::
    BeamContactAssemblyManagerInDirect(
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>&
            assembly_contact_elepairs,
        const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& beam_to_solid_params)
    : BeamContactAssemblyManager()
{
  // Create the mortar manager. We add 1 to the MaxAllGID since this gives the maximum GID and NOT
  // the length of the GIDs.
  mortar_manager_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidMortarManager>(
      new BEAMINTERACTION::BeamToSolidMortarManager(
          discret, beam_to_solid_params, discret->dof_row_map()->MaxAllGID() + 1));

  // Setup the mortar manager.
  mortar_manager_->Setup();
  mortar_manager_->SetLocalMaps(assembly_contact_elepairs);
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::evaluate_force_stiff(
    Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<CORE::LINALG::SparseMatrix> fe_sysmat)
{
  // Evaluate the global mortar matrices.
  mortar_manager_->evaluate_global_coupling_contributions(data_state->GetDisColNp());

  // Add the global mortar matrices to the force vector and stiffness matrix.
  mortar_manager_->add_global_force_stiffness_penalty_contributions(
      data_state, fe_sysmat, fe_sysvec);
}


double BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::get_energy(
    const Teuchos::RCP<const Epetra_Vector>& disp) const
{
  return mortar_manager_->get_energy();
}

FOUR_C_NAMESPACE_CLOSE
