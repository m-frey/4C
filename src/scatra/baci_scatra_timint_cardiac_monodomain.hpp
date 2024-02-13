/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for cardiac monodomain

\level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_CARDIAC_MONODOMAIN_HPP
#define BACI_SCATRA_TIMINT_CARDIAC_MONODOMAIN_HPP

#include "baci_config.hpp"

#include "baci_scatra_timint_implicit.hpp"

BACI_NAMESPACE_OPEN


/*==========================================================================*/
// forward declarations
/*==========================================================================*/


namespace SCATRA
{
  class TimIntCardiacMonodomain : public virtual ScaTraTimIntImpl
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntCardiacMonodomain(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! setup algorithm
    void Setup() override;

    //! time update of time-dependent materials
    void ElementMaterialTimeUpdate();

    void OutputState() override;

    //! Set ep-specific parameters
    void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const override;

    //! return pointer to old scatra dofrowmap for coupled problems
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapScatra() { return dofmap_; };

    /*========================================================================*/
    //! @name electrophysiology
    /*========================================================================*/

    //! activation_time at times n+1
    Teuchos::RCP<Epetra_Vector> activation_time_np_;

    //! activation threshold for postprocessing
    double activation_threshold_;

    //! maximum expected number of material internal state variables
    int nb_max_mat_int_state_vars_;

    //! material internal state at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_internal_state_np_;

    //! one component of the material internal state at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_internal_state_np_component_;

    //! maximum expected number of material ionic currents
    int nb_max_mat_ionic_currents_;

    //! material ionic currents at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_ionic_currents_np_;

    //! one component of the material ionic currents at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_ionic_currents_np_component_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> ep_params_;

    //! row map of dofs (one dof for each node)
    Teuchos::RCP<Epetra_Map> dofmap_;
  };

};  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif
