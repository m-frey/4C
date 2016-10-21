/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_beaminteraction.cpp

\brief Evaluation of all beam interaction terms

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_beaminteraction.H"
#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"

#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_beamcontact/beam3contact_manager.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteraction::BeamInteraction()
    : disnp_ptr_(Teuchos::null),
      stiff_beaminteract_ptr_(Teuchos::null),
      f_beaminteract_np_ptr_(Teuchos::null),
      beamcman_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // setup the pointers for displacement and stiffness
  disnp_ptr_ = GState().GetMutableDisNp();
  stiff_beaminteract_ptr_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *GState().DofRowMapView(), 81, true, true));
  f_beaminteract_np_ptr_ =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // create beam contact manager
  beamcman_ = Teuchos::rcp(new CONTACT::Beam3cmanager(Discret(),0.0));

  // gmsh output at beginning of simulation
#ifdef GMSHTIMESTEPS
  beamcman_->GmshOutput(*disnp_ptr_, 0, 0, true);
#endif

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  disnp_ptr_ = GState().GetDisNp();

  // Zero out force and stiffness contributions
  f_beaminteract_np_ptr_->PutScalar(0.0);
  stiff_beaminteract_ptr_->Zero();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateForce()
{
  CheckInitSetup();

  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", EvalData().GetNlnIter());
  beamcontactparams.set("dt", EvalData().GetDeltaTime());
  beamcontactparams.set("numstep", EvalData().GetStepNp());

  beamcman_->Evaluate(*stiff_beaminteract_ptr_,*f_beaminteract_np_ptr_,
      *disnp_ptr_,beamcontactparams,true,EvalData().GetTotalTime());

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateStiff()
{
  CheckInitSetup();


  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", EvalData().GetNlnIter());
  beamcontactparams.set("dt", EvalData().GetDeltaTime());
  beamcontactparams.set("numstep", EvalData().GetStepNp());

  beamcman_->Evaluate(*stiff_beaminteract_ptr_,*f_beaminteract_np_ptr_,
      *disnp_ptr_,beamcontactparams,true,EvalData().GetTotalTime());

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateForceStiff()
{
  CheckInitSetup();

  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", EvalData().GetNlnIter());
  beamcontactparams.set("dt", EvalData().GetDeltaTime());
  beamcontactparams.set("numstep", EvalData().GetStepNp());

  beamcman_->Evaluate(*stiff_beaminteract_ptr_,*f_beaminteract_np_ptr_,
      *disnp_ptr_,beamcontactparams,true,EvalData().GetTotalTime());

  // visualization of current Newton step
#ifdef GMSHNEWTONSTEPS
  beamcman_->GmshOutput(*disnp_ptr_,EvalData().GetStepNp(),EvalData().GetNlnIter());
  beamcman_->ConsoleOutput();
#endif // #ifdef GMSHNEWTONSTEPS

  // update constraint norm
  beamcman_->UpdateConstrNorm();    // ToDo

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  STR::AssembleVector(1.0,f,-timefac_np,*f_beaminteract_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_beaminteract_ptr_,false,timefac_np,1.0);
  // no need to keep it
  stiff_beaminteract_ptr_->Zero();
  // nothing to do
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{
  beamcman_->WriteRestart(iowriter);    // ToDo

  // since the global OutputStepState() routine is not called, if the
  // restart is written, we have to do it here manually.
  OutputStepState(iowriter);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  beamcman_->ReadRestart(ioreader);   //ToDo
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  // empty ToDo
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateStepState(
    const double& timefac_n)
{
  beamcman_->Update(*disnp_ptr_, EvalData().GetStepNp(), EvalData().GetNlnIter());

  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr =
      GState().GetMutableFstructureOld();

  fstructold_ptr->Update(-timefac_n,*f_beaminteract_np_ptr_,1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateStepElement()
{
  // empty
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::DetermineEnergy()
{
  // ToDo
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::ResetStepState()
{

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::BeamInteraction::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::BeamInteraction::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::BeamInteraction::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::PostOutput()
{
  CheckInitSetup();
 // empty

  return;
} // PostOutput()
