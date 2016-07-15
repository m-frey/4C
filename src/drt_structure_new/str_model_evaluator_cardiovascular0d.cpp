/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_cardiovascular0d.cpp

\brief Evaluation and assembly of all 0D cardiovascular model terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#include "str_model_evaluator_cardiovascular0d.H"

#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Cardiovascular0D::Cardiovascular0D()
    : n_conds_(0),
      disnp_ptr_(Teuchos::null),
      stiff_ptr_(Teuchos::null)
{
  std::cout << "STR::MODELEVALUATOR::Cardiovascular0D::Cardiovascular0D() ############## ..." << std::endl;
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Setup()
{
  std::cout << "Setup ############## ..." << std::endl;

  CheckInit();

  Teuchos::RCP<DRT::Discretization> dis = GState().GetMutableDiscret();

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();

  Teuchos::RCP<LINALG::Solver> dummysolver;

  // initialize 0D cardiovascular manager
  cardvasc0dman_ = Teuchos::rcp(new UTILS::Cardiovascular0DManager(dis,
                                                        disnp_ptr_,
                                                        DRT::Problem::Instance()->StructuralDynamicParams(),
                                                        DRT::Problem::Instance()->Cardiovascular0DStructuralParams(),
                                                        *dummysolver));

  cardvasc0dman_->PrintNewtonHeader();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::ApplyForce(
    const Epetra_Vector& x,
    Epetra_Vector& f)
{
  //std::cout << "ApplyForce ############## ..." << std::endl;
//  CheckInitSetup();
//  Reset(x);
//
//  double time_np = GState().GetTimeNp();
//  Teuchos::ParameterList pwindk;
//  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_, pwindk);
//
//  // --- assemble right-hand-side ---------------------------------------
//  AssembleRhs(f);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::ApplyStiff(const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  //std::cout << "ApplyStiff ############## ..." << std::endl;
//  CheckInitSetup();
//  Reset(x, jac);
//
//  double time_np = GState().GetTimeNp();
//  Teuchos::ParameterList pwindk;
//  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_, pwindk);
//
//  // --- assemble jacobian matrix ---------------------------------------
//  AssembleJacobian(jac);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  //std::cout << "ApplyForceStiff ############## ..." << std::endl;
  CheckInitSetup();
  Reset(x, jac);

  // get structural displacement DOFs

  Teuchos::RCP<Epetra_Vector> block_vec_ptr =
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // get structural stiffness matrix
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd =
      GState().ExtractModelBlock(jac,INPAR::STR::model_cardiovascular0d,
          STR::block_displ_displ);

  //std::cout << "vorher " << *jac_dd << std::endl;
  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pwindk;
  pwindk.set("time_step_size", (*GState().GetDeltaTime())[0]);

  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_, block_vec_ptr, jac_dd, pwindk);
//  cardvasc0dman_->EvaluateForceStiff2(time_np, disnp_ptr_, *block_vec_ptr, *jac_dd, pwindk);

//  jac_dd->Add(*block_ptr,false,1.0,1.0);

  // assemble 0D model contribution to structural rhs (Neumann-like load)
  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

  double foo;
  block_vec_ptr->Norm2(&foo);
//  std::cout << std::setprecision(10)<< "########################## STRUCTUR res-L2-norm: " << foo << std::endl;

//  std::cout << "NACHHER " << *block_ptr << std::endl;
  // reset the block pointer, just to be on the safe side
  block_vec_ptr = Teuchos::null;

  // --- assemble right-hand-side ---------------------------------------
  AssembleRhs(f);
  // --- assemble jacobian matrix ---------------------------------------
  AssembleJacobian(jac);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::AssembleRhs(Epetra_Vector& f) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  // assemble 0D model rhs
  block_vec_ptr = cardvasc0dman_->GetCardiovascular0DRHS();

  if (block_vec_ptr.is_null())
    dserror("The 0D cardiovascular model vector is a NULL pointer, although \n"
        "the structural part indicates, that 0D cardiovascular model contributions \n"
        "are present!");

  //std::cout << f << std::endl;
  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);
  //std::cout << f << std::endl;

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::AssembleJacobian(
    LINALG::SparseOperator& jac) const
{

  Teuchos::RCP<const LINALG::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetMatDstructDcv0ddof();
//  std::cout << *block_ptr << std::endl;
  GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_cardiovascular0d,
      STR::block_displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzd - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetMatDcardvasc0dDd()->Transpose();
//  std::cout << *block_ptr << std::endl;
  GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_cardiovascular0d,
      STR::block_lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetCardiovascular0DStiffness();
//  std::cout << *block_ptr << std::endl;
  GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_cardiovascular0d,
      STR::block_lm_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{

  iowriter.WriteVector("cvdof",
                        cardvasc0dman_->Get0DDofVector());
  iowriter.WriteVector("refvolval",
                        cardvasc0dman_->GetRefVolValue());
  iowriter.WriteVector("reffluxval",
                        cardvasc0dman_->GetRefFluxValue());
  iowriter.WriteVector("refdfluxval",
                        cardvasc0dman_->GetRefDFluxValue());
  iowriter.WriteVector("refddfluxval",
                        cardvasc0dman_->GetRefDDFluxValue());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::ReadRestart(
    IO::DiscretizationReader& ioreader)
{

  double time_n = GState().GetTimeN();
  cardvasc0dman_->ReadRestart(ioreader,time_n);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{

  Teuchos::RCP<Epetra_Vector> cv0d_incr =
      Teuchos::rcp(new Epetra_Vector(*GetBlockDofRowMapPtr()));

  LINALG::Export(dir,*cv0d_incr);

  Teuchos::RCP<Epetra_Vector> dis_incr = GState().ExportModelEntries(INPAR::STR::model_structure,dir);

  cardvasc0dman_->UpdateCv0DDof(cv0d_incr);

  // store for manager-internal monitoring...
  cardvasc0dman_->StoreCv0dDofIncrement(cv0d_incr);
  cardvasc0dman_->StoreStructuralDisplIncrement(dis_incr);

  double norm_0d_disincr;
  cv0d_incr->Norm2(&norm_0d_disincr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::UpdateStepState()
{
  cardvasc0dman_->UpdateTimeStep();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::UpdateStepElement()
{
  // nothing to do
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::DetermineEnergy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Reset(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  //stiff_ptr_ = GState().ExtractDisplBlock(jac);

  CheckInitSetup();
  Reset(x);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  //disnp_ptr_ = GState().GetDisNp();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Cardiovascular0D::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return cardvasc0dman_->GetCardiovascular0DMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Cardiovascular0D::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Cardiovascular0D::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}
