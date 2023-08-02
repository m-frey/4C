/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_poro.H"

#include "baci_contact_interface.H"
#include "baci_contact_nitsche_utils.H"
#include "baci_contact_paramsinterface.H"
#include "baci_coupling_adapter.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_so3_plast_ssn.H"

#include <Epetra_FEVector.h>
#include <Epetra_Operator.h>

void CONTACT::CoNitscheStrategyPoro::ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
    const int iter, bool predictor)
{
  if (predictor) return;

  CONTACT::CoNitscheStrategy::ApplyForceStiffCmt(dis, kt, f, step, iter, predictor);

  // Evaluation for all interfaces
  fp_ = CreateRhsBlockPtr(DRT::UTILS::VecBlockType::porofluid);
  kpp_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::porofluid_porofluid);
  kpd_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::porofluid_displ);
  kdp_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_porofluid);
  //    for (int i = 0; i < (int) interface_.size(); ++i)
  //    {
  //      for (int e=0;e<interface_[i]->Discret().ElementColMap()->NumMyElements();++e)
  //      {
  //        MORTAR::MortarElement* mele
  //        =dynamic_cast<MORTAR::MortarElement*>(interface_[i]->Discret().gElement(
  //            interface_[i]->Discret().ElementColMap()->GID(e)));
  //        mele->GetNitscheContainer().ClearAll();
  //      }
  //    }
}

void CONTACT::CoNitscheStrategyPoro::SetState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  if (statename == MORTAR::state_svelocity)
  {
    SetParentState(statename, vec);
  }
  else
    CONTACT::CoNitscheStrategy::SetState(statename, vec);
}

void CONTACT::CoNitscheStrategyPoro::SetParentState(
    const enum MORTAR::StateType& statename, const Epetra_Vector& vec)
{
  //
  if (statename == MORTAR::state_fvelocity || statename == MORTAR::state_fpressure)
  {
    Teuchos::RCP<DRT::Discretization> dis = DRT::Problem::Instance()->GetDis("porofluid");
    if (dis == Teuchos::null) dserror("didn't get my discretization");

    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(), true));
    CORE::LINALG::Export(vec, *global);


    // set state on interfaces
    for (const auto& interface : interface_)
    {
      DRT::Discretization& idiscret = interface->Discret();

      for (int j = 0; j < interface->Discret().ElementColMap()->NumMyElements(); ++j)
      {
        const int gid = interface->Discret().ElementColMap()->GID(j);

        auto* ele = dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(gid));

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;

        if (ele->ParentSlaveElement())  // if this pointer is NULL, this parent is impermeable
        {
          // this gets values in local order
          ele->ParentSlaveElement()->LocationVector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*global, myval, lm);

          std::vector<double> vel;
          std::vector<double> pres;

          for (int n = 0; n < ele->ParentSlaveElement()->NumNode(); ++n)
          {
            for (unsigned dim = 0; dim < 3; ++dim)
            {
              vel.push_back(myval[n * 4 + dim]);
            }
            pres.push_back(myval[n * 4 + 3]);
          }

          ele->MoData().ParentPFPres() = pres;
          ele->MoData().ParentPFVel() = vel;
          ele->MoData().ParentPFDof() = lm;
        }
      }
    }
  }
  else
    CONTACT::CoNitscheStrategy::SetParentState(statename, vec);
}

Teuchos::RCP<Epetra_FEVector> CONTACT::CoNitscheStrategyPoro::SetupRhsBlockVec(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  switch (bt)
  {
    case DRT::UTILS::VecBlockType::porofluid:
      return Teuchos::rcp(
          new Epetra_FEVector(*DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap()));
    default:
      return CONTACT::CoNitscheStrategy::SetupRhsBlockVec(bt);
  }
}

Teuchos::RCP<const Epetra_Vector> CONTACT::CoNitscheStrategyPoro::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::VecBlockType::porofluid:
      return Teuchos::rcp(new Epetra_Vector(Copy, *(fp_), 0));
    default:
      return CONTACT::CoNitscheStrategy::GetRhsBlockPtr(bp);
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::CoNitscheStrategyPoro::SetupMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_porofluid:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *DRT::Problem::Instance()->GetDis("structure")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    case DRT::UTILS::MatBlockType::porofluid_displ:
    case DRT::UTILS::MatBlockType::porofluid_porofluid:
      return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *Teuchos::rcpFromRef<const Epetra_Map>(
              *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap()),
          100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
    default:
      return CONTACT::CoNitscheStrategy::SetupMatrixBlockPtr(bt);
  }
}

void CONTACT::CoNitscheStrategyPoro::CompleteMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc)
{
  switch (bt)
  {
    case DRT::UTILS::MatBlockType::displ_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::porofluid_displ:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix())
              .GlobalAssemble(
                  *DRT::Problem::Instance()->GetDis("structure")->DofRowMap(),  // col map
                  *DRT::Problem::Instance()->GetDis("porofluid")->DofRowMap(),  // row map
                  true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    case DRT::UTILS::MatBlockType::porofluid_porofluid:
      if (dynamic_cast<Epetra_FECrsMatrix&>(*kc->EpetraMatrix()).GlobalAssemble(true, Add))
        dserror("GlobalAssemble(...) failed");
      break;
    default:
      CONTACT::CoNitscheStrategy::CompleteMatrixBlockPtr(bt, kc);
      break;
  }
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::CoNitscheStrategyPoro::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bp) const
{
  if (!curr_state_eval_) dserror("you didn't evaluate this contact state first");

  switch (bp)
  {
    case DRT::UTILS::MatBlockType::porofluid_porofluid:
      return kpp_;
    case DRT::UTILS::MatBlockType::porofluid_displ:
      return kpd_;
    case DRT::UTILS::MatBlockType::displ_porofluid:
      return kdp_;
    default:
      return CONTACT::CoNitscheStrategy::GetMatrixBlockPtr(bp, nullptr);
  }
}
