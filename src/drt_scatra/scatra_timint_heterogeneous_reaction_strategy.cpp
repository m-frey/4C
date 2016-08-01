/*----------------------------------------------------------------------*/
/*!
 \file scatra_timint_heterogeneous_reaction_strategy.cpp

 \brief Solution strategy for heterogeneous reactions. This is not meshtying!!!

 <pre>
    \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#include "scatra_timint_heterogeneous_reaction_strategy.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_lib/drt_dofset_merged_proxy.H"

#include "../drt_particle/binning_strategy.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"


/*----------------------------------------------------------------------*
 | constructor                                               vuong 06/16 |
 *----------------------------------------------------------------------*/
SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy(
    SCATRA::ScaTraTimIntImpl* scatratimint
    ):
    MeshtyingStrategyStd(scatratimint)
{
  return;
} // SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy


/*------------------------------------------------------------------------*
 | evaluate heterogeneous reactions (actually no mesh tying    vuong 06/16 |
 *------------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::EvaluateMeshtying()
{

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::calc_heteroreac_mat_and_rhs);

  // provide element parameter list with numbers of dofsets associated with displacement and velocity dofs on scatra discretization
  condparams.set<int>("ndsdisp",scatratimint_->NdsDisp());
  condparams.set<int>("ndsvel",scatratimint_->NdsVel());

  // set global state vectors according to time-integration scheme
  discret_->ClearState();
  discret_->SetState("phinp",scatratimint_->Phiafnp());
  discret_->SetState("hist",scatratimint_->Hist());

  // provide scatra discretization with convective velocity
  discret_->SetState(
      scatratimint_->NdsVel(),
      "convective velocity field",
      scatratimint_->Discretization()->GetState(scatratimint_->NdsVel(),"convective velocity field"));

  // provide scatra discretization with velocity
  discret_->SetState(
      scatratimint_->NdsVel(),
      "velocity field",
      scatratimint_->Discretization()->GetState(scatratimint_->NdsVel(),"velocity field"));

  if(scatratimint_->IsALE())
  {
    discret_->SetState(
        scatratimint_->NdsDisp(),
        "dispnp",
        scatratimint_->Discretization()->GetState(scatratimint_->NdsDisp(),"dispnp"));
  }

  discret_->Evaluate(condparams,scatratimint_->SystemMatrix(),scatratimint_->Residual());
  return;
} // SCATRA::HeterogeneousReactionStrategy::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | initialize meshtying objects                              vuong 06/16 |
 *----------------------------------------------------------------------*/
void SCATRA::HeterogeneousReactionStrategy::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  InitConvCheckStrategy();

  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( scatratimint_->Discretization()->Comm().Clone());

  // standard case
  discret_ = Teuchos::rcp(new DRT::Discretization(scatratimint_->Discretization()->Name(), com));

  // call complete without assigning degrees of freedom
  discret_->FillComplete(false,true,true);

  Teuchos::RCP<DRT::Discretization> scatradis = scatratimint_->Discretization();

  // create scatra elements if the scatra discretization is empty
  {
    // fill scatra discretization by cloning fluid discretization
    DRT::UTILS::CloneDiscretizationFromCondition<SCATRA::ScatraReactionCloneStrategy>(
        *scatradis,
        *discret_,
        "ScatraHeteroReactionSlave");

    // set implementation type of cloned scatra elements
    for(int i=0; i<discret_->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(discret_->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");

      element->SetImplType(INPAR::SCATRA::impltype_advreac);
    }
  }

  // redistribute discr. with help of binning strategy
  if(scatradis->Comm().NumProc()>1)
  {
    // create vector of discr.
    std::vector<Teuchos::RCP<DRT::Discretization> > dis;
    dis.push_back(scatradis);
    dis.push_back(discret_);

    std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

    /// binning strategy is created and parallel redistribution is performed
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
      Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
  }

  {
    // build a dofset that merges the DOFs from both sides
    Teuchos::RCP<DRT::DofSet> newdofset =
        Teuchos::rcp(new DRT::DofSetMergedProxy(
            scatradis->GetDofSetSubProxy(),
            scatradis,
            "ScatraHeteroReactionMaster",
            "ScatraHeteroReactionSlave"));

    // assign the dofset to the reaction discretization
    discret_->ReplaceDofSet(newdofset,false);

    // add all secondary dofsets as sub proxies
    for(int ndofset=1;ndofset<scatratimint_->Discretization()->NumDofSets();++ndofset)
      discret_->AddDofSet(scatratimint_->Discretization()->GetDofSetSubProxy(ndofset));

    // done. Rebuild all maps.
    discret_->FillComplete(true,true,true);
  }

  return;
}






