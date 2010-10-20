/*----------------------------------------------------------------------*/
/*!
\file thrtimint_impl_contact.cpp
\brief Thermal contact routines for Implicit time integration for
       spatial discretised thermal dynamics

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                                mgit 10/10 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                    mgit 10/10 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_impl.H"
#include "thr_aux.H"

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"

/*----------------------------------------------------------------------*
 |  Modify thermal system of equation towards thermal contact mgit 09/10|
 *----------------------------------------------------------------------*/
void THR::TimIntImpl::ApplyThermoContact(Teuchos::RCP<LINALG::SparseMatrix>& tang,
                         Teuchos::RCP<Epetra_Vector>& feff,
                         Teuchos::RCP<Epetra_Vector>& temp)
{
  // only in the case of contact
  if(cmtman_==Teuchos::null)
    return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  tang->Complete();

  // convert maps (from structure discretization to thermo discretization)
  // slave-, active-, inactive-, master-, activemaster-, n- smdofs
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  ConvertMaps (sdofs,adofs,mdofs);

  // map of active and master dofs
  amdofs = LINALG::MergeMap(adofs,mdofs,false);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);

  // row map of thermal problem
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(discret_->DofRowMap())));

  // split problemrowmap in n+am
  ndofs = LINALG::SplitMap(*problemrowmap,*smdofs);

  // modifications only for active nodes
  if (adofs->NumGlobalElements()==0)
    return;

  // assemble Mortar Matrices D and M in thermo dofs for active nodes
  RCP<LINALG::SparseMatrix> dmatrix = rcp(new LINALG::SparseMatrix(*sdofs,10));
  RCP<LINALG::SparseMatrix> mmatrix = rcp(new LINALG::SparseMatrix(*sdofs,100));

  AssembleDM(*dmatrix,*mmatrix);
  
  // FillComplete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*mdofs,*sdofs);

  // assemble Matrix A
  RCP<LINALG::SparseMatrix> amatrix = rcp(new LINALG::SparseMatrix(*sdofs,10));
  AssembleA(*amatrix);
  
  // Fill Complete
  amatrix->Complete();
  
  // active part of dmatrix and mmatrix
  RCP<Epetra_Map> tmp;
  RCP<LINALG::SparseMatrix> dmatrixa,mmatrixa,amatrixa,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
  LINALG::SplitMatrix2x2(dmatrix,adofs,idofs,adofs,idofs,dmatrixa,tmp1,tmp2,tmp3);
  LINALG::SplitMatrix2x2(mmatrix,adofs,idofs,mdofs,tmp,mmatrixa,tmp4,tmp5,tmp6);
  LINALG::SplitMatrix2x2(amatrix,adofs,idofs,sdofs,tmp,amatrixa,tmp4,tmp5,tmp6);

  // assemble mechanical dissipation
  RCP<Epetra_Vector> mechdissrate = LINALG::CreateVector(*mdofs,true);
  AssembleMechDissRate(*mechdissrate);

  // matrices from linearized thermal contact condition
  RCP<LINALG::SparseMatrix>  thermcontLM = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<LINALG::SparseMatrix>  thermcontTEMP = rcp(new LINALG::SparseMatrix(*adofs,3));
  RCP<Epetra_Vector>         thermcontRHS = LINALG::CreateVector(*adofs,true);

  // assemble thermal contact contition
  AssembleThermContCondition(*thermcontLM,*thermcontTEMP,*thermcontRHS,*dmatrixa,*mmatrixa,*amatrixa,adofs,mdofs);

  // complete the matrices
  thermcontLM->Complete(*sdofs,*adofs);
  thermcontTEMP->Complete(*smdofs,*adofs);

  /**********************************************************************/
  /* Modification of the stiff matrix and rhs towards thermo contact    */
  /**********************************************************************/

  /**********************************************************************/
  /* Create inv(D)                                                      */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*sdofs,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication M^ = inv(D) * M
  RCP<LINALG::SparseMatrix> mhatmatrix;
  mhatmatrix = LINALG::MLMultiply(*invd,false,*mmatrix,false,false,false,true);

  /**********************************************************************/
  /* Split tang into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  RCP<LINALG::SparseMatrix> tangmatrix = rcp(new LINALG::SparseMatrix(*tang));
  LINALG::SplitMatrix2x2(tangmatrix,smdofs,ndofs,smdofs,ndofs,ksmsm,ksmn,knsm,knn);

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,sdofs,mdofs,sdofs,mdofs,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,sdofs,mdofs,ndofs,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,ndofs,tempmap,sdofs,mdofs,kns,knm,tempmtx1,tempmtx2);

  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*problemrowmap,*feff,smdofs,fsm,ndofs,fn);
  LINALG::SplitVector(*smdofs,*fsm,sdofs,fs,mdofs,fm);

  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,adofs,idofs,adofs,idofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,adofs,idofs,ndofs,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,adofs,idofs,mdofs,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,mdofs,tempmap,adofs,idofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*adofs));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*idofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*sdofs,*fs,adofs,fa,idofs,fi);

  // abbreviations for active and inactive set
  int aset = adofs->NumGlobalElements();
  int iset = idofs->NumGlobalElements();

  // active part of invd and mhatmatrix
  RCP<Epetra_Map> tmpmap;
  RCP<LINALG::SparseMatrix> invda,mhata;
  LINALG::SplitMatrix2x2(invd,sdofs,tmpmap,adofs,idofs,invda,tmp1,tmp2,tmp3);
  LINALG::SplitMatrix2x2(mhatmatrix,adofs,idofs,mdofs,tmpmap,mhata,tmp1,tmp2,tmp3);

  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmnmod->Add(*kmn,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhata,true,*kan,false,false,false,true);
  kmnmod->Add(*kmnadd,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = rcp(new LINALG::SparseMatrix(*mdofs,100));
  kmmmod->Add(*kmm,false,1.0,1.0);
  RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhata,true,*kam,false,false,false,true);
  kmmmod->Add(*kmmadd,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmimod->Add(*kmi,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhata,true,*kai,false,false,false,true);
    kmimod->Add(*kmiadd,false,1.0,1.0);
    kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  }

  // kmi: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmamod->Add(*kma,false,1.0,1.0);
    RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhata,true,*kaa,false,false,false,true);
    kmamod->Add(*kmaadd,false,1.0,1.0);
    kmamod->Complete(kma->DomainMap(),kma->RowMap());
  }

  // kan: thermcontlm*invd*kan
  RCP<LINALG::SparseMatrix> kanmod;
  if (aset)
  {
    kanmod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
    kanmod->Complete(kan->DomainMap(),kan->RowMap());
  }

  // kam: thermcontlm*invd*kam
  RCP<LINALG::SparseMatrix> kammod;
  if (aset)
  {
    kammod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
    kammod->Complete(kam->DomainMap(),kam->RowMap());
  }

  // kai: thermcontlm*invd*kai
  RCP<LINALG::SparseMatrix> kaimod;
  if (aset && iset)
  {
    kaimod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
    kaimod->Complete(kai->DomainMap(),kai->RowMap());
  }

  // kaa: thermcontlm*invd*kaa
  RCP<LINALG::SparseMatrix> kaamod;
  if (aset)
  {
    kaamod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
    kaamod->Complete(kaa->DomainMap(),kaa->RowMap());
  }

  // Modifications towards rhs
  // FIXGIT: pay attention to genalpha
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*mdofs));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  // fa: mutliply with thermcontlm
  RCP<Epetra_Vector> famod;
  {
    famod = rcp(new Epetra_Vector(*adofs));
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
    temp->Multiply(false,*fa,*famod);
  }

  /**********************************************************************/
  /* Global setup of tangnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> tangnew = rcp(new LINALG::SparseMatrix(*problemrowmap,81,true,false,tangmatrix->GetMatrixtype()));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap);

  // add n submatrices to tangnew
  tangnew->Add(*knn,false,1.0,1.0);
  tangnew->Add(*knm,false,1.0,1.0);
  tangnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to tangnew
  tangnew->Add(*kmnmod,false,1.0,1.0);
  tangnew->Add(*kmmmod,false,1.0,1.0);
  if (iset) tangnew->Add(*kmimod,false,1.0,1.0);
  if (aset) tangnew->Add(*kmamod,false,1.0,1.0);

  // add i submatrices to tangnew
  if (iset) tangnew->Add(*kin,false,1.0,1.0);
  if (iset) tangnew->Add(*kim,false,1.0,1.0);
  if (iset) tangnew->Add(*kii,false,1.0,1.0);
  if (iset) tangnew->Add(*kia,false,1.0,1.0);

  // add a submatrices to tangnew
  if (aset) tangnew->Add(*kanmod,false,1.0,1.0);
  if (aset) tangnew->Add(*kammod,false,1.0,1.0);
  if (aset && iset) tangnew->Add(*kaimod,false,1.0,1.0);
  if (aset) tangnew->Add(*kaamod,false,1.0,1.0);

  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  // add mechanical dissipation to feffnew
  RCP<Epetra_Vector> mechdissrateexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*mechdissrate,*mechdissrateexp);
  feffnew->Update(-1.0,*mechdissrateexp,1.0);

  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp;
  if (iset)
  {
    fiexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*fi,*fiexp);
    feffnew->Update(1.0,*fiexp,1.0);
  }

  // add a subvector to feffnew
  RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*famod,*famodexp);
    feffnew->Update(1.0,*famodexp,+1.0);
  }

  // add linearized thermo contact condition
  tangnew->Add(*thermcontTEMP,false,-1.0,+1.0);

  // add rhs of thermal contact condition to feffnew
  RCP<Epetra_Vector> thermcontRHSexp = rcp(new Epetra_Vector(*problemrowmap));
  LINALG::Export(*thermcontRHS,*thermcontRHSexp);
  feffnew->Update(-1.0,*thermcontRHSexp,1.0);

  // FillComplete tangnew (square)
  tangnew->Complete();

  /**********************************************************************/
  /* Replace tang and feff by tangnew and feffnew                     */
  /**********************************************************************/
  tang = tangnew;
  feff = feffnew;

  // leave this place
  return;
}

/*----------------------------------------------------------------------*
 | convert maps form structure dofs to thermo dofs            mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::ConvertMaps(RCP<Epetra_Map>& slavedofs,
                         RCP<Epetra_Map>& activedofs,
                         RCP<Epetra_Map>& masterdofs)
{

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes/dofs
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // define local variables
    int slavecountnodes = 0;
    vector<int> myslavegids(slavenodes->NumMyElements());

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      myslavegids[slavecountnodes] = (discretstruct_->Dof(1,node))[0];
      ++slavecountnodes;
    }

    // resize the temporary vectors
    myslavegids.resize(slavecountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gslavecountnodes;
    Comm().SumAll(&slavecountnodes,&gslavecountnodes,1);

    // create active node map and active dof map
    slavedofs = rcp(new Epetra_Map(gslavecountnodes,slavecountnodes,&myslavegids[0],0,Comm()));

    // active nodes/dofs
    const RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // define local variables
    int countnodes = 0;
    vector<int> mynodegids(activenodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mynodegids[countnodes] = (discretstruct_->Dof(1,node))[0];
      ++countnodes;
    }

    // resize the temporary vectors
    mynodegids.resize(countnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gcountnodes;
    Comm().SumAll(&countnodes,&gcountnodes,1);

    // create active node map and active dof map
    activedofs = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));

    // master nodes/dofs
    const RCP<Epetra_Map> masternodes = interface[m]->MasterRowNodes();

    // define local variables
    int mastercountnodes = 0;
    vector<int> mymastergids(masternodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mymastergids[mastercountnodes] = (discretstruct_->Dof(1,node))[0];
      ++mastercountnodes;
    }

    // resize the temporary vectors
    mymastergids.resize(mastercountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gmastercountnodes;
    Comm().SumAll(&mastercountnodes,&gmastercountnodes,1);

    // create active node map and active dof map
    masterdofs = rcp(new Epetra_Map(gmastercountnodes,mastercountnodes,&mymastergids[0],0,Comm()));
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mortar matrices in thermo dofs                    mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleDM(LINALG::SparseMatrix& dmatrix,
                         LINALG::SparseMatrix& mmatrix)

{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // This is a little bit complicated and a lot of parallel stuff has to
  // be done here. The point is that, when assembling the mortar matrix
  // M, we need the temperature dof from the master node which can lie on
  // a complete different proc. For this reason, we have tho keep all procs
  // around

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveFullNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      /************************************************** D-matrix ******/
      if (Comm().MyPID()==cnode->Owner())
      {
        if ((cnode->MoData().GetD()).size()>0)
        {
          vector<map<int,double> > dmap = cnode->MoData().GetD();
          int rowdisp = cnode->Dofs()[0];
          double val = (dmap[0])[rowdisp];
          dmatrix.Assemble(val, rowtemp, rowtemp);
        }
      }

      /************************************************** M-matrix ******/
      set<int> mnodes;
      int mastergid=0;
      set<int>::iterator mcurr;
      int mastersize = 0;
      vector<map<int,double> > mmap;

      if (Comm().MyPID()==cnode->Owner())
      {
        mmap = cnode->MoData().GetM();
        mnodes = cnode->FriData().GetMNodes();
        mastersize = mnodes.size();
        mcurr = mnodes.begin();
      }

      // commiunicate number of master nodes
      Comm().Broadcast(&mastersize,1,cnode->Owner());

      // loop over all according master nodes
      for (int l=0;l<mastersize;++l)
      {
        if (Comm().MyPID()==cnode->Owner())
          mastergid=*mcurr;

        // communicate GID of masternode
        Comm().Broadcast(&mastergid,1,cnode->Owner());

        DRT::Node* mnode = (interface[m]->Discret()).gNode(mastergid);
        DRT::Node* mnodeges = discretstruct_->gNode(mastergid);

        // temperature and displacement dofs
        int coltemp = 0;
        int coldis = 0;
        if(Comm().MyPID()==mnode->Owner())
        {
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          coltemp = discretstruct_->Dof(1,mnodeges)[0];
          coldis = (cmnode->Dofs())[0];
        }

        // communicate temperature and displacement dof
        Comm().Broadcast(&coltemp,1,mnode->Owner());
        Comm().Broadcast(&coldis,1,mnode->Owner());

        // do the assembly
        if (Comm().MyPID()==cnode->Owner())
        {
          double val = mmap[0][coldis];
          if (abs(val)>1e-12) mmatrix.Assemble(val, rowtemp, coltemp);
          ++mcurr;
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble A matrix in thermo dofs                           mgit 10/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleA(LINALG::SparseMatrix& amatrix)

{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // This is a little bit complicated and a lot of parallel stuff has to
  // be done here. The point is that, when assembling the mortar matrix
  // M, we need the temperature dof from the master node which can lie on
  // a complete different proc. For this reason, we have tho keep all procs
  // around

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveFullNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      /************************************************** A-matrix ******/
      set<int> anodes;
      int slavegid=0;
      set<int>::iterator scurr;
      int slavesize = 0;
      vector<map<int,double> > amap;

      if (Comm().MyPID()==cnode->Owner())
      {
        amap = cnode->FriData().GetA();
        anodes = cnode->FriData().GetANodes();
        slavesize = anodes.size();
        scurr = anodes.begin();
      }

      // commiunicate number of master nodes
      Comm().Broadcast(&slavesize,1,cnode->Owner());

      // loop over all according master nodes
      for (int l=0;l<slavesize;++l)
      {
        if (Comm().MyPID()==cnode->Owner())
          slavegid=*scurr;

        // communicate GID of masternode
        Comm().Broadcast(&slavegid,1,cnode->Owner());

        DRT::Node* mnode = (interface[m]->Discret()).gNode(slavegid);
        DRT::Node* mnodeges = discretstruct_->gNode(slavegid);

        // temperature and displacement dofs
        int coltemp = 0;
        int coldis = 0;
        if(Comm().MyPID()==mnode->Owner())
        {
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          coltemp = discretstruct_->Dof(1,mnodeges)[0];
          coldis = (cmnode->Dofs())[0];
        }

        // communicate temperature and displacement dof
        Comm().Broadcast(&coltemp,1,mnode->Owner());
        Comm().Broadcast(&coldis,1,mnode->Owner());

        // do the assembly
        if (Comm().MyPID()==cnode->Owner())
        {
          double val = amap[0][coldis];
          if (abs(val)>1e-12) amatrix.Assemble(val, rowtemp, coltemp);
          ++scurr;
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mechanical dissipation for master nodes            mgit 08/10|
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleMechDissRate(Epetra_Vector& mechdissrate)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");
  
  // time step size
  double dt = GetTimeStepSize();

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // loop over master full nodes
    // master nodes are redundant on all procs and the entry of the
    // mechanical dissipation lies on proc which did the evaluation of 
    // mortar integrals and mechanical dissipation 

    // master nodes
    const RCP<Epetra_Map> masternodes = interface[m]->MasterFullNodes();

    // loop over all masternodes nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // mechanical dissipation to be assembled     
      double mechdissglobal = 0;
      
      // mechanical dissipation on proc
      double mechdissproc = 1/dt*cnode->MechDiss();
            
      // sum all entries to mechdissglobal
      Comm().SumAll(&mechdissproc,&mechdissglobal,1);
      
      // check if entry is only from one processor
      if(mechdissproc!=mechdissglobal)
      {  
        if (abs(mechdissproc)>1e-12)
          dserror ("Error in AssembleMechDissRate: Entries from more than one proc");
      }
      
      // owner of master node does the assembly
      if(Comm().MyPID()==cnode->Owner())
      {
        // row dof of temperature
        int rowtemp = discretstruct_->Dof(1,nodeges)[0];

        Epetra_SerialDenseVector mechdissiprate(1);
        vector<int> dof(1);
        vector<int> owner(1);

        mechdissiprate(0) = mechdissglobal;
        dof[0] = rowtemp;
        owner[0] = cnode->Owner();

        // do assembly
        if(abs(mechdissiprate(0))>1e-12)
          LINALG::Assemble(mechdissrate, mechdissiprate, dof, owner);
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble the thermal contact conditions for slave nodes    mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::TimIntImpl::AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                                    LINALG::SparseMatrix& thermcontTEMP,
                                                    Epetra_Vector& thermcontRHS,
                                                    LINALG::SparseMatrix& dmatrix,
                                                    LINALG::SparseMatrix& mmatrix,
                                                    LINALG::SparseMatrix& amatrix,
                                                    RCP<Epetra_Map> activedofs,
                                                    RCP<Epetra_Map> masterdofs)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet and for one heat
  // transfer coefficient
  // FIXGIT: The heat transfer coefficient should be a condition on
  // the single interfaces!!
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::AssembleThermContCondition: Only for one interface yet.");

  // heat transfer coefficient for slave and master surface
  double heattranss = interface[0]->IParams().get<double>("HEATTRANSSLAVE");
  double heattransm = interface[0]->IParams().get<double>("HEATTRANSMASTER");

  if (heattranss <= 0 or heattransm <= 0)
   dserror("Error: Choose realistic heat transfer parameter");

  // time step size
  double dt = GetTimeStepSize();

  double beta = heattranss*heattransm/(heattranss+heattransm);
  double delta = heattranss/(heattranss+heattransm);

  // with respect to Lagrange multipliers
  thermcontLM.Add(amatrix,false,1.0,1.0);

  // with respect to temperature
  thermcontTEMP.Add(dmatrix,false,-beta,1.0);
  thermcontTEMP.Add(mmatrix,false,+beta,1.0);

  RCP<Epetra_Vector> fa, fm, rest1, rest2;

  // row map of thermal problem
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(discret_->DofRowMap())));

  LINALG::SplitVector(*problemrowmap,*tempn_,activedofs,fa,masterdofs,fm);

  RCP <Epetra_Vector> DdotTemp = rcp(new Epetra_Vector(*activedofs));
  dmatrix.Multiply(false,*fa,*DdotTemp);
  thermcontRHS.Update(-beta,*DdotTemp,1.0);

  RCP <Epetra_Vector> MdotTemp = rcp(new Epetra_Vector(*activedofs));
  mmatrix.Multiply(false,*fm,*MdotTemp);
  thermcontRHS.Update(+beta,*MdotTemp,1.0);

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes (full map)
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(Comm().MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      Epetra_SerialDenseVector mechdissiprate(1);
      vector<int> dof(1);
      vector<int> owner(1);

      mechdissiprate(0) = delta/dt*cnode->MechDiss();
      dof[0] = rowtemp;
      owner[0] = cnode->Owner();

      // do assembly
      if(abs(mechdissiprate(0)>1e-12 and cnode->Active()))
        LINALG::Assemble(thermcontRHS, mechdissiprate, dof, owner);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
