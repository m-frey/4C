/*!----------------------------------------------------------------------
\file drt_contact_manager_base.cpp
\brief Main class to control all contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
<\pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialComm.h"
#include "drt_contact_manager_base.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::ManagerBase::ManagerBase() :
dim_(0),
alphaf_(0.0),
activesetconv_(false),
activesetsteps_(0),
isincontact_(false)
{
  //**********************************************************************
  // empty constructor
  //**********************************************************************
  // Setup of the contact library has to be done by a derived class. This
  // derived class is specific to the FEM code into which the contact
  // library is meant to be integrated. For BACI this is realized via the
  // CONTACT::Manager class! There the following actions are performed:
  //**********************************************************************
  // 1) get problem dimension (2D or 3D) and store into dim_
  // 2) read and check contact input parameters
  // 3) read and check contact boundary conditions
  // 4) build contact interfaces
  //**********************************************************************
  
  // create a simple serial communicator
  comm_ = rcp(new Epetra_SerialComm());
  
  return;
}

/*----------------------------------------------------------------------*
 |  set current deformation state (public)                    popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::SetState(const string& statename,
                                    const RCP<Epetra_Vector> vec)
{
  // set state on interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->SetState(statename,vec);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  initialize Mortar stuff for next Newton step (public)     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::InitializeMortar()
{
  // initialize / reset interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Initialize();
  }

  // intitialize Dold and Mold if not done already
   if (dold_==null)
   {
     dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
     dold_->Zero();
     dold_->Complete();
   }
   if (mold_==null)
   {
     mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
     mold_->Zero();
     mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
   }
   
   // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
   dmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
   mmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
   mhatmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
   g_          = LINALG::CreateVector(*gsnoderowmap_,true);
   
   // (re)setup global matrices containing fc derivatives
   lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
   linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));

  return;
}

/*----------------------------------------------------------------------*
 |  initialize contact for next Newton step (public)          popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::Initialize()
{
  // (re)setup global normal and tangent matrices
  nmatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  tmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  
  // (re)setup global Tresca friction matrix L and vector R
  string ftype   = scontact_.get<string>("friction type","none");
  if (ftype=="tresca")
  {
    lmatrix_ = rcp(new LINALG::SparseMatrix(*gslipt_,10));
    r_       = LINALG::CreateVector(*gslipt_,true);
  }
  
  // (re)setup global matrices containing derivatives
  smatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  pmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate Mortar matrices D,M & gap g~ only (public)       popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::EvaluateMortar()
{
  /**********************************************************************/
  /* evaluate interfaces                                                */
  /* (nodal normals, projections, Mortar integration, Mortar assembly)  */
  /**********************************************************************/
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }
  
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::Evaluate(RCP<LINALG::SparseMatrix> kteff,
                                    RCP<Epetra_Vector> feff)
{ 
  // check if Tresca friction and/or basis transformation should be applied
  string ftype   = scontact_.get<string>("friction type","none");
  bool btrafo = scontact_.get<bool>("basis transformation",false);
  
  // Tresca friction cases
  if (ftype=="tresca")
  {
    if (btrafo)
      EvaluateTrescaBasisTrafo(kteff,feff);
    else
      EvaluateTrescaNoBasisTrafo(kteff,feff);
  }
  
  // Other cases (Frictionless, Stick, MeshTying)
  else
  {
  if (btrafo)
    EvaluateBasisTrafo(kteff,feff);
  else
    EvaluateNoBasisTrafo(kteff,feff);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Tresca friction with basis trafo (public)        mgit 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::EvaluateTrescaBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                                    RCP<Epetra_Vector> feff)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes,         */
  /* global matrix t with tangent vectors of active nodes               */
  /* and global matrix l and vector r for frictional contact            */
  /**********************************************************************/
  
  // read tresca friction bound
  double frbound = scontact_.get<double>("friction bound",0.0);
  
  // read weighting factor ct
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both firctional conditions are violated. Here we have to weigh
  // the two violations via ct!
  double ct = scontact_.get<double>("semismooth ct",0.0);

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleTresca(*lmatrix_,*r_,frbound,ct); 
  }
  
  // FillComplete() global matrices N, T and L
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
  if(gslipt_->NumGlobalElements()) lmatrix_->Complete(*gslipt_,*gslipt_);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
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
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!) 
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
  // we also need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;
  
  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
  
  /**********************************************************************/
  /* Apply basis transformation to k                                    */
  /**********************************************************************/
  // define temporary RCP mod
  RCP<LINALG::SparseMatrix> mod;
  
  // kss: nothing to do
  RCP<LINALG::SparseMatrix> kssmod = kss;
  
  // ksm: add kss*T(mbar)
  RCP<LINALG::SparseMatrix> ksmmod = LINALG::Multiply(*kss,false,*mhatmatrix_,false,false);
  ksmmod->Add(*ksm,false,1.0,1.0);
  ksmmod->Complete(ksm->DomainMap(),ksm->RowMap());
  
  // ksn: nothing to do
  RCP<LINALG::SparseMatrix> ksnmod = ksn;
  
  // kms: add T(mbar)*kss
  RCP<LINALG::SparseMatrix> kmsmod = LINALG::Multiply(*mhatmatrix_,true,*kss,false,false);
  kmsmod->Add(*kms,false,1.0,1.0);
  kmsmod->Complete(kms->DomainMap(),kms->RowMap());
  
  // kmm: add kms*T(mbar) + T(mbar)*ksm + T(mbar)*kss*mbar
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*kms,false,*mhatmatrix_,false,false);
  mod = LINALG::Multiply(*mhatmatrix_,true,*ksm,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  mod = LINALG::Multiply(*mhatmatrix_,true,*kss,false);
  mod = LINALG::Multiply(*mod,false,*mhatmatrix_,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  kmmmod->Add(*kmm,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmn: add T(mbar)*ksn
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhatmatrix_,true,*ksn,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kns: nothing to do
  RCP<LINALG::SparseMatrix> knsmod = kns;
  
  // knm: add kns*mbar
  RCP<LINALG::SparseMatrix> knmmod = LINALG::Multiply(*kns,false,*mhatmatrix_,false,false);
  knmmod->Add(*knm,false,1.0,1.0);
  knmmod->Complete(knm->DomainMap(),knm->RowMap());
  
  // knn: nothing to do
  RCP<LINALG::SparseMatrix> knnmod = knn;
  
  /**********************************************************************/
  /* Apply basis transformation to f                                    */
  /**********************************************************************/
  // fs: nothing to be done
  RCP<Epetra_Vector> fsmod = fs;
  
  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhatmatrix_->Multiply(true,*fs,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
  
  // fn: nothing to be done
  RCP<Epetra_Vector> fnmod = fn;
  
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaamod, kaimod, kiamod, kiimod;
  
  // we want to split ksnmod / ksmmod into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kanmod, kinmod, kammod, kimmod;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kssmod,gactivedofs_,gidofs,gactivedofs_,gidofs,kaamod,kaimod,kiamod,kiimod);
  LINALG::SplitMatrix2x2(ksnmod,gactivedofs_,gidofs,gndofrowmap_,tempmap,kanmod,tempmtx1,kinmod,tempmtx2);
  LINALG::SplitMatrix2x2(ksmmod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kammod,tempmtx1,kimmod,tempmtx2);
  

  /**********************************************************************/
  /* Split active quantities into slip / stick                          */
  /**********************************************************************/

  // we want to split kaamod into 2 groups sl,st = 4 blocks
  RCP<LINALG::SparseMatrix> kslslmod, kslstmod, kstslmod, kststmod;
  
  // we want to split kanmod / kammod / kaimod into 2 groups sl,st = 2 blocks
  RCP<LINALG::SparseMatrix> kslnmod, kstnmod, kslmmod, kstmmod, kslimod, kstimod; 
 
  // some temporary RCPs
  RCP<Epetra_Map> temp1map;
  RCP<LINALG::SparseMatrix> temp1mtx1;
  RCP<LINALG::SparseMatrix> temp1mtx2;
  RCP<LINALG::SparseMatrix> temp1mtx3;
  
  // we will get the stick rowmap as a by-product
  RCP<Epetra_Map> gstdofs;
   
  // do the splitting
  LINALG::SplitMatrix2x2(kaamod,gslipdofs_,gstdofs,gslipdofs_,gstdofs,kslslmod,kslstmod,kstslmod,kststmod);
  LINALG::SplitMatrix2x2(kanmod,gslipdofs_,gstdofs,gndofrowmap_,temp1map,kslnmod,temp1mtx1,kstnmod,temp1mtx2);
  LINALG::SplitMatrix2x2(kammod,gslipdofs_,gstdofs,gmdofrowmap_,temp1map,kslmmod,temp1mtx1,kstmmod,temp1mtx2);
  LINALG::SplitMatrix2x2(kaimod,gslipdofs_,gstdofs,gidofs,temp1map,kslimod,temp1mtx1,kstimod,temp1mtx2);
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> famod, fimod;
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    famod = rcp(new Epetra_Vector(*fsmod));
  else if (!gactivedofs_->NumGlobalElements())
    fimod = rcp(new Epetra_Vector(*fsmod));
  else
  {
    LINALG::SplitVector(*fsmod,*gactivedofs_,famod,*gidofs,fimod);
  }
  
  /**********************************************************************/
  /* Isolate active and slip part from invd and dold                              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invda, invdsl;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invd_,gslipdofs_,gstdofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));
  invdsl->Scale(1/(1-alphaf_));
  
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
    
  /**********************************************************************/
  /* Gen-alpha modifications                                            */
  /**********************************************************************/
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fimod->Update(-alphaf_,*tempveci,1.0);
  }
   
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    famod->Update(-alphaf_,*tempvec,1.0);
  }
  
  // we want to split famod into 2 groups sl,st
  RCP<Epetra_Vector> fslmod, fstmod;
  
  // do the vector splitting a -> sl+st
  if(gactivedofs_->NumGlobalElements())
  { 
    if (!gstdofs->NumGlobalElements())
      fslmod = rcp(new Epetra_Vector(*famod));
    else if (!gslipdofs_->NumGlobalElements())
      fstmod = rcp(new Epetra_Vector(*famod));
    else
    {
      LINALG::SplitVector(*famod,*gslipdofs_,fslmod,*gstdofs,fstmod);
    }
  }
  // we will get the stickt rowmap as a by-product
  RCP<Epetra_Map> gstickt;
  
  // temporary RCPs
  RCP<Epetra_Map> tmap;
  
  // temporary RCPs
  RCP<LINALG::SparseMatrix> tm1, tm2;
  
  // we want to split the tmatrix_ into 2 groups
  RCP<LINALG::SparseMatrix> tslmatrix, tstmatrix;
  LINALG::SplitMatrix2x2(tmatrix_,gslipt_,gstickt,gslipdofs_,tmap,tslmatrix,tm1,tm2,tstmatrix); 
  
  // do the multiplications with t matrix
  RCP<LINALG::SparseMatrix> tkslnmod, tkslmmod, tkslimod, tkslslmod, tlmatrix, tkslstmod;
  RCP<Epetra_Vector> tfslmod;
  
  if(gslipdofs_->NumGlobalElements())
  {
    // kslnmod: multiply with tslmatrix
    tkslnmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
    tkslnmod = LINALG::Multiply(*tkslnmod,false,*kslnmod,false,true);
      
    // kslmmod: multiply with tslmatrix
    tkslmmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
    tkslmmod = LINALG::Multiply(*tkslmmod,false,*kslmmod,false,true);
      
    // friction
    // lmatrix: multiply with tslmatrix 
    tlmatrix = LINALG::Multiply(*lmatrix_,false,*tslmatrix,false,true);
    
    // kslslmod: multiply with tslmatrix
    tkslslmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
    tkslslmod = LINALG::Multiply(*tkslslmod,false,*kslslmod,false,false);
      
    // add tlmatirx to tkslslmod
    tkslslmod->Add(*tlmatrix,false,1.0,1.0);
    tkslslmod->Complete(kslslmod->DomainMap(),kslslmod->RowMap());
    
    if (gidofs->NumGlobalElements())
    {
      //kslimod: multiply with tslmatrix
      tkslimod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
      tkslimod = LINALG::Multiply(*tkslimod,false,*kslimod,false,true);
    }       
      
    if (gstdofs->NumGlobalElements())
    {
      //kslstmod: multiply with tslmatrix
      tkslstmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
      tkslstmod = LINALG::Multiply(*tkslstmod,false,*kslstmod,false,true);
    }       
          
    // fslmod: multiply with tmatrix    
    tfslmod = rcp(new Epetra_Vector(*gslipt_));
    RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
    temp->Multiply(false,*fslmod,*tfslmod); 
    
    // friction
    // add r to famod
    tfslmod->Update(1.0,*r_,1.0);
  }
  
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);
  
  // add n / m submatrices to kteffnew
  kteffnew->Add(*knnmod,false,1.0,1.0);
  kteffnew->Add(*knmmod,false,1.0,1.0);
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  
  // add a / i submatrices to kteffnew, if existing
  if (knsmod!=null) kteffnew->Add(*knsmod,false,1.0,1.0);
  if (kmsmod!=null) kteffnew->Add(*kmsmod,false,1.0,1.0);
  if (kinmod!=null) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (kimmod!=null) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (kiimod!=null) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (kiamod!=null) kteffnew->Add(*kiamod,false,1.0,1.0);
  
  // add matrix of normals to kteffnew
  kteffnew->Add(*nmatrix_,false,1.0,1.0);
  // add matrix of tangents of sticky nodes to kteffnew
  if(tstmatrix!=null) kteffnew->Add(*tstmatrix,false,1.0,1.0);
  
  // add submatrices with tangents to kteffnew, if existing
  if (tkslnmod!=null) kteffnew->Add(*tkslnmod,false,1.0,1.0);
  if (tkslmmod!=null) kteffnew->Add(*tkslmmod,false,1.0,1.0);
  if (tkslimod!=null) kteffnew->Add(*tkslimod,false,1.0,1.0);
  if (tkslslmod!=null) kteffnew->Add(*tkslslmod,false,1.0,1.0);
  if (tkslstmod!=null) kteffnew->Add(*tkslstmod,false,1.0,1.0);
    
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n / m subvectors to feffnew
  RCP<Epetra_Vector> fnmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fnmod,*fnmodexp);
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fnmodexp,1.0,*fmmodexp,1.0);
  
  // add i / ta subvectors to feffnew, if existing
  RCP<Epetra_Vector> fimodexp = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fstmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> tfslmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  
  if (fimod!=null) LINALG::Export(*fimod,*fimodexp);
  if (tfslmod!=null) LINALG::Export(*tfslmod,*tfslmodexp);
  
    
  feffnew->Update(1.0,*fimodexp,1.0);
  feffnew->Update(1.0,*tfslmodexp,1.0);
  
  // add weighted gap vector to feffnew, if existing
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  if (gact->GlobalLength()) LINALG::Export(*gact,*gexp);
  feffnew->Update(1.0,*gexp,1.0);
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                   */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  //LINALG::PrintSparsityToPostscript(*(kteff->EpetraMatrix()));
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate trecsa friction without basis trafo (public) gitterle 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::EvaluateTrescaNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                                      RCP<Epetra_Vector> feff)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  
  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
    
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix l and vector r for frictional contact            */
  /**********************************************************************/
  // here and for the splitting later, we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  // read tresca friction bound
  double frbound = scontact_.get<double>("friction bound",0.0);
    
  // read weighting factor ct
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both firctional conditions are violated. Here we have to weigh
  // the two violations via ct!
  double ct = scontact_.get<double>("semismooth ct",0.0);

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
//    interface_[i]->AssembleS(*smatrix_);
//    interface_[i]->AssembleP(*pmatrix_);
//    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
    interface_[i]->AssembleTresca(*lmatrix_,*r_,frbound,ct);
  }
    
  // FillComplete() global matrices N and T and L
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
   
  // FillComplete() global matrix  L
  lmatrix_->Complete(*gslipt_,*gslipt_);  
  
  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofs,*gactiven_);
  
  // FillComplete() global matrix P
  // (actually gsdofrowmap_ is in general sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  pmatrix_->Complete(*gsmdofs,*gactivet_);
  
  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofs,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofs,*gmdofrowmap_);
  
  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
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
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
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
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
    
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;
  
  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);
  
  /**********************************************************************/
  /* Split active quantities into slip / stick                          */
  /**********************************************************************/

  // we want to split kaa into 2 groups sl,st = 4 blocks
  RCP<LINALG::SparseMatrix> kslsl, kslst, kstsl, kstst;
    
  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  RCP<LINALG::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti; 
  
  // some temporary RCPs
  RCP<Epetra_Map> temp1map;
  RCP<LINALG::SparseMatrix> temp1mtx1;
  RCP<LINALG::SparseMatrix> temp1mtx2;
  RCP<LINALG::SparseMatrix> temp1mtx3;
  
  // we will get the stick rowmap as a by-product
  RCP<Epetra_Map> gstdofs;
     
  // do the splitting
  LINALG::SplitMatrix2x2(kaa,gslipdofs_,gstdofs,gslipdofs_,gstdofs,kslsl,kslst,kstsl,kstst);
  LINALG::SplitMatrix2x2(kan,gslipdofs_,gstdofs,gndofrowmap_,temp1map,ksln,temp1mtx1,kstn,temp1mtx2);
  LINALG::SplitMatrix2x2(kam,gslipdofs_,gstdofs,gmdofrowmap_,temp1map,kslm,temp1mtx1,kstm,temp1mtx2);
  LINALG::SplitMatrix2x2(kai,gslipdofs_,gstdofs,gidofs,temp1map,ksli,temp1mtx1,ksti,temp1mtx2);
  
  // we want to split fs into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*gactivedofs_));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*gidofs));
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    *fa = *fs;
  else if (!gactivedofs_->NumGlobalElements())
    *fi = *fs;
  else
  {
    LINALG::SplitVector(*fs,*gactivedofs_,fa,*gidofs,fi);
  }
  
  /**********************************************************************/
  /* Isolate active and slip part from mhat, invd and dold              */
  /* Isolate slip part from T                                           */
  /**********************************************************************/
  
  RCP<LINALG::SparseMatrix> mhata, mhatsl, mhatst;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(mhatmatrix_,gslipdofs_,gstdofs,gmdofrowmap_,tempmap,mhatsl,tempmtx2,mhatst,tempmtx3);
  mhata_=mhata;
  
  RCP<LINALG::SparseMatrix> invda, invdsl;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invd_,gslipdofs_,gstdofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));
  invdsl->Scale(1/(1-alphaf_));
  
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
  
  // we will get the stickt rowmap as a by-product
  RCP<Epetra_Map> gstickt;
    
  // temporary RCPs
  RCP<Epetra_Map> tmap;
    
  // temporary RCPs
  RCP<LINALG::SparseMatrix> tm1, tm2;
    
  // we want to split the tmatrix_ into 2 groups
  RCP<LINALG::SparseMatrix> tslmatrix, tstmatrix;
  LINALG::SplitMatrix2x2(tmatrix_,gslipt_,gstickt,gslipdofs_,tmap,tslmatrix,tm1,tm2,tstmatrix);
  
    
  /**********************************************************************/
  /* Split LinD and LinM into blocks                                    */
  /**********************************************************************/
  // we want to split lindmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> lindai, lindaa, lindam, lindas;
  
  // we want to split linmmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> linmmi, linmma, linmmm, linmms;
  
  if (fulllin)
  {
    // do the splitting
    LINALG::SplitMatrix2x2(lindmatrix_,gactivedofs_,gidofs,gmdofrowmap_,gsdofrowmap_,lindam,lindas,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(lindas,gactivedofs_,tempmap,gactivedofs_,gidofs,lindaa,lindai,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmmatrix_,gmdofrowmap_,tempmap,gmdofrowmap_,gsdofrowmap_,linmmm,linmms,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,linmma,linmmi,tempmtx1,tempmtx2);
  
    // modification of kai, kaa, kam
    // (this has to be done first as they are needed below)
    // (note, that kai, kaa, kam have to be UNcompleted again first!!!)
    kai->UnComplete();
    kaa->UnComplete();
    kam->UnComplete();
    kai->Add(*lindai,false,1.0-alphaf_,1.0);
    kaa->Add(*lindaa,false,1.0-alphaf_,1.0);
    kam->Add(*lindam,false,1.0-alphaf_,1.0);
    kai->Complete(*gidofs,*gactivedofs_);
    kaa->Complete();
    kam->Complete(*gmdofrowmap_,*gactivedofs_);
  }
  
  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do
  
  // knm: nothing to do
  
  // kns: nothing to do
  
  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhata,true,*kan,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhata,true,*kam,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  if (fulllin) kmmmod->Add(*linmmm,false,1.0-alphaf_,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod = LINALG::Multiply(*mhata,true,*kai,false,false);
  kmimod->Add(*kmi,false,1.0,1.0);
  if (fulllin) kmimod->Add(*linmmi,false,1.0-alphaf_,1.0);
  kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  
  // kma: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod = LINALG::Multiply(*mhata,true,*kaa,false,false);
  kmamod->Add(*kma,false,1.0,1.0);
  if (fulllin) kmamod->Add(*linmma,false,1.0-alphaf_,1.0);
  kmamod->Complete(kma->DomainMap(),kma->RowMap());
  
  // kin: nothing to do
  
  // kim: nothing to do
  
  // kii: nothing to do
  
  // kisl: nothing to do
  
  // kist: nothing to do

  // n*mbaractive: do the multiplication 
  RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata,false,true);
  
  // t*mbarstick: do the multiplication
  RCP<LINALG::SparseMatrix> tmhatst = LINALG::Multiply(*tstmatrix,false,*mhatst,false,true);
  
  // nmatrix: nothing to do
  
  // ksln: multiply with tslmatrix
  RCP<LINALG::SparseMatrix> kslnmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  kslnmod = LINALG::Multiply(*kslnmod,false,*ksln,false,true);
  
  // kslm: multiply with tslmatrix
  RCP<LINALG::SparseMatrix> kslmmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  kslmmod = LINALG::Multiply(*kslmmod,false,*kslm,false,false);
  
  // friction 
  // lmatrix: multiply with tslmatrix, also multiply with mbarslip 
  RCP<LINALG::SparseMatrix> ltmatrix = LINALG::Multiply(*lmatrix_,false,*tslmatrix,false,true);
  RCP<LINALG::SparseMatrix> ltmatrixmb = LINALG::Multiply(*ltmatrix,false,*mhatsl,false,true);
  
   // subtract ltmatrixmb from kslmmod
  kslmmod->Add(*ltmatrixmb,false,-1.0,1.0);
  kslmmod->Complete(kslm->DomainMap(),kslm->RowMap());
  
  // ksli: multiply with tslmatrix
  RCP<LINALG::SparseMatrix> kslimod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  kslimod = LINALG::Multiply(*kslimod,false,*ksli,false,true);
  
  // kslsl: multiply with tslmatrix
  RCP<LINALG::SparseMatrix> kslslmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  kslslmod = LINALG::Multiply(*kslslmod,false,*kslsl,false,true);
  
  // add tlmatirx to kslslmod
  kslslmod->Add(*ltmatrix,false,1.0,1.0);
  kslslmod->Complete(kslsl->DomainMap(),kslsl->RowMap());
  
  // slstmod: multiply with tslmatrix
  RCP<LINALG::SparseMatrix> kslstmod = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  kslstmod = LINALG::Multiply(*kslstmod,false,*kslst,false,true);
    
  // fn: nothing to do
  
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fi->Update(-alphaf_,*tempveci,1.0);
  }
  
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    fa->Update(-alphaf_,*tempvec,1.0);
  }
 
  // we want to split famod into 2 groups sl,st
  RCP<Epetra_Vector> fsl, fst;
    
  // do the vector splitting a -> sl+st
  if(gactivedofs_->NumGlobalElements())
  { 
    if (!gstdofs->NumGlobalElements())
      fsl = rcp(new Epetra_Vector(*fa));
    else if (!gslipdofs_->NumGlobalElements())
      fst = rcp(new Epetra_Vector(*fa));
    else
    {
      LINALG::SplitVector(*fa,*gslipdofs_,fsl,*gstdofs,fst);
    }
  }
  
  // fm: add alphaf * old contact forces (t_n)
  RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mold_->Multiply(true,*zold_,*tempvecm);
  fm->Update(alphaf_,*tempvecm,1.0);
    
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
  
  // fsl: mutliply with tmatrix
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> fslmod = rcp(new Epetra_Vector(*gslipt_));
  RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*tslmatrix,false,*invdsl,false,true);
  
  if(gslipdofs_->NumGlobalElements())
  {
  temp->Multiply(false,*fsl,*fslmod);

  // friction
  // add r to fslmod
  fslmod->Update(1.0,*r_,1.0);
  }  
  // gactive: nothing to do
  
  // add jump from stick nodes to r.h.s.
  // mostly nonzero, not when changing a slip node to a stick one within
  // a time step in the semi-smooth Newton
  
  // temporary RCPs
  RCP<Epetra_Map> tmap1;
  
  // temporary vector
  RCP<Epetra_Vector> restjump;
  
  RCP<Epetra_Vector> stjump = rcp(new Epetra_Vector(*gstdofs));
  if (gstdofs->NumGlobalElements()) LINALG::Export(*jump_,*stjump);
  
  RCP<Epetra_Vector> tstjump = rcp(new Epetra_Vector(*gstickt));
  
  tstmatrix->Multiply(false,*stjump,*tstjump);
      
 #ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> deriv = rcp(new LINALG::SparseMatrix(*gactiven_,81));
    deriv->Add(*nmatrix_,false,1.0,1.0);
    deriv->Add(*smatrix_,false,1.0,1.0);
    deriv->Add(*nmhata,false,-1.0,1.0);
    deriv->Complete(*gsmdofs,*gactiven_);
    cout << *deriv << endl;
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP
  
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i=0; i<(int)interface_.size();++i)
  {
    cout << *pmatrix_ << endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM
    
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);
  
  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmimod,false,1.0,1.0);
  kteffnew->Add(*kmamod,false,1.0,1.0);
  
  // add i submatrices to kteffnew
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kin,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kim,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kii,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kia,false,1.0,1.0);
  
  // add matrix nmhata to keteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmhata,false,-1.0,1.0);

  // add matrix n to kteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmatrix_,false,1.0,1.0);

  // add matrix t to kteffnew
  if(tstmatrix!=null) kteffnew->Add(*tstmatrix,false,1.0,1.0);
  
  // add matrix tmhatst to kteffnew
  if(tmhatst!=null) kteffnew->Add(*tmhatst,false,-1.0,1.0);

  // add full linearization terms to kteffnew
  if (fulllin)
  {
   if (gactiven_->NumGlobalElements()) kteffnew->Add(*smatrix_,false,1.0,1.0);
   if (gactivet_->NumGlobalElements()) kteffnew->Add(*pmatrix_,false,-1.0,1.0);
  }
  
  // add a submatrices to kteffnew
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslnmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslmmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslimod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslslmod,false,1.0,1.0);
  if (gslipt_->NumGlobalElements()) kteffnew->Add(*kslstmod,false,1.0,1.0);
    
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);
  
  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);
  
  // add i and sl subvector to feffnew
  RCP<Epetra_Vector> fiexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fi,*fiexp);
  if (gidofs->NumGlobalElements()) feffnew->Update(1.0,*fiexp,1.0);
  
  // add a subvector to feffnew
  RCP<Epetra_Vector> fslmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fslmod,*fslmodexp);
  if (gslipnodes_->NumGlobalElements())feffnew->Update(1.0,*fslmodexp,1.0);
  
  //stick nodes: add tstjump to r.h.s
  if(gstdofs->NumGlobalElements())
  { 
  RCP<Epetra_Vector> tstjumpexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*tstjump,*tstjumpexp);
  feffnew->Update(-1.0,*tstjumpexp,+1.0);
  }
  
  // add weighted gap vector to feffnew, if existing
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*gact,*gexp);
  if (gact->GlobalLength()) feffnew->Update(1.0,*gexp,1.0);
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact with basis transformation (public)       popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::EvaluateBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                              RCP<Epetra_Vector> feff)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /**********************************************************************/
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    
  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
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
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!) 
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
  // we also need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  
  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
    
  /**********************************************************************/
  /* Apply basis transformation to k                                    */
  /**********************************************************************/
  // define temporary RCP mod
  RCP<LINALG::SparseMatrix> mod;
  
  // kss: nothing to do
  RCP<LINALG::SparseMatrix> kssmod = kss;
  
  // ksm: add kss*T(mbar)
  RCP<LINALG::SparseMatrix> ksmmod = LINALG::Multiply(*kss,false,*mhatmatrix_,false,false);
  ksmmod->Add(*ksm,false,1.0,1.0);
  ksmmod->Complete(ksm->DomainMap(),ksm->RowMap());
  
  // ksn: nothing to do
  RCP<LINALG::SparseMatrix> ksnmod = ksn;
  
  // kms: add T(mbar)*kss
  RCP<LINALG::SparseMatrix> kmsmod = LINALG::Multiply(*mhatmatrix_,true,*kss,false,false);
  kmsmod->Add(*kms,false,1.0,1.0);
  kmsmod->Complete(kms->DomainMap(),kms->RowMap());
  
  // kmm: add kms*T(mbar) + T(mbar)*ksm + T(mbar)*kss*mbar
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*kms,false,*mhatmatrix_,false,false);
  mod = LINALG::Multiply(*mhatmatrix_,true,*ksm,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  mod = LINALG::Multiply(*mhatmatrix_,true,*kss,false);
  mod = LINALG::Multiply(*mod,false,*mhatmatrix_,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  kmmmod->Add(*kmm,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmn: add T(mbar)*ksn
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhatmatrix_,true,*ksn,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kns: nothing to do
  RCP<LINALG::SparseMatrix> knsmod = kns;
  
  // knm: add kns*mbar
  RCP<LINALG::SparseMatrix> knmmod = LINALG::Multiply(*kns,false,*mhatmatrix_,false,false);
  knmmod->Add(*knm,false,1.0,1.0);
  knmmod->Complete(knm->DomainMap(),knm->RowMap());
  
  // knn: nothing to do
  RCP<LINALG::SparseMatrix> knnmod = knn;
  
  /**********************************************************************/
  /* Apply basis transformation to f                                    */
  /**********************************************************************/
  // fs: nothing to be done
  RCP<Epetra_Vector> fsmod = fs;
  
  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhatmatrix_->Multiply(true,*fs,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
  
  // fn: nothing to be done
  RCP<Epetra_Vector> fnmod = fn;
  
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaamod, kaimod, kiamod, kiimod;
  
  // we want to split ksnmod / ksmmod into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kanmod, kinmod, kammod, kimmod;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kssmod,gactivedofs_,gidofs,gactivedofs_,gidofs,kaamod,kaimod,kiamod,kiimod);
  LINALG::SplitMatrix2x2(ksnmod,gactivedofs_,gidofs,gndofrowmap_,tempmap,kanmod,tempmtx1,kinmod,tempmtx2);
  LINALG::SplitMatrix2x2(ksmmod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kammod,tempmtx1,kimmod,tempmtx2);
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> famod, fimod;
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    famod = rcp(new Epetra_Vector(*fsmod));
  else if (!gactivedofs_->NumGlobalElements())
    fimod = rcp(new Epetra_Vector(*fsmod));
  else
  {
    LINALG::SplitVector(*fsmod,*gactivedofs_,famod,*gidofs,fimod);
  }
  
  /**********************************************************************/
  /* Isolate active / inactive part from dold                           */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
  
  /**********************************************************************/
  /* Gen-alpha modifications                                            */
  /**********************************************************************/
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fimod->Update(-alphaf_,*tempveci,1.0);
  }
   
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    famod->Update(-alphaf_,*tempvec,1.0);
  }
   
  // do the multiplications with t matrix
  RCP<LINALG::SparseMatrix> tkanmod, tkammod, tkaimod, tkaamod;
  RCP<Epetra_Vector> tfamod;
  
  if(gactivedofs_->NumGlobalElements())
  {
    tkanmod = LINALG::Multiply(*tmatrix_,false,*kanmod,false);
    tkammod = LINALG::Multiply(*tmatrix_,false,*kammod,false);
    tkaamod = LINALG::Multiply(*tmatrix_,false,*kaamod,false);
    
    if (gidofs->NumGlobalElements())
      tkaimod = LINALG::Multiply(*tmatrix_,false,*kaimod,false);
    
    tfamod = rcp(new Epetra_Vector(tmatrix_->RowMap()));
    tmatrix_->Multiply(false,*famod,*tfamod);
  }
  
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);
  
  // add n / m submatrices to kteffnew
  kteffnew->Add(*knnmod,false,1.0,1.0);
  kteffnew->Add(*knmmod,false,1.0,1.0);
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  
  // add a / i submatrices to kteffnew, if existing
  if (knsmod!=null) kteffnew->Add(*knsmod,false,1.0,1.0);
  if (kmsmod!=null) kteffnew->Add(*kmsmod,false,1.0,1.0);
  if (kinmod!=null) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (kimmod!=null) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (kiimod!=null) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (kiamod!=null) kteffnew->Add(*kiamod,false,1.0,1.0);
  
  // add matrix of normals to kteffnew
  kteffnew->Add(*nmatrix_,false,1.0,1.0);

  if (ftype=="none")
  {
    // add submatrices with tangents to kteffnew, if existing
    if (tkanmod!=null) kteffnew->Add(*tkanmod,false,1.0,1.0);
    if (tkammod!=null) kteffnew->Add(*tkammod,false,1.0,1.0);
    if (tkaimod!=null) kteffnew->Add(*tkaimod,false,1.0,1.0);
    if (tkaamod!=null) kteffnew->Add(*tkaamod,false,1.0,1.0);
  }
  else if (ftype=="stick")
  {
    // add matrix of tangents to kteffnew
    if (tmatrix_!=null) kteffnew->Add(*tmatrix_,false,1.0,1.0);
  }
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n / m subvectors to feffnew
  RCP<Epetra_Vector> fnmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fnmod,*fnmodexp);
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fnmodexp,1.0,*fmmodexp,1.0);
  
  // add i / ta subvectors to feffnew, if existing
  RCP<Epetra_Vector> fimodexp = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> tfamodexp = rcp(new Epetra_Vector(*problemrowmap_));
  if (fimod!=null) LINALG::Export(*fimod,*fimodexp);
  if (tfamod!=null) LINALG::Export(*tfamod,*tfamodexp);
  
  if (ftype=="none")
   feffnew->Update(1.0,*fimodexp,1.0,*tfamodexp,1.0);
  else if (ftype=="stick")
    feffnew->Update(1.0,*fimodexp,0.0,*tfamodexp,1.0);
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  if (ctype!="meshtying")
  {
    // add weighted gap vector to feffnew, if existing
    RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
    if (gact->GlobalLength()) LINALG::Export(*gact,*gexp);
    feffnew->Update(1.0,*gexp,1.0);
  }
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                   */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  //LINALG::PrintSparsityToPostscript(*(kteff->EpetraMatrix()));
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact without basis transformation (public)    popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::EvaluateNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                                RCP<Epetra_Vector> feff)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  
  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /**********************************************************************/
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /* and global matrix s with normal derivatives of active nodes        */
  /**********************************************************************/
  // here and for the splitting later, we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleP(*pmatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
  }
    
  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
  
  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofs,*gactiven_);
  
  // FillComplete() global matrix P
  // (actually gsdofrowmap_ is in general sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  pmatrix_->Complete(*gsmdofs,*gactivet_);
  
  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofs,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofs,*gmdofrowmap_);
  
  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
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
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
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
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
    
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;
  
  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*gactivedofs_));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*gidofs));
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    *fa = *fs;
  else if (!gactivedofs_->NumGlobalElements())
    *fi = *fs;
  else
  {
    LINALG::SplitVector(*fs,*gactivedofs_,fa,*gidofs,fi);
  }
  
  /**********************************************************************/
  /* Isolate active part from mhat, invd and dold                       */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> mhata;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);
  mhata_=mhata;
  
  RCP<LINALG::SparseMatrix> invda;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));
  
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
  
  /**********************************************************************/
  /* Split LinD and LinM into blocks                                    */
  /**********************************************************************/
  // we want to split lindmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> lindai, lindaa, lindam, lindas;
  
  // we want to split linmmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> linmmi, linmma, linmmm, linmms;
  
  if (fulllin)
  {
    // do the splitting
    LINALG::SplitMatrix2x2(lindmatrix_,gactivedofs_,gidofs,gmdofrowmap_,gsdofrowmap_,lindam,lindas,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(lindas,gactivedofs_,tempmap,gactivedofs_,gidofs,lindaa,lindai,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmmatrix_,gmdofrowmap_,tempmap,gmdofrowmap_,gsdofrowmap_,linmmm,linmms,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,linmma,linmmi,tempmtx1,tempmtx2);
  
    // modification of kai, kaa, kam
    // (this has to be done first as they are needed below)
    // (note, that kai, kaa, kam have to be UNcompleted again first!!!)
    kai->UnComplete();
    kaa->UnComplete();
    kam->UnComplete();
    kai->Add(*lindai,false,1.0-alphaf_,1.0);
    kaa->Add(*lindaa,false,1.0-alphaf_,1.0);
    kam->Add(*lindam,false,1.0-alphaf_,1.0);
    kai->Complete(*gidofs,*gactivedofs_);
    kaa->Complete();
    kam->Complete(*gmdofrowmap_,*gactivedofs_);
  }
  
  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do
  
  // knm: nothing to do
  
  // kns: nothing to do
  
  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhata,true,*kan,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhata,true,*kam,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  if (fulllin) kmmmod->Add(*linmmm,false,1.0-alphaf_,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod = LINALG::Multiply(*mhata,true,*kai,false,false);
  kmimod->Add(*kmi,false,1.0,1.0);
  if (fulllin) kmimod->Add(*linmmi,false,1.0-alphaf_,1.0);
  kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  
  // kma: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod = LINALG::Multiply(*mhata,true,*kaa,false,false);
  kmamod->Add(*kma,false,1.0,1.0);
  if (fulllin) kmamod->Add(*linmma,false,1.0-alphaf_,1.0);
  kmamod->Complete(kma->DomainMap(),kma->RowMap());
  
  // kin: nothing to do
  
  // kim: nothing to do
  
  // kii: nothing to do
  
  // kia: nothing to do
  
  // n*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata,false,true);
  
  // nmatrix: nothing to do
  
  // kan: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kanmod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kanmod = LINALG::Multiply(*kanmod,false,*kan,false,true);
  
  // kam: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kammod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kammod = LINALG::Multiply(*kammod,false,*kam,false,true);
  
  // kai: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaimod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaimod = LINALG::Multiply(*kaimod,false,*kai,false,true);
  
  // kaa: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaamod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaamod = LINALG::Multiply(*kaamod,false,*kaa,false,true);
  
  // t*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> tmhata = LINALG::Multiply(*tmatrix_,false,*mhata,false,true);
  
  // fn: nothing to do
  
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fi->Update(-alphaf_,*tempveci,1.0);
  }
  
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    fa->Update(-alphaf_,*tempvec,1.0);
  }
    
  // fm: add alphaf * old contact forces (t_n)
  RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mold_->Multiply(true,*zold_,*tempvecm);
  fm->Update(alphaf_,*tempvecm,1.0);
    
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
    
  // fa: mutliply with tmatrix
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> famod = rcp(new Epetra_Vector(*gactivet_));
  RCP<LINALG::SparseMatrix> tinvda = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  tinvda->Multiply(false,*fa,*famod);
  
  // gactive: nothing to do
  
#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> deriv = rcp(new LINALG::SparseMatrix(*gactiven_,81));
    deriv->Add(*nmatrix_,false,1.0,1.0);
    deriv->Add(*smatrix_,false,1.0,1.0);
    deriv->Add(*nmhata,false,-1.0,1.0);
    deriv->Complete(*gsmdofs,*gactiven_);
    cout << *deriv << endl;
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP
  
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i=0; i<(int)interface_.size();++i)
  {
    cout << *pmatrix_ << endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM
    
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*problemrowmap_,81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap_);
  
  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmimod,false,1.0,1.0);
  kteffnew->Add(*kmamod,false,1.0,1.0);
  
  // add i submatrices to kteffnew
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kin,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kim,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kii,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kia,false,1.0,1.0);
  
  // add matrix nmhata to keteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmhata,false,-1.0,1.0);

  // add matrix n to kteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmatrix_,false,1.0,1.0);
  
  // add full linearization terms to kteffnew
  if (fulllin)
  {
   if (gactiven_->NumGlobalElements()) kteffnew->Add(*smatrix_,false,1.0,1.0);
   if (gactivet_->NumGlobalElements()) kteffnew->Add(*pmatrix_,false,-1.0,1.0);
  }
  
  if (ftype=="none")
  {
    // add a submatrices to kteffnew
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kanmod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kammod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaimod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaamod,false,1.0,1.0);
  }
  else if (ftype=="stick")
  {
    // add matrices t and tmhata to kteffnew
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*tmatrix_,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*tmhata,false,-1.0,1.0);
  }
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);
  
  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);
  
  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fi,*fiexp);
  if (gidofs->NumGlobalElements()) feffnew->Update(1.0,*fiexp,1.0);
  
  if (ctype!="meshtying")
  {
    // add weighted gap vector to feffnew, if existing
    RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*gact,*gexp);
    if (gact->GlobalLength()) feffnew->Update(1.0,*gexp,1.0);
  }
  
  if (ftype=="none")
  {
    // add a subvector to feffnew
    RCP<Epetra_Vector> famodexp = rcp(new Epetra_Vector(*problemrowmap_));
    LINALG::Export(*famod,*famodexp);
    if (gactivenodes_->NumGlobalElements())feffnew->Update(1.0,*famodexp,1.0);
  }
  else if (ftype=="stick")
  {
    // do nothing here
  }
  else
    dserror("ERROR: Invalid type of friction law");  
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;

  return;
}

/*----------------------------------------------------------------------*
 |  Recovery method for displacements and LM (public)         popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::Recover(RCP<Epetra_Vector> disi)
{ 
  // check if basis transformation should be applied
  bool btrafo = scontact_.get<bool>("basis transformation",false);
  
  if (btrafo)
    RecoverBasisTrafo(disi);
  else
    RecoverNoBasisTrafo(disi);
    
  return;
}

/*----------------------------------------------------------------------*
 |  Recovery method (basis trafo case)                        popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::RecoverBasisTrafo(RCP<Epetra_Vector> disi)
{ 
  // extract incremental jump from disi (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*incrjump_);
  
  // friction
  // sum up incremental jumps from active set nodes
  jump_->Update(1.0,*incrjump_,1.0);
  
  // friction
  // store updaded jumps to nodes
  StoreNodalQuantities(ManagerBase::jump);
  
  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);
  
  // recover slave displacement increments
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(mhatmatrix_->RowMap()));
  mhatmatrix_->Multiply(false,*disim,*mod);
  
  RCP<Epetra_Vector> modexp = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*mod,*modexp);
  disi->Update(1.0,*modexp,1.0);
  
  /**********************************************************************/
  /* Update Lagrange multipliers                                        */
  /**********************************************************************/
  // approximate update
  //invd_->Multiply(false,*fs_,*z_);
  
  // full update
  z_->Update(1.0,*fs_,0.0);
  RCP<Epetra_Vector> mod2 = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> slavedisp = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*slavedisp);
  kss_->Multiply(false,*slavedisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  RCP<Epetra_Vector> masterdisp = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*masterdisp);
  ksm_->Multiply(false,*masterdisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  RCP<Epetra_Vector> innerdisp = rcp(new Epetra_Vector(*gndofrowmap_));
  LINALG::Export(*disi,*innerdisp);
  ksn_->Multiply(false,*innerdisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  dold_->Multiply(false,*zold_,*mod);
  z_->Update(-alphaf_,*mod,1.0);
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  z_->Scale(1/(1-alphaf_));
  
  // store updated LM into nodes
  StoreNodalQuantities(ManagerBase::lmupdate);
  
  /*
  // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
#ifdef DEBUG
  //debugging (check for z_i = 0)
  RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*z_,*zinactive);
    cout << *zinactive << endl;
  }
  
  // debugging (check for N*[d_a] = g_a and T*z_a = 0)
  if (gactivedofs_->NumGlobalElements())
  { 
    RCP<Epetra_Vector> activejump = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> gtest = rcp(new Epetra_Vector(*gactiven_));
    LINALG::Export(*incrjump_,*activejump);
    nmatrix_->Multiply(false,*activejump,*gtest);
    cout << *gtest << endl << *g_ << endl;
    
    RCP<Epetra_Vector> zactive = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> zerotest = rcp(new Epetra_Vector(*gactivet_));
    LINALG::Export(*z_,*zactive);
    tmatrix_->Multiply(false,*zactive,*zerotest);
    cout << *zerotest << endl;
    cout << *zactive << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 | Recovery method (no basis trafo case) (public)             popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::RecoverNoBasisTrafo(RCP<Epetra_Vector> disi)
{
  // extract slave displacements from disi
  RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*disis);
  
  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);
  
  // recover incremental jump (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  mhatmatrix_->Multiply(false,*disim,*incrjump_);
  incrjump_->Update(1.0,*disis,-1.0);
  
  // friction
  // sum up incremental jumps from active set nodes
  jump_->Update(1.0,*incrjump_,1.0);
  // friction
  // store updaded jumps to nodes
  StoreNodalQuantities(ManagerBase::jump);
    
  /**********************************************************************/
  /* Update Lagrange multipliers z_n+1                                  */
  /**********************************************************************/
  // approximate update
  //invd_->Multiply(false,*fs_,*z_);
    
  // full update
  z_->Update(1.0,*fs_,0.0);
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gsdofrowmap_));
  kss_->Multiply(false,*disis,*mod);
  z_->Update(-1.0,*mod,1.0);
  ksm_->Multiply(false,*disim,*mod);
  z_->Update(-1.0,*mod,1.0);
  RCP<Epetra_Vector> disin = rcp(new Epetra_Vector(*gndofrowmap_));
  LINALG::Export(*disi,*disin);
  ksn_->Multiply(false,*disin,*mod);
  z_->Update(-1.0,*mod,1.0);
  dold_->Multiply(false,*zold_,*mod);
  z_->Update(-alphaf_,*mod,1.0);
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  z_->Scale(1/(1-alphaf_));
  
  // store updated LM into nodes
  StoreNodalQuantities(ManagerBase::lmupdate);
    
  /* 
  // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
#ifdef DEBUG
  //debugging (check for z_i = 0)
  RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*z_,*zinactive);
    cout << *zinactive << endl;
  }
  
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  
  // debugging (check for N*[d_a] = g_a and T*z_a = 0)
  if (gactivedofs_->NumGlobalElements())
  { 
    RCP<Epetra_Vector> activejump = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> gtest = rcp(new Epetra_Vector(*gactiven_));
    RCP<Epetra_Vector> gtest2 = rcp(new Epetra_Vector(*gactiven_));
    LINALG::Export(*incrjump_,*activejump);
    nmatrix_->Multiply(false,*activejump,*gtest);
    
    RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
    RCP<Epetra_Vector> disism = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*disi,*disism);
    if (fulllin)
    {
      smatrix_->Multiply(false,*disism,*gtest2);
      gtest->Update(1.0,*gtest2,1.0);
    }
    cout << *gtest << endl << *g_ << endl;
    
    RCP<Epetra_Vector> zactive = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> zerotest = rcp(new Epetra_Vector(*gactivet_));
    RCP<Epetra_Vector> zerotest2 = rcp(new Epetra_Vector(*gactivet_));
    LINALG::Export(*z_,*zactive);
    tmatrix_->Multiply(false,*zactive,*zerotest);
    if (fulllin)
    {
    pmatrix_->Multiply(false,*disis,*zerotest2);
    zerotest->Update(1.0,*zerotest2,1.0);
    }
    cout << *zerotest << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::UpdateActiveSet()
{
  // get input parameter ctype
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  // assume that active set has converged and check for opposite
  activesetconv_=true;
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // friction
      double tz = 0.0;
      double tjump = 0.0;
      
      if(ftype=="tresca")
      { 
        // compute tangential part of Lagrange multiplier
        tz = cnode->txi()[0]*cnode->lm()[0] + cnode->txi()[1]*cnode->lm()[1];
        
        // compute tangential part of Lagrange multiplier
        tjump = cnode->txi()[0]*cnode->jump()[0] + cnode->txi()[1]*cnode->jump()[1];
      }
      
      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << endl;
        
        // check for penetration
        if (wgap < 0)
        {
          cnode->Active() = true;
          activesetconv_ = false;
        }
      }
      
      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << endl;
        
        // check for tensile contact forces
        if (nz <= 0) // no averaging of Lagrange multipliers
        //if (0.5*nz+0.5*nzold <= 0) // averaging of Lagrange multipliers
        {
          if (ctype!="meshtying")
          {
            cnode->Active() = false;
            
            // friction
            if(ftype=="tresca")
            {
              cnode->Slip() = false;    
            }
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        }

        // friction
        else
        {
          if(ftype=="tresca")
          {
            double frbound = scontact_.get<double>("friction bound",0.0);
            double ct = scontact_.get<double>("semismooth ct",0.0);
            
            if(cnode->Slip() == false)  
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if(abs(tz+ct*tjump)-frbound > 0) {}
                // do nothing (slip was correct)
              else
              {
               cnode->Slip() = false;
               activesetconv_ = false;
              }
            }
          } // if(ftype=="tresca")
        } // if (nz <= 0)
      } // if (cnode->Active()==false)
    } // loop over all slave nodes
  } // loop over all interfaces
  
  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);
  
  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }
  
  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=null) zigzagthree_  = rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=null) zigzagtwo_    = rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=null) zigzagone_ = rcp(new Epetra_Map(*gactivenodes_));

  // update zig-zagging history for slip nodes (shift by one)
  if (zigzagsliptwo_!=null) zigzagslipthree_  = rcp(new Epetra_Map(*zigzagsliptwo_));
  if (zigzagslipone_!=null) zigzagsliptwo_    = rcp(new Epetra_Map(*zigzagslipone_));
  if (gslipnodes_!=null) zigzagslipone_ = rcp(new Epetra_Map(*gslipnodes_));

  
  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;
  
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
  }
  
  // CHECK FOR ZIG-ZAGGING / JAMMING OF THE ACTIVE SET
  // *********************************************************************
  // A problem of the active set strategy which sometimes arises is known
  // from optimization literature as jamming or zig-zagging. This means
  // that within a load/time-step the algorithm can have more than one
  // solution due to the fact that the active set is not unique. Hence the
  // algorithm jumps between the solutions of the active set. The non-
  // uniquenesss results either from highly curved contact surfaces or
  // from the FE discretization, Thus the uniqueness of the closest-point-
  // projection cannot be guaranteed.
  // *********************************************************************
  // To overcome this problem we monitor the development of the active
  // set scheme in our contact algorithms. We can identify zig-zagging by
  // comparing the current active set with the active set of the second-
  // and third-last iteration. If an identity occurs, we consider the
  // active set strategy as converged instantly, accepting the current
  // version of the active set and proceeding with the next time/load step.
  // This very simple approach helps stabilizing the contact algorithm!
  // *********************************************************************
  bool zigzagging = false;

  if(ftype!="tresca") // FIXGIT: For tresca friction zig-zagging is not eliminated 
  {
    if (ActiveSetSteps()>2)
    {
      if (zigzagtwo_!=null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_) and zigzagsliptwo_->SameAs(*gslipnodes_)) 
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;
        
          // output to screen
          if (Comm().MyPID()==0)
            cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << endl;
        }
      }
      
      if (zigzagthree_!=null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_) and zigzagslipthree_->SameAs(*gslipnodes_))
        {
          // set active set converged
          activesetconv_ = true;
          zigzagging = true;
        
          // output to screen
          if (Comm().MyPID()==0)
            cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << endl;
        }
      }
    }
  } // if (ftype!="tresca")
  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = null;
    zigzagtwo_  = null;
    zigzagthree_= null;
  }
  
  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET ITERATION " << ActiveSetSteps()-1
         << " NOT CONVERGED - REPEAT TIME STEP................." << endl;
  else if (Comm().MyPID()==0 && activesetconv_==true)
    cout << "ACTIVE SET CONVERGED IN " << ActiveSetSteps()-zigzagging
         << " STEP(S)................." << endl;
    
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  
  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::UpdateActiveSetSemiSmooth()
{
  // FIXME: Here we do not consider zig-zagging yet!
  
  // get input parameter ctype
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weigh the two violations via cn!
  double cn = scontact_.get<double>("semismooth cn",0.0);
        
  // assume that active set has converged and check for opposite
  activesetconv_=true;
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // friction
      double tz;
      double tjump;
      
      if(ftype=="tresca")
      { 
        // compute tangential part of Lagrange multiplier
        tz = cnode->txi()[0]*cnode->lm()[0] + cnode->txi()[1]*cnode->lm()[1];
        
        // compute tangential part of Lagrange multiplier
        tjump = cnode->txi()[0]*cnode->jump()[0] + cnode->txi()[1]*cnode->jump()[1];
      }  

      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << endl;
                
        // check for penetration and/or tensile contact forces
        if (nz - cn*wgap > 0)
        {
          cnode->Active() = true;
//          cnode->Slip() = true;
          activesetconv_ = false;
        }
      }
      
      // check nodes of active set ***************************************
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << endl;
                  
        // check for tensile contact forces and/or penetration
        if (nz - cn*wgap <= 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
        {
          if (ctype!="meshtying")
          {
            cnode->Active() = false;
            
            // friction
            if(ftype=="tresca")
            {
              cnode->Slip() = false;    
            }
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        } 
        
        // friction
        else
        {
          if(ftype=="tresca")
          {
            double frbound = scontact_.get<double>("friction bound",0.0);
            double ct = scontact_.get<double>("semismooth ct",0.0);
            
            if(cnode->Slip() == false)  
            {
              // check (tz+ct*tjump)-frbound <= 0
              if(abs(tz+ct*tjump)-frbound <= 0) {}
                // do nothing (stick was correct)
              else
              {
                 cnode->Slip() = true;
                 activesetconv_ = false;
              }
            }
            else
            {
              // check (tz+ct*tjump)-frbound > 0
              if(abs(tz+ct*tjump)-frbound > 0) {}
               // do nothing (slip was correct)
              else
              {
               cnode->Slip() = false;
               activesetconv_ = false;
              }
            }
          } // if (fytpe=="tresca")
        } // if (nz - cn*wgap <= 0)
      } // if (cnode->Active()==false)
    } // loop over all slave nodes
  } // loop over all interfaces

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);
  
  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }
  
  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;
  
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
  }
  
  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET HAS CHANGED... CHANGE No. " << ActiveSetSteps()-1 << endl;
  
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  
  return;
}
/*----------------------------------------------------------------------*
 |  Compute contact forces (public)                           popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::ContactForces(RCP<Epetra_Vector> fresm)
{
  // Note that we ALWAYS use a TR-like approach to compute the contact
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}
  
  // FIXME: fresm is only here for debugging purposes!
  // compute two subvectors of fc each via Lagrange multipliers z_n+1, z_n
  RCP<Epetra_Vector> fcslavetemp  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  RCP<Epetra_Vector> fcslavetempend  = rcp(new Epetra_Vector(dold_->RowMap()));
  RCP<Epetra_Vector> fcmastertempend = rcp(new Epetra_Vector(mold_->DomainMap()));
  dmatrix_->Multiply(false,*z_,*fcslavetemp);
  mmatrix_->Multiply(true,*z_,*fcmastertemp);
  dold_->Multiply(false,*zold_,*fcslavetempend);
  mold_->Multiply(true,*zold_,*fcmastertempend);
  
  // export the contact forces to full dof layout
  RCP<Epetra_Vector> fcslave  = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcslaveend  = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmasterend = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcslavetemp,*fcslave);
  LINALG::Export(*fcmastertemp,*fcmaster);
  LINALG::Export(*fcslavetempend,*fcslaveend);
  LINALG::Export(*fcmastertempend,*fcmasterend);
  
  // build total contact force vector (TR-like!!!)
  fc_=fcslave;
  fc_->Update(-(1.0-alphaf_),*fcmaster,1.0-alphaf_);
  fc_->Update(alphaf_,*fcslaveend,1.0);
  fc_->Update(-alphaf_,*fcmasterend,1.0);
  
  /*
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
#ifdef DEBUG
  RCP<Epetra_Vector> fresmslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fresmmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  LINALG::Export(*fresm,*fresmslave);
  LINALG::Export(*fresm,*fresmmaster);
  
  vector<double> gfcs(3);
  vector<double> ggfcs(3);
  vector<double> gfcm(3);
  vector<double> ggfcm(3);
  int dimcheck = (gsdofrowmap_->NumGlobalElements())/(gsnoderowmap_->NumGlobalElements());
  if (dimcheck!=2 && dimcheck!=3) dserror("ERROR: ContactForces: Debugging for 3D not implemented yet");
  
  for (int i=0;i<fcslavetemp->MyLength();++i)
  {
    if ((i+dimcheck)%dimcheck == 0) gfcs[0]+=(*fcslavetemp)[i];
    else if ((i+dimcheck)%dimcheck == 1) gfcs[1]+=(*fcslavetemp)[i];
    else if ((i+dimcheck)%dimcheck == 2) gfcs[2]+=(*fcslavetemp)[i];
    else dserror("ERROR: Contact Forces: Dim. error in debugging part!");
  }
  
  for (int i=0;i<fcmastertemp->MyLength();++i)
  {
    if ((i+dimcheck)%dimcheck == 0) gfcm[0]-=(*fcmastertemp)[i];
    else if ((i+dimcheck)%dimcheck == 1) gfcm[1]-=(*fcmastertemp)[i];
    else if ((i+dimcheck)%dimcheck == 2) gfcm[2]-=(*fcmastertemp)[i];
    else dserror("ERROR: Contact Forces: Dim. error in debugging part!");
  }
  
  for (int i=0;i<3;++i)
  {
    Comm().SumAll(&gfcs[i],&ggfcs[i],1);
    Comm().SumAll(&gfcm[i],&ggfcm[i],1);
  }
  
  double slavenorm = 0.0;
  fcslavetemp->Norm2(&slavenorm);
  double slavenormend = 0.0;
  fcslavetempend->Norm2(&slavenormend);
  double fresmslavenorm = 0.0;
  fresmslave->Norm2(&fresmslavenorm);
  if (Comm().MyPID()==0)
  {
    cout << "Slave Contact Force Norm (n+1):  " << slavenorm << endl;
    cout << "Slave Contact Force Norm (n):  " << slavenormend << endl;
    cout << "Slave Residual Force Norm: " << fresmslavenorm << endl;
    cout << "Slave Contact Force Vector: " << ggfcs[0] << " " << ggfcs[1] << " " << ggfcs[2] << endl;
  }
  double masternorm = 0.0;
  fcmastertemp->Norm2(&masternorm);
  double masternormend = 0.0;
  fcmastertempend->Norm2(&masternormend);
  double fresmmasternorm = 0.0;
  fresmmaster->Norm2(&fresmmasternorm);
  if (Comm().MyPID()==0)
  {
    cout << "Master Contact Force Norm (n+1): " << masternorm << endl;
    cout << "Master Contact Force Norm (n): " << masternormend << endl;
    cout << "Master Residual Force Norm " << fresmmasternorm << endl;
    cout << "Master Contact Force Vector: " << ggfcm[0] << " " << ggfcm[1] << " " << ggfcm[2] << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::StoreNodalQuantities(ManagerBase::QuantityType type,
                                                RCP<Epetra_Vector> vec)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n interfaces!");
    
    // get global quantity to be stored in nodes
    RCP<Epetra_Vector> vectorglobal = null;
    switch(type)
    {
    case ManagerBase::lmcurrent:
    {
      vectorglobal = LagrMult();
      break;
    }
    case ManagerBase::lmold:
    {
      vectorglobal = LagrMultOld();
      break;
    }
    case ManagerBase::lmupdate:
    {
      vectorglobal = LagrMult();
      break;
    }
    case ManagerBase::jump:
    {
      vectorglobal = Jump(); 
      break;
    }
    case ManagerBase::dirichlet:
    {
      if (vec==null) dserror("Dirichtoggle vector has to be applied on input");
      vectorglobal = vec;
      break;
    }
    default:
      dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
    } // switch
    
    // export global quantity to current interface slave dof row map
    RCP<Epetra_Map> sdofrowmap = interface_[i]->SlaveRowDofs();
    RCP<Epetra_Vector> vectorinterface = rcp(new Epetra_Vector(*sdofrowmap));
    LINALG::Export(*vectorglobal,*vectorinterface);
    
    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // be aware of problem dimension
      int dim = Dim();
      
      // index for first DOF of current node in Epetra_Vector
      int locindex = vectorinterface->Map().LID(dim*gid);
      
      // extract this node's quantity from vectorinterface
      for (int k=0;k<dim;++k)
      {
        switch(type)
        {
        case ManagerBase::lmcurrent:
        {
          cnode->lm()[k] = (*vectorinterface)[locindex+k];
          break;
        }
        case ManagerBase::lmold:
        {
          cnode->lmold()[k] = (*vectorinterface)[locindex+k];
          break;
        }
        case ManagerBase::lmupdate:
        {
          // print a warning if a non-DBC inactive dof has a non-zero value
          // (only in semi-smooth Newton case, of course!)
          bool semismooth = scontact_.get<bool>("semismooth newton",false);
          if (semismooth && !cnode->Dbc()[k] && !cnode->Active() && abs((*vectorinterface)[locindex+k])>1.0e-8)
            cout << "***WARNING***: Non-D.B.C. inactive node " << cnode->Id() << " has non-zero Lag. Mult.: dof "
                 << cnode->Dofs()[k] << " lm " << (*vectorinterface)[locindex+k] << endl;
          
          // throw a dserror if node is Active and DBC
          if (cnode->Dbc()[k] && cnode->Active())
            dserror("ERROR: Slave Node %i is active and at the same time carries D.B.C.s!", cnode->Id());
          
          // explicity set global Lag. Mult. to zero for D.B.C nodes
          if (cnode->IsDbc())
            (*vectorinterface)[locindex+k] = 0.0;
          
          // store updated LM into node
          cnode->lm()[k] = (*vectorinterface)[locindex+k];
          break;
        }
        case ManagerBase::jump:
        {
          cnode->jump()[k] = (*vectorinterface)[locindex+k];
          break;
        }
        case ManagerBase::dirichlet:
        {
          cnode->Dbc()[k] = (*vectorinterface)[locindex+k];
          break;
        }
        default:
          dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
        } // switch 
      }
    }
  }
    
  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::StoreDM(const string& state)
{
  //store Dold and Mold matrix in D and M
  if (state=="current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  } 
    
  // store D and M matrix in Dold and Mold
  else if (state=="old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
  }
  
  // unknown conversion
  else
    dserror("ERROR: StoreDM: Unknown conversion requested!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Print current active set to screen                        popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::PrintActiveSet()
{
  
  // get input parameter ctype
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
    
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // friction
      double tz;
      double tjump;
           
      if(ftype=="tresca")
      {     
        // compute tangential part of Lagrange multiplier
        tz = cnode->txi()[0]*cnode->lm()[0] + cnode->txi()[1]*cnode->lm()[1];
           
        // compute tangential part of Lagrange multiplier
        tjump = cnode->txi()[0]*cnode->jump()[0] + cnode->txi()[1]*cnode->jump()[1];
      }  
      
      // get D.B.C. status of current node
      bool dbc = cnode->IsDbc();
      
      // print nodes of inactive set *************************************
      if (cnode->Active()==false)
        cout << "INACTIVE: " << dbc << " " << gid << " " << wgap << " " << nz << endl;
      
      // print nodes of active set ***************************************
      else
      {
        if (ctype != "frictional")
          cout << "ACTIVE:   " << dbc << " " << gid << " " << nz <<  " " << wgap << endl;
        
        else
        {
          if (cnode->Slip() == false)
            cout << "ACTIVE:   " << dbc << " " << gid << " " << nz <<  " " << wgap << " STICK" << " " << tz << endl;
          else
            cout << "ACTIVE:   " << dbc << " " << gid << " " << nz << " " << wgap << " SLIP" << " " << tjump << endl;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::ManagerBase::VisualizeGmsh(const int step, const int iter)
{
  //check for frictional contact
  bool fric = false;
  string ftype   = scontact_.get<string>("friction type","none");
  if (ftype == "tresca" || ftype == "coulomb") fric=true;
  
  // visualization with gmsh
  for (int i=0;i<(int)interface_.size();++i)
    interface_[i]->VisualizeGmsh(interface_[i]->CSegs(),step,iter,fric);
}

#endif  // #ifdef CCADISCRET
