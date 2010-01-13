/*!----------------------------------------------------------------------
\file meshtying_abstract_strategy.cpp

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
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "Epetra_SerialComm.h"
#include "meshtying_abstract_strategy.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::MtAbstractStrategy::MtAbstractStrategy(DRT::Discretization& discret, RCP<Epetra_Map> problemrowmap,
                                                Teuchos::ParameterList params,
                                                vector<RCP<CONTACT::MtInterface> > interface,
                                                int dim, RCP<Epetra_Comm> comm, double alphaf) :
MORTAR::StrategyBase(problemrowmap,params,dim,comm,alphaf),
probdiscret_(discret),
interface_(interface)
{
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------                     

  // merge interface maps to global maps
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_, interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_, interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, interface_[i]->MasterRowDofs());
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map) 
  gndofrowmap_ = LINALG::SplitMap(*problemrowmap_, *gsdofrowmap_);
  gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);

  
  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------   

  // setup Lagrange muliplier vectors
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::MtAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 | set current deformation state                             popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::SetState(const string& statename,
                                           const RCP<Epetra_Vector> vec)
{
  if (statename=="displacement")
  {
    // set state on interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->SetState(statename, vec);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration             popp 12/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::MortarCoupling(const RCP<Epetra_Vector> dis)
{ 
  // set state
  SetState("displacement",dis);
  
  //********************************************************************
  // initialize and evaluate interfaces
  //********************************************************************
  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();
    
    // evaluate interfaces with 'nonlinear'==false
    // (this is meshyting, thus linear mortar coupling is sufficient)
    interface_[i]->Evaluate(false);
  }
  
  //********************************************************************
  // initialize and evaluate global mortar stuff
  //********************************************************************
  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_       = LINALG::CreateVector(*gsdofrowmap_, true);

  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDM(*dmatrix_,*mmatrix_);
    interface_[i]->AssembleG(*g_);
  }
  
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Evaluate(RCP<LINALG::SparseOperator>& kteff,
                                           RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis)
{
  // trivial (no choice as for contact)
  EvaluateMeshtying(kteff,feff,dis);
  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers into MtNode                    popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n interfaces!");

    // get global quantity to be stored in nodes
    RCP<Epetra_Vector> vectorglobal = null;
    switch (type)
    {
      case MORTAR::StrategyBase::lmcurrent:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case MORTAR::StrategyBase::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      default:
        dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
    } // switch

    // export global quantity to current interface slave dof row map
    RCP<Epetra_Map> sdofrowmap = interface_[i]->SlaveRowDofs();
    RCP<Epetra_Vector> vectorinterface = rcp(new Epetra_Vector(*sdofrowmap));

    if (vectorglobal != null) LINALG::Export(*vectorglobal, *vectorinterface);
    else dserror("ERROR: StoreNodalQuantities: Null vector handed in!");

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = mtnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      vector<int> locindex(dim);

      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(mtnode->Dofs()[dof]);
        if (locindex[dof]<0) dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch(type)
        {
        case MORTAR::StrategyBase::lmcurrent:
        {
          mtnode->lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmold:
        {
          mtnode->lmold()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmuzawa:
        {
          mtnode->lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmupdate:
        {
          // throw a dserror if node is Active and DBC
          if (mtnode->IsDbc())
            dserror("ERROR: Slave Node %i is active and at the same time carries D.B.C.s!", mtnode->Id());
     
          // store updated LM into node
          mtnode->lm()[dof] = (*vectorinterface)[locindex[dof]];
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
 |  Store dirichlet B.C. status into MtNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k=0;k<mtnode->NumDof();++k)
      {
        int currdof = mtnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid>=0)
        {
          mtnode->SetDbc() = true;
          break;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update meshtying at end of time step                       popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Update(int istep)
{
  // store Lagrange multipliers
  // (we need this for interpolation of the next generalized mid-point)
  zold_->Update(1.0,*z_,0.0);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for meshtying                    popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::DoReadRestart(IO::DiscretizationReader& reader,
                                                RCP<Epetra_Vector> dis)
{
  // set restart displacement state
  // evaluate interface and restart mortar quantities
  MortarCoupling(dis);
 
  // read restart information on Lagrange multipliers
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  reader.ReadVector(LagrMult(),"lagrmultold");
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  reader.ReadVector(LagrMultOld(),"lagrmultold");
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);
  
  // only for Augmented strategy
  INPAR::CONTACT::SolvingStrategy st = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
  if (st == INPAR::CONTACT::solution_auglag)
  {
    zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    reader.ReadVector(LagrMultUzawa(),"lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }
    
  return;
}

/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::InterfaceForces(RCP<Epetra_Vector> fresm)
{
  // Note that we ALWAYS use a TR-like approach to compute the interface
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}

  /*// CHECK OF INTERFACE FORCE EQUILIBRIUM
  
  // compute two subvectors of fc each via Lagrange multipliers z_n+1, z_n
  RCP<Epetra_Vector> fcslavetemp = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  RCP<Epetra_Vector> fcslavetempend = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertempend = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  dmatrix_->Multiply(false, *z_, *fcslavetemp);
  mmatrix_->Multiply(true, *z_, *fcmastertemp);
  dmatrix_->Multiply(false, *zold_, *fcslavetempend);
  mmatrix_->Multiply(true, *zold_, *fcmastertempend);
  
  // export the interface forces to full dof layout
  RCP<Epetra_Vector> fcslave = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcslaveend = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmasterend = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);
  LINALG::Export(*fcslavetempend, *fcslaveend);
  LINALG::Export(*fcmastertempend, *fcmasterend);

  // build slave and master interface force vector
  RCP<Epetra_Vector> fresmslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fresmmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  LINALG::Export(*fresm,*fresmslave);
  LINALG::Export(*fresm,*fresmmaster);

  // interface forces and moments
  vector<double> gfcs(3);
  vector<double> ggfcs(3);
  vector<double> gfcm(3);
  vector<double> ggfcm(3);
  vector<double> gmcs(3);
  vector<double> ggmcs(3);
  vector<double> gmcm(3);
  vector<double> ggmcm(3);

  vector<double> gmcsnew(3);
  vector<double> ggmcsnew(3);
  vector<double> gmcmnew(3);
  vector<double> ggmcmnew(3);

  // weighted gap vector
  RCP<Epetra_Vector> gapslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: InterfaceForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = mtnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = mtnode->xspatial()[d];
        lm[d] = mtnode->Dofs()[d];
        lmowner[d] = mtnode->Owner();
      }
      LINALG::Assemble(*gapslave,posnode,lm,lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j=0;j<interface_[i]->MasterRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcmastertemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: InterfaceForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = mtnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = mtnode->xspatial()[d];
        lm[d] = mtnode->Dofs()[d];
        lmowner[d] = mtnode->Owner();
      }
      LINALG::Assemble(*gapmaster,posnode,lm,lmowner);
    }
  }

  // weighted gap
  RCP<Epetra_Vector> gapslavefinal  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmasterfinal = rcp(new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false,*gapslave,*gapslavefinal);
  mmatrix_->Multiply(false,*gapmaster,*gapmasterfinal);
  RCP<Epetra_Vector> gapfinal = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0,*gapslavefinal,0.0);
  gapfinal->Update(-1.0,*gapmasterfinal,1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      vector<double> lm(3);
      vector<double> nodegaps(3);
      vector<double> nodegapm(3);

      // LMs and gaps
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(mtnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: InterfaceForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = mtnode->lm()[d];
      }

      // moments
      vector<double> nodemoments(3);
      vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1]*lm[2]-nodegaps[2]*lm[1];
      nodemoments[1] = nodegaps[2]*lm[0]-nodegaps[0]*lm[2];
      nodemoments[2] = nodegaps[0]*lm[1]-nodegaps[1]*lm[0];
      nodemomentm[0] = nodegapm[1]*lm[2]-nodegapm[2]*lm[1];
      nodemomentm[1] = nodegapm[2]*lm[0]-nodegapm[0]*lm[2];
      nodemomentm[2] = nodegapm[0]*lm[1]-nodegapm[1]*lm[0];
      for (int d=0;d<3;++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }

      //cout << "NORMAL: " << mtnode->n()[0] << " " << mtnode->n()[1] << " " << mtnode->n()[2] << endl;
      //cout << "LM:     " << lm[0] << " " << lm[1] << " " << lm[2] << endl;
      //cout << "GAP:    " << nodegaps[0]-nodegapm[0] << " " << nodegaps[1]-nodegapm[1] << " " << nodegaps[2]-nodegapm[2] << endl;
    }
  }

  // summing up over all processors
  for (int i=0;i<3;++i)
  {
    Comm().SumAll(&gfcs[i],&ggfcs[i],1);
    Comm().SumAll(&gfcm[i],&ggfcm[i],1);
    Comm().SumAll(&gmcs[i],&ggmcs[i],1);
    Comm().SumAll(&gmcm[i],&ggmcm[i],1);
    Comm().SumAll(&gmcsnew[i],&ggmcsnew[i],1);
    Comm().SumAll(&gmcmnew[i],&ggmcmnew[i],1);
  }

  // output
  if (Comm().MyPID()==0)
  {
    cout << "Slave Meshtying Force Vector:       " << ggfcs[0] << " " << ggfcs[1] << " " << ggfcs[2] << endl;
    cout << "Slave Meshtying Moment Vector:      " << ggmcs[0] << " " << ggmcs[1] << " " << ggmcs[2] << endl;
    cout << "Slave Meshtying Moment Vector (v2): " << ggmcsnew[0] << " " << ggmcsnew[1] << " " << ggmcsnew[2] << endl;
  }

  if (Comm().MyPID()==0)
  {
    cout << "Master Meshtying Force Vector:       " << ggfcm[0] << " " << ggfcm[1] << " " << ggfcm[2] << endl;
    cout << "Master Meshtying Moment Vector:      " << ggmcm[0] << " " << ggmcm[1] << " " << ggmcm[2] << endl;
    cout << "Master Meshtying Moment Vector (v2): " << ggmcmnew[0] << " " << ggmcmnew[1] << " " << ggmcmnew[2] << endl;
  }
  // CHECK OF INTERFACE FORCE EQUILIBRIUM ----------------------------------
  */

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------- CONTACT::MtAbstractStrategy\n"
       << "Meshtying interfaces: " << (int)interface_.size() << endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    cout << *(interface_[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtAbstractStrategy::PrintActiveSet()
{
  // print message
  if (Comm().MyPID()==0)
    cout << "Meshtying interface-------------------------------------------------------------\n";
  Comm().Barrier();

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MtNode* mtnode = static_cast<MtNode*>(node);

      // compute Lagrange multiplier
      double lm[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k) lm[k] = mtnode->lmold()[k];
      
      // print nodes of active set *************************************
      printf("ACTIVE: %d \t lm[0]: %e \t lm[1]: %e \t lm[2]: %e \n",gid,lm[0],lm[1],lm[2]);
      fflush(stdout);
    }
  }

  Comm().Barrier();

  return;
}

#endif // CCADISCRET
