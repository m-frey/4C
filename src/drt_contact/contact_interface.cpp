/*!----------------------------------------------------------------------
\file contact_interface.cpp
\brief One contact interface

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

#ifndef PARALLEL
#include "Epetra_SerialComm.h"
#endif
#include <Epetra_CrsMatrix.h>
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_integrator.H"
#include "contact_coupling2d.H"
#include "contact_coupling3d.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "selfcontact_binarytree.H"
#include "../drt_mortar/mortar_binarytree.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoInterface::CoInterface(const int id, const Epetra_Comm& comm,
                                  const int dim,
                                  const Teuchos::ParameterList& icontact,
                                  bool selfcontact,
                                  bool redundant) :
MORTAR::MortarInterface(id,comm,dim,icontact,redundant),
selfcontact_(selfcontact),
friction_(false)
{
  // set frictional contact status
  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(icontact,"FRICTION");
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;

  // check for redundant slave storage
  // (needed for self contact but not wanted for general contact)
  if (selfcontact_ && !redundant)
    dserror("ERROR: We need redundant interface storage for self contact");
  if (!selfcontact_ && redundant)
    dserror("ERROR: We do not want redundant interface storage for contact");

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CoInterface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
    os << "Contact ";
  MORTAR::MortarInterface::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  add contact node (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddCoNode(RCP<CONTACT::CoNode> cnode)
{
  idiscret_->AddNode(cnode);
  return;
}

/*----------------------------------------------------------------------*
 |  add contact element (public)                             mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddCoElement(RCP<CONTACT::CoElement> cele)
{
  // check for quadratic 3d slave elements to be modified
  if (cele->IsSlave() && (cele->Shape()==DRT::Element::quad8 || cele->Shape()==DRT::Element::tri6))
    quadslave3d_=true;

  idiscret_->AddElement(cele);
  return;
}

/*----------------------------------------------------------------------*
 |  redistribute contact interface (public)                   popp 08/10|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::Redistribute(int index)
{
  // we need PARALLEL and PARMETIS defined for this
#if !defined(PARALLEL) || !defined(PARMETIS)
  dserror("ERROR: Redistribution of mortar interface needs PARMETIS");
#endif

  // make sure we are supposed to be here
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(IParams(),"PARALLEL_REDIST")==INPAR::MORTAR::parredist_none)
    dserror("ERROR: You are not supposed to be here...");
  
  // some local variables
  RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  const int myrank  = comm->MyPID();
  const int numproc = comm->NumProc();
  Epetra_Time time(*comm);
  set<int>::const_iterator iter;

  // vector containing all proc ids
  vector<int> allproc(numproc);
  for (int i=0; i<numproc; ++i) allproc[i] = i;

  //**********************************************************************
  // (1) SLAVE splitting in close / non-close parts
  //**********************************************************************
  // perform contact search (still with non-optimal distribution)
  Initialize();
  if (SearchAlg()==INPAR::MORTAR::search_bfele)           EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg()==INPAR::MORTAR::search_binarytree) EvaluateSearchBinarytree();
  else                                                    dserror("ERROR: Invalid search algorithm");

  // split slave element row map and build redundant vector of
  // all close / non-close slave node ids on all procs
  vector<int> closeele, noncloseele;
  vector<int> localcns, localfns;

  // loop over all row elements to gather the local information
  for (int i=0; i<SlaveRowElements()->NumMyElements(); ++i)
  {
    // get element
    int gid = SlaveRowElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
    MORTAR::MortarElement* cele = static_cast<MORTAR::MortarElement*>(ele);

    // store element id and adjacent node ids
    int close = cele->MoData().NumSearchElements();
    if (close > 0)
    {
      closeele.push_back(gid);
      for (int k=0;k<cele->NumNode();++k) localcns.push_back(cele->NodeIds()[k]);
    }
    else
    {
      noncloseele.push_back(gid);
      for (int k=0;k<cele->NumNode();++k) localfns.push_back(cele->NodeIds()[k]);
    }
  }

  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  for (int i=0;i<SlaveColElements()->NumMyElements();++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

    mele->MoData().SearchElements().resize(0);
  }

  // we need an arbitrary preliminary element row map
  RCP<Epetra_Map> scroweles  = rcp(new Epetra_Map(-1,(int)closeele.size(),&closeele[0],0,Comm()));
  RCP<Epetra_Map> sncroweles = rcp(new Epetra_Map(-1,(int)noncloseele.size(),&noncloseele[0],0,Comm()));
  RCP<Epetra_Map> mroweles   = rcp(new Epetra_Map(*MasterRowElements()));

  // check for consistency
  if (scroweles->NumGlobalElements()==0 && sncroweles->NumGlobalElements()==0)
    dserror("ERROR: Redistribute: Both slave sets (close/non-close) are empty");

  //**********************************************************************
  // (2) SPECIAL CASES and output to screen
  //**********************************************************************
  // print element overview
  if (!myrank)
  {
    int cl = scroweles->NumGlobalElements();
    int ncl = sncroweles->NumGlobalElements();
     int ma = mroweles->NumGlobalElements();
    cout << "Element overview: " << cl << " / " << ncl << " / " << ma << "  (close-S / non-close-S / M)";
  }

  // print old parallel distribution
  PrintParallelDistribution(index);

  // use simple base class method if there are ONLY close elements
  // (return value TRUE, because redistribution performed)
  if (scroweles->NumGlobalElements()==0 || sncroweles->NumGlobalElements()==0)
  {
    MORTAR::MortarInterface::Redistribute();
    return true;
  }

  //**********************************************************************
  // (3a) PREPARATIONS decide how many procs are used
  //**********************************************************************
  // first we assume that all procs will be used
  int scproc = numproc;
  int sncproc = numproc;
  int mproc = numproc;
  
  // minimum number of elements per proc
  int minele = IParams().get<int>("MIN_ELEPROC");
  
  // calculate real number of procs to be used
  if (minele > 0)
  {
    scproc  = static_cast<int>((scroweles->NumGlobalElements()) / minele);
    sncproc = static_cast<int>((sncroweles->NumGlobalElements()) / minele);
    mproc   = static_cast<int>((mroweles->NumGlobalElements()) / minele);
    if (scroweles->NumGlobalElements() < 2*minele)  scproc = 1;
    if (sncroweles->NumGlobalElements() < 2*minele) sncproc = 1;
    if (mroweles->NumGlobalElements() < 2*minele)   mproc = 1;
    if (scproc > numproc)  scproc = numproc;
    if (sncproc > numproc) sncproc = numproc;
    if (mproc > numproc)   mproc = numproc;
  }

  // print message
  if (!myrank)
  {
    cout << "\nProcs used for redistribution: " << scproc << " / " << sncproc << " / " << mproc << " (close-S / non-close-S / M)";
    cout << "\nRedistributing interface using 3-PARMETIS.......";
  }

  //**********************************************************************
  // (3b) PREPARATIONS build initial node graph
  //**********************************************************************
  // create graph object
  Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*SlaveRowNodes(),108,false));

  // loop over all row nodes to fill graph
  for (int k=0;k<SlaveRowNodes()->NumMyElements();++k)
  {
    int gid = SlaveRowNodes()->GID(k);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);

    // find adjacent elements first
    for (int k=0;k<node->NumElement();++k)
    {
      // store adjacent nodes
      DRT::Element* ele = node->Elements()[k];
      int numnode = ele->NumNode();
      vector<int> nodeids(numnode);
      for (int n=0;n<numnode;++n) nodeids[n] = ele->NodeIds()[n];

      int err = graph->InsertGlobalIndices(gid,numnode,&nodeids[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
      if (err==1) dserror("graph->InsertGlobalIndices returned %d",err);
    }
  }

  // fill graph and optimize storage
  graph->FillComplete();
  graph->OptimizeStorage();

  //**********************************************************************
  // (4) CLOSE SLAVE redistribution
  //**********************************************************************
  RCP<Epetra_Map> scrownodes = Teuchos::null;
  RCP<Epetra_Map> sccolnodes = Teuchos::null;

  // build redundant vector of all close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  vector<int> globalcns;
  set<int> setglobalcns;
  vector<int> scnids;
  LINALG::Gather<int>(localcns,globalcns,numproc,&allproc[0],Comm());
  for (int i=0;i<(int)globalcns.size();++i) setglobalcns.insert(globalcns[i]);
  for (iter=setglobalcns.begin();iter!=setglobalcns.end();++iter) scnids.push_back(*iter);

  //**********************************************************************
  // call PARMETIS (again with #ifdef to be on the safe side)
#if defined(PARALLEL) && defined(PARMETIS)
  // old version
  //DRT::UTILS::PartUsingParMetis(idiscret_,scroweles,scrownodes,sccolnodes,scnids,numproc,scproc,comm,time,false);
  // new version
  DRT::UTILS::PartUsingParMetis(idiscret_,scroweles,scrownodes,sccolnodes,comm,false);
#endif
  //**********************************************************************

  //**********************************************************************
  // (5) NON-CLOSE SLAVE redistribution
  //**********************************************************************
  RCP<Epetra_Map> sncrownodes = Teuchos::null;
  RCP<Epetra_Map> snccolnodes = Teuchos::null;

  // build redundant vector of all non-close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  vector<int> globalfns;
  set<int> setglobalfns;
  vector<int> sncnids;
  LINALG::Gather<int>(localfns,globalfns,numproc,&allproc[0],Comm());
  for (int i=0;i<(int)globalfns.size();++i) setglobalfns.insert(globalfns[i]);
  for (iter=setglobalfns.begin();iter!=setglobalfns.end();++iter) sncnids.push_back(*iter);

  //**********************************************************************
  // call PARMETIS (again with #ifdef to be on the safe side)
#if defined(PARALLEL) && defined(PARMETIS)
  // old version
  //DRT::UTILS::PartUsingParMetis(idiscret_,sncroweles,sncrownodes,snccolnodes,sncnids,numproc,sncproc,comm,time,false);
  // new version
  DRT::UTILS::PartUsingParMetis(idiscret_,sncroweles,sncrownodes,snccolnodes,comm,false);
#endif
  //**********************************************************************

  //**********************************************************************
  // (6) MASTER redistribution
  //**********************************************************************
  RCP<Epetra_Map> mrownodes = Teuchos::null;
  RCP<Epetra_Map> mcolnodes = Teuchos::null;

  // build redundant vector of all master node ids on all procs
  // (do not include crosspoints / boundary nodes if there are any)
  vector<int> mnids;
  vector<int> mnidslocal(MasterRowNodesNoBound()->NumMyElements());
  for (int i=0; i<MasterRowNodesNoBound()->NumMyElements(); ++i)
    mnidslocal[i] = MasterRowNodesNoBound()->GID(i);
  LINALG::Gather<int>(mnidslocal,mnids,numproc,&allproc[0],Comm());

  //**********************************************************************
  // call PARMETIS (again with #ifdef to be on the safe side)
#if defined(PARALLEL) && defined(PARMETIS)
  // old version
  //DRT::UTILS::PartUsingParMetis(idiscret_,mroweles,mrownodes,mcolnodes,mnids,numproc,mproc,comm,time,false);
  // new version
  DRT::UTILS::PartUsingParMetis(idiscret_,mroweles,mrownodes,mcolnodes,comm,false);
#endif
  //**********************************************************************

  //**********************************************************************
  // (7) Merge global interface node row and column map
  //**********************************************************************
  // merge slave node row map from close and non-close parts
  RCP<Epetra_Map> srownodes = Teuchos::null;

  //----------------------------------CASE 1: ONE OR BOTH SLAVE SETS EMPTY
  if (scrownodes==Teuchos::null || sncrownodes==Teuchos::null)
  {
     dserror("ERROR: Redistribute: You should not be here");
  }
  //-------------------------------------CASE 2: BOTH SLAVE SETS NON-EMPTY
  else
  {
    // find intersection set of close and non-close nodes
    set<int> intersec;
    for (iter=setglobalcns.begin();iter!=setglobalcns.end();++iter)
    {
      set<int>::const_iterator found = setglobalfns.find(*iter);
      if (found!=setglobalfns.end()) intersec.insert(*found);
    }

    // build slave node row map
    vector<int> mygids(scrownodes->NumMyElements() + sncrownodes->NumMyElements());
    int count = scrownodes->NumMyElements();

    // first get GIDs of input scrownodes
    for (int i=0;i<count;++i) mygids[i] = scrownodes->GID(i);

    // then add GIDs of input sncrownodes (only new ones)
    for (int i=0;i<sncrownodes->NumMyElements();++i)
    {
      // check for intersection gid
      // don't do anything for intersection gids (scrownodes dominates!!!)
      set<int>::const_iterator found = intersec.find(sncrownodes->GID(i));
      if (found!=intersec.end()) continue;

      // check for overlap
      if (scrownodes->MyGID(sncrownodes->GID(i)))
        dserror("LINALG::MergeMap: Result map is overlapping");

      // add new GIDs to mygids
      mygids[count]=sncrownodes->GID(i);
      ++count;
    }
    mygids.resize(count);
    sort(mygids.begin(),mygids.end());
    srownodes = rcp(new Epetra_Map(-1,(int)mygids.size(),&mygids[0],0,scrownodes->Comm()));
  }

  // merge interface node row map from slave and master parts
  RCP<Epetra_Map> rownodes = LINALG::MergeMap(srownodes,mrownodes,false);

  // IMPORTANT NOTE:
  // While merging from the two different slave parts of the discretization
  // (close slave, non-close slave) is feasible for the node row map,
  // this is not possible for the node column map. Some necessary
  // information on ghosting at the transition between close and non-close
  // slave region would always be missed! Thus, we reconstruct a
  // suitable slave node column map "by hand" here. This is quite simply
  // done by exporting the initial node graph to the new distribution
  // and by then asking for its column map.

   // create the output graph (with new slave node row map) and export to it
   RCP<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(Copy,*srownodes,108,false));
   Epetra_Export exporter(graph->RowMap(),*srownodes);
   int err = outgraph->Export(*graph,exporter,Add);
   if (err<0) dserror("Graph export returned err=%d",err);

  // trash old graph
  graph=null;

  // call fill complete and optimize storage
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  // get column map from the graph -> build slave node column map
  // (do stupid conversion from Epetra_BlockMap to Epetra_Map)
  const Epetra_BlockMap& bcol = outgraph->ColMap();
  RCP<Epetra_Map> scolnodes = rcp(new Epetra_Map(bcol.NumGlobalElements(),bcol.NumMyElements(),bcol.MyGlobalElements(),0,Comm()));

  // trash new graph
  outgraph=null;
  
  // merge interface node column map from slave and master parts
  RCP<Epetra_Map> colnodes = LINALG::MergeMap(scolnodes,mcolnodes,false);

  //**********************************************************************
  // (8) Get partitioning information into discretization
  //**********************************************************************
  // build reasonable element maps from the already valid and final node maps
  // (note that nothing is actually redistributed in here)
  RCP<Epetra_Map> roweles  = Teuchos::null;
  RCP<Epetra_Map> coleles  = Teuchos::null;
  Discret().BuildElementRowColumn(*rownodes,*colnodes,roweles,coleles);

  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodes);
  Discret().ExportRowElements(*roweles);

  // export nodes and elements to the column map (create ghosting)
  Discret().ExportColumnNodes(*colnodes);
  Discret().ExportColumnElements(*coleles);

  // print message
  if (!myrank) cout << "done!" << endl;

  return true;
}

/*----------------------------------------------------------------------*
 | collect distribution data (public)                         popp 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::CollectDistributionData(int& loadele, int& crowele)
{
  // loop over proc's column slave elements of the interface
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CoElement* selement = static_cast<CoElement*>(ele1);
    
    // bool indicating coupling partners
    bool add = (selement->MoData().NumSearchElements()>0);
    
    // check if this element has any coupling partners and add
    // element ID to input variable loadele if so
    if (add) loadele += 1;
    
    // check if - in addition - the active proc owns this element
    // and add element ID to input variable rowele if so
    if (add && selement->Owner()==Comm().MyPID()) crowele += 1;
  }
  
  return;  
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::CreateSearchTree()
{
  // ***WARNING:*** This is commented out here, as idiscret_->SetState()
  // needs all the procs around, not only the interface local ones!
  // if (!lComm()) return;
  
  // warning
#ifdef MORTARGMSHCTN
  if (Dim()==3 && Comm().MyPID()==0)
  {
    cout << "\n******************************************************************\n";
    cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    cout << "******************************************************************\n";
  }
#endif // #ifdef MORTARGMSHCTN

  // binary tree search
  if (SearchAlg()==INPAR::MORTAR::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (SelfContact())
    {
      // set state in interface to intialize all kinds of quantities
      RCP<Epetra_Vector> zero =rcp(new Epetra_Vector(*idiscret_->DofRowMap()));
      SetState("displacement",zero);

      // create fully overlapping map of all contact elements
      RCP<Epetra_Map> elefullmap = LINALG::AllreduceEMap(*idiscret_->ElementRowMap());

      // create binary tree object for self contact search
      // (TODO: NOTE THAT SELF CONTACT SEARCH IS NOT YET PARALLELIZED!)
      binarytreeself_ = rcp(new CONTACT::SelfBinaryTree(Discret(),lComm(),elefullmap,Dim(),SearchParam()));

    }
    //*****TWO BODY CONTACT*****
    else
    {
      // get out of here if not participating in interface
      if (!lComm()) return;

      // create fully overlapping map of all master elements
      RCP<Epetra_Map> melefullmap = LINALG::AllreduceEMap(*melerowmap_);
      
      // create binary tree object for contact search and setup tree
      binarytree_ = rcp(new MORTAR::BinaryTree(Discret(),selecolmap_,melefullmap,Dim(),SearchParam()));

      // initialize active contact nodes via binarytree
      // binarytree_->SearchContactInit(binarytree_->Sroot(), binarytree_->Mroot());
    }
  }
  
  // no binary tree search
  else
  {
    if (SelfContact()) dserror("ERROR: Binarytree search needed for self contact");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)
  
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CoNode* node = static_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));

    // reset feasible projection status
    node->HasProj() = false;
    
    if (friction_)
    {  
      FriNode* frinode = static_cast<FriNode*>(node);
  
      // reset nodal mechanical dissipation
      frinode->MechDiss() = 0.0;
      
      // reset matrix B quantities
      frinode->GetBNodes().clear();
      
      // reset nodal B maps
      for (int j=0;j<(int)((frinode->GetB()).size());++j)
        (frinode->GetB())[j].clear();
      
      (frinode->GetB()).resize(0);

    }
  }
  
  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i=0;i<SlaveColNodesBound()->NumMyElements();++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    
    //reset nodal normal and tangents and jumps
    for (int j=0;j<3;++j)
    {
      cnode->MoData().n()[j]=0.0;
      cnode->CoData().txi()[j]=0.0;
      cnode->CoData().teta()[j]=0.0;
    }

    // reset nodal Mortar maps
    for (int j=0;j<(int)((cnode->MoData().GetD()).size());++j)
      (cnode->MoData().GetD())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetM()).size());++j)
      (cnode->MoData().GetM())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetMmod()).size());++j)
      (cnode->MoData().GetMmod())[j].clear();
    
    (cnode->MoData().GetD()).resize(0);
    (cnode->MoData().GetM()).resize(0);
    (cnode->MoData().GetMmod()).resize(0);
    
    // reset derivative maps of normal vector
    for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
      (cnode->CoData().GetDerivN())[j].clear();
    (cnode->CoData().GetDerivN()).resize(0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
      (cnode->CoData().GetDerivTxi())[j].clear();
    (cnode->CoData().GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
      (cnode->CoData().GetDerivTeta())[j].clear();
    (cnode->CoData().GetDerivTeta()).resize(0);
    
    // reset derivative map of Mortar matrices
    (cnode->CoData().GetDerivD()).clear();
    (cnode->CoData().GetDerivM()).clear();
    
    // reset nodal weighted gap and derivative
    cnode->CoData().Getg() = 1.0e12;
    (cnode->CoData().GetDerivG()).clear();
    
    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((cnode->CoData().GetDerivZ()).size()); ++j)
      (cnode->CoData().GetDerivZ())[j].clear();
    (cnode->CoData().GetDerivZ()).resize(0);

    if (friction_)
    {  
      FriNode* frinode = static_cast<FriNode*>(cnode);
  
      // reset SNodes and Mnodes
      frinode->FriData().GetSNodes().clear();
      frinode->FriData().GetMNodes().clear();

      // reset matrix A quantities
      frinode->FriData().GetANodes().clear();
      
      // reset nodal A maps
      for (int j=0;j<(int)((frinode->FriData().GetA()).size());++j)
        (frinode->FriData().GetA())[j].clear();
      
      (frinode->FriData().GetA()).resize(0);
      
      // reset wear increment
      frinode->FriData().DeltaWear() = 0.0;
      
    }
  }

  //**********************************************************************
  // In general, it is sufficient to reset search candidates only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need to reset the search candidates of
  // all slave elements in the fully overlapping column map there. This
  // is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to reset candidates / search lists
    // (use fully overlapping column map of S+M elements)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      DRT::Element* ele = idiscret_->lColElement(i);
      MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);
    }
  }
  else
  {
    // loop over all elements to reset candidates / search lists
    // (use standard slave column map)
    for (int i=0;i<SlaveColElements()->NumMyElements();++i)
    {
      int gid = SlaveColElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
      MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::SetElementAreas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      MORTAR::MortarElement* element = static_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
      element->InitializeDataContainer();
      element->MoData().Area()=element->ComputeArea();
    }
  }
  else
  {
    // refer call back to base class version
    MORTAR::MortarInterface::SetElementAreas();
  }

  return;
}

/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ExportNodalNormals()
{
  // create empty data objects
  map<int,RCP<Epetra_SerialDenseMatrix> > triad;

  map<int,vector<int> > n_x_key;
  map<int,vector<int> > n_y_key;
  map<int,vector<int> > n_z_key;
  map<int,vector<int> > txi_x_key;
  map<int,vector<int> > txi_y_key;
  map<int,vector<int> > txi_z_key;
  map<int,vector<int> > teta_x_key;
  map<int,vector<int> > teta_y_key;
  map<int,vector<int> > teta_z_key;

  map<int,vector<double> > n_x_val;
  map<int,vector<double> > n_y_val;
  map<int,vector<double> > n_z_val;
  map<int,vector<double> > txi_x_val;
  map<int,vector<double> > txi_y_val;
  map<int,vector<double> > txi_z_val;
  map<int,vector<double> > teta_x_val;
  map<int,vector<double> > teta_y_val;
  map<int,vector<double> > teta_z_val;

  map<int,double>::iterator iter;

  // build info on row map
  for(int i=0; i<snoderowmapbound_->NumMyElements();++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    // fill nodal matrix
    RCP<Epetra_SerialDenseMatrix> loc = rcp(new Epetra_SerialDenseMatrix(3,3));
    (*loc)(0,0) = cnode->MoData().n()[0];
    (*loc)(1,0) = cnode->MoData().n()[1];
    (*loc)(2,0) = cnode->MoData().n()[2];
    (*loc)(0,1) = cnode->CoData().txi()[0];
    (*loc)(1,1) = cnode->CoData().txi()[1];
    (*loc)(2,1) = cnode->CoData().txi()[2];
    (*loc)(0,2) = cnode->CoData().teta()[0];
    (*loc)(1,2) = cnode->CoData().teta()[1];
    (*loc)(2,2) = cnode->CoData().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    vector<map<int,double> >& derivn    = cnode->CoData().GetDerivN();
    vector<map<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    vector<map<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for(iter=derivn[0].begin();iter!=derivn[0].end();++iter)
    {
      n_x_key[gid].push_back(iter->first);
      n_x_val[gid].push_back(iter->second);
    }
    for(iter=derivn[1].begin();iter!=derivn[1].end();++iter)
    {
      n_y_key[gid].push_back(iter->first);
      n_y_val[gid].push_back(iter->second);
    }
    for(iter=derivn[2].begin();iter!=derivn[2].end();++iter)
    {
      n_z_key[gid].push_back(iter->first);
      n_z_val[gid].push_back(iter->second);
    }

    for(iter=derivtxi[0].begin();iter!=derivtxi[0].end();++iter)
    {
      txi_x_key[gid].push_back(iter->first);
      txi_x_val[gid].push_back(iter->second);
    }
    for(iter=derivtxi[1].begin();iter!=derivtxi[1].end();++iter)
    {
      txi_y_key[gid].push_back(iter->first);
      txi_y_val[gid].push_back(iter->second);
    }
    for(iter=derivtxi[2].begin();iter!=derivtxi[2].end();++iter)
    {
      txi_z_key[gid].push_back(iter->first);
      txi_z_val[gid].push_back(iter->second);
    }

    for(iter=derivteta[0].begin();iter!=derivteta[0].end();++iter)
    {
      teta_x_key[gid].push_back(iter->first);
      teta_x_val[gid].push_back(iter->second);
    }
    for(iter=derivteta[1].begin();iter!=derivteta[1].end();++iter)
    {
      teta_y_key[gid].push_back(iter->first);
      teta_y_val[gid].push_back(iter->second);
    }
    for(iter=derivteta[2].begin();iter!=derivteta[2].end();++iter)
    {
      teta_z_key[gid].push_back(iter->first);
      teta_z_val[gid].push_back(iter->second);
    }
  }

  // communicate from slave node row to column map
  DRT::Exporter ex(*snoderowmapbound_,*snodecolmapbound_,Comm());
  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    // extract info
    RCP<Epetra_SerialDenseMatrix> loc = triad[gid];
    cnode->MoData().n()[0]    = (*loc)(0,0);
    cnode->MoData().n()[1]    = (*loc)(1,0);
    cnode->MoData().n()[2]    = (*loc)(2,0);
    cnode->CoData().txi()[0]  = (*loc)(0,1);
    cnode->CoData().txi()[1]  = (*loc)(1,1);
    cnode->CoData().txi()[2]  = (*loc)(2,1);
    cnode->CoData().teta()[0] = (*loc)(0,2);
    cnode->CoData().teta()[1] = (*loc)(1,2);
    cnode->CoData().teta()[2] = (*loc)(2,2);

    // extract derivative info
    vector<map<int,double> >& derivn    = cnode->CoData().GetDerivN();
    vector<map<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    vector<map<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for (int k=0;k<(int)(derivn.size());++k)
      derivn[k].clear();
    derivn.resize(3);
    for (int k=0;k<(int)(derivtxi.size());++k)
      derivtxi[k].clear();
    derivtxi.resize(3);
    for (int k=0;k<(int)(derivteta.size());++k)
      derivteta[k].clear();
    derivteta.resize(3);

    for (int k=0;k<(int)(n_x_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k=0;k<(int)(n_y_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k=0;k<(int)(n_z_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k=0;k<(int)(txi_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k=0;k<(int)(txi_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k=0;k<(int)(txi_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k=0;k<(int)(teta_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k=0;k<(int)(teta_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k=0;k<(int)(teta_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  // print nodal normals
  /*for (int p=0;p<Comm().NumProc();++p)
  {
    // one proc after the other
    if (p==Comm().MyPID())
    {
      cout << "\n*****\nPROC " << p << "\n*****" << endl;
      for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
      {
        int gid = snodecolmapbound_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cnode = static_cast<CoNode*>(node);

        // print normal and tangents at each slave node
        cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " Normal: " << cnode->MoData().n()[0]
             << " " << cnode->MoData().n()[1] << " " << cnode->MoData().n()[2] << endl;
        cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TXi: " << cnode->CoData().txi()[0]
             << " " << cnode->CoData().txi()[1] << " " << cnode->CoData().txi()[2] << endl;
        cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TEta: " << cnode->CoData().teta()[0]
             << " " << cnode->CoData().teta()[1] << " " << cnode->CoData().teta()[2] << endl;

        // print linearizations at each slave node
        cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinN: ";
        for (iter=cnode->CoData().GetDerivN()[0].begin();iter!=cnode->CoData().GetDerivN()[0].end();++iter)
          cout << "\n" << iter->first << "\t" << iter->second;
        cout << endl;
        cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinTxi: ";
        for (iter=cnode->CoData().GetDerivTxi()[0].begin();iter!=cnode->CoData().GetDerivTxi()[0].end();++iter)
          cout << "\n" << iter->first << "\t" << iter->second;
        cout << endl;
        cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinTeta: ";
        for (iter=cnode->CoData().GetDerivteta()[0].begin();iter!=cnode->CoData().GetDerivTeta()[0].end();++iter)
          cout << "\n" << iter->first << "\t" << iter->second;
        cout << endl;
      }
      cout << endl << endl;
    }

    // barrier
    Comm().Barrier();
  }*/

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::EvaluateSearchBinarytree()
{
  // ***WARNING:*** This is commented out here, as UpdateMasterSlaveSets()
  // needs all the procs around, not only the interface local ones!
  // if (!lComm()) return true;

  // *********************************************************************
  // Possible versions for self contact:
  // *********************************************************************
  //
  // 1) Combined Update and Search
  // -> In this case we have to call SearchContactCombined(), which
  //    does both top-down update (where necessary) and search. Then
  //    the dynamics master/slave assignment routine UpdateMasterSlaveSets()
  //    is called and the new slave nodes' data containers are initialized.
  //
  // 2) Separate Update and Search
  // -> In this case we have to call SearchContactSeparate(), which
  //    does both bottom-up update (on whole interface) and search. Then
  //    the dynamics master/slave assignment routine UpdateMasterSlaveSets()
  //    is called and the new slave nodes' data containers are initialized.
  //
  // *********************************************************************
  if (SelfContact())
  {
    // calculate minimal element length
    binarytreeself_->SetEnlarge(false);

    // update and search for contact with a combined algorithm
    //binarytreeself_->SearchContactCombined();

    // update and search for contact with separate algorithms
    binarytreeself_->SearchContactSeparate();

    // update master/slave sets of interface
    UpdateMasterSlaveSets();
    
    // initialize node data container
    // (include slave side boundary nodes / crosspoints)
    for (int i=0; i<SlaveColNodesBound()->NumMyElements(); ++i)
    {
      int gid = SlaveColNodesBound()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
      MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(node);

      // initialize container if not yet initialized before
      mnode->InitializeDataContainer();
    }

    // no initialization of element data container as this would
    // possibly destroy the information on search elements again
    // (this was already done in SetElementAreas())
  }

  // *********************************************************************
  // Possible versions for 2-body contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchContactCombined(), which
  //    does both top-down update (where necessary) and search.
  //
  // 2) Separate Update and Contact Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchContactSeparate(). Of course, the bottom-up
  //    update makes more sense here. For very large contact problems,
  //    this version is preferable and thus chosen as default.
  //
  // *********************************************************************
  else
  {
    // get out of here if not participating in interface
    if (!lComm()) return true;
    
    // calculate minimal element length
    binarytree_->SetEnlarge(false);

    // update tree in a top down way
    //binarytree_->UpdateTreeTopDown();

    // update tree in a bottom up way
    binarytree_->UpdateTreeBottomUp();

#ifdef MORTARGMSHCTN
    for (int i=0;i<(int)(binarytree_->CouplingMap().size());i++)
      binarytree_->CouplingMap()[i].clear();
    binarytree_->CouplingMap().clear();
    binarytree_->CouplingMap().resize(2);
#endif // #ifdef MORTARGMSHCTN

    // search for contact with a separate algorithm
    binarytree_->SearchSeparate();

    // search for contact with a combined algorithm
    //binarytree_->SearchCombined();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate Mortar matrix D on slave element (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateSlave(MORTAR::MortarElement& sele)
{
  //**********************************************************************
  dserror("ERROR: IntegrateSlave method is outdated!");
  //**********************************************************************

  // create a CONTACT integrator instance with correct NumGP and Dim
  CONTACT::CoIntegrator integrator(shapefcn_,sele.Shape());

  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // do the element integration (integrate and linearize D)
  int nrow = sele.NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  integrator.IntegrateDerivSlave2D3D(sele,sxia,sxib,dseg);

  // do the assembly into the slave nodes
  integrator.AssembleD(Comm(),sele,*dseg);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlaps     popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateCoupling(MORTAR::MortarElement* sele,
                                             vector<MORTAR::MortarElement*> mele)
{
  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim()==2)
  {
    // new coupling method -> loop over all master elements is hidden
    // inside the Coupling2d class (in EvaluateCoupling method)

    // *************************************************** linear 2D ***
    // ************************************************ quadratic 2D ***
    // neither quadratic interpolation nor mixed linear and quadratic
    // interpolation need any special treatment in the 2d case
    
    // create instance of coupling class
    CONTACT::CoCoupling2d coup(shapefcn_,Discret(),Dim(),sele,mele);
    
    // do coupling
    coup.EvaluateCoupling();
  }
  // ************************************************************** 3D ***
  else if (Dim()==3)
  {
    bool auxplane = DRT::INPUT::IntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

    // loop over all master elements associated with this slave element
    for (int m=0;m<(int)mele.size();++m)
    {
      // ************************************************** quadratic 3D ***
      // also treats the mixed linear and quadratic interpolation case
      if (sele->IsQuad3d() || mele[m]->IsQuad3d())
      {
        // only for auxiliary plane 3D version
        if (!auxplane) dserror("ERROR: Quadratic 3D coupling only for AuxPlane case!");
  
        // build linear integration elements from quadratic MortarElements
        vector<RCP<MORTAR::IntElement> > sauxelements(0);
        vector<RCP<MORTAR::IntElement> > mauxelements(0);
        SplitIntElements(*sele,sauxelements);
        SplitIntElements(*mele[m],mauxelements);
  
        // get LM interpolation and testing type
        INPAR::MORTAR::LagMultQuad3D lmtype =
          DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad3D>(IParams(),"LAGMULT_QUAD3D");
              
        // loop over all IntElement pairs for coupling
        for (int i=0;i<(int)sauxelements.size();++i)
        {
          for (int j=0;j<(int)mauxelements.size();++j)
          {
            // create instance of coupling class
            CONTACT::CoCoupling3dQuad coup(shapefcn_,Discret(),Dim(),true,auxplane,
                          *sele,*mele[m],*sauxelements[i],*mauxelements[j],lmtype);
            // do coupling
            coup.EvaluateCoupling();
          }
        }
      }
  
      // ***************************************************** linear 3D ***
      else
      {
        // create instance of coupling class
        CONTACT::CoCoupling3d coup(shapefcn_,Discret(),Dim(),false,auxplane,
                                   *sele,*mele[m]);
        // do coupling
        coup.EvaluateCoupling();
      }
    }
  }
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateKappaPenalty(CONTACT::CoElement& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // check for auxiliary plane 3D version
  bool auxplane = DRT::INPUT::IntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

  // ************************************************** quadratic 3D ***
  if (Dim()==3 && sele.IsQuad3d())
  {
    // only for auxiliary plane 3D version
    if (!auxplane) dserror("ERROR: Quadratic 3D contact only for AuxPlane case!");

    // get LM interpolation and testing type
    INPAR::MORTAR::LagMultQuad3D lmtype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad3D>(IParams(),"LAGMULT_QUAD3D");
          
    // build linear integration elements from quadratic CElements
    vector<RCP<MORTAR::IntElement> > sauxelements(0);
    SplitIntElements(sele,sauxelements);

    // different options for mortar integration
    if (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin)
    {
      // do the element integration of kappa and store into gap
      int nrow = sele.NumNode();
      RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

      // create a CONTACT integrator instance with correct NumGP and Dim
      CONTACT::CoIntegrator integrator(shapefcn_,sele.Shape());
      integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

      // do the assembly into the slave nodes
      integrator.AssembleG(Comm(),sele,*gseg);
    }
    
    else if (lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      // integrate each int element seperately
      for (int i=0;i<(int)sauxelements.size();++i)
      {
        // do the int element integration of kappa and store into gap
        int nrow = sauxelements[i]->NumNode();
        RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

        // create a CONTACT integrator instance with correct NumGP and Dim
        CONTACT::CoIntegrator integrator(shapefcn_,sauxelements[i]->Shape());
        integrator.IntegrateKappaPenalty(sele,*(sauxelements[i]),sxia,sxib,gseg,lmtype);

        // do the assembly into the slave nodes
        integrator.AssembleG(Comm(),*(sauxelements[i]),*gseg);
      }
    }

    else
    {
      dserror("ERROR: IntegrateKappaPenalty: Invalid case for 3D mortar contact LM interpolation");
    }
  }

  // *************************************************** other cases ***
  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.NumNode();
    RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

    // create a CONTACT integrator instance with correct NumGP and Dim
    CONTACT::CoIntegrator integrator(shapefcn_,sele.Shape());
    integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

    // do the assembly into the slave nodes
    integrator.AssembleG(Comm(),sele,*gseg);
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateRelMov(const RCP<Epetra_Vector> xsmod,
                                          const RCP<LINALG::SparseMatrix> dmatrixmod,
                                          const RCP<LINALG::SparseMatrix> doldmod)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
 
  if (friction_ == false)
    dserror("Error in CoInterface::EvaluateRelMov(): Only evaluated for frictional contact");

  // parameters
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k=0;k<3;++k)
      nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

    vector <double> jump(dim);
    for(int dim=0;dim<Dim();dim++)
      jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    double kappa = cnode->CoData().Kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_penalty)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             DRT::INPUT::IntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")!=1)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             DRT::INPUT::IntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")==1)
    {
      if((nz - cn*gap > 0) or cnode->Active()) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_auglag)
    {
      if(lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
      dserror("Error in Interface::EvaluateRelMov(): Solution strategy not known!");

    if(activeinfuture==true)
    {
      vector<map<int,double> > dmap = cnode->MoData().GetD();
      vector<map<int,double> > dmapold = cnode->FriData().GetDOld();

      set <int> snodes = cnode->FriData().GetSNodes();

      // check if there are entries in the old D map
      if(dmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old D-Map!");

      map<int,double>::iterator colcurr;
      set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = static_cast<CoNode*>(snode);
        const int* sdofs = csnode->Dofs();

        double dik = (dmap[0])[sdofs[0]];
        double dikold = (dmapold[0])[sdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<csnode->NumDof();++dim)
        {
          int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
          jump[dim]-=(dik-dikold)*(*xsmod)[locid];
        }
      } //  loop over adjacent slave nodes

      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->FriData().GetMOld();

      set <int> mnodescurrent = cnode->FriData().GetMNodes();
      set <int> mnodesold = cnode->FriData().GetMNodesOld();

      // check if there are entries in the M map
      if(mmap.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No M-Map!");
      
      // check if there are entries in the old M map
      if(mmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old M-Map!");

      if(mnodesold.size() <1)
        dserror ("Error in Interface::EvaluateRelMov(): No old M-Set!"); 
      
      set <int> mnodes;
      set<int>::iterator mcurr;

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
           jump[dim]+= (mik-mikold)*(cmnode->xspatial()[dim]);
        }
      } //  loop over master nodes

      // write it to nodes
      for(int dim=0;dim<Dim();dim++)
        cnode->FriData().jump()[dim] = jump[dim];
 
      // linearization of jump vector

      // reset derivative map of jump
      for (int j=0; j<(int)((cnode->FriData().GetDerivJump()).size()); ++j)
        (cnode->FriData().GetDerivJump())[j].clear();
      (cnode->FriData().GetDerivJump()).resize(0);
      
      /*** 01  **********************************************************/
      
      if(dmatrixmod==null)
      {  
        // loop over according slave nodes 
        for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
        {
          int gid = *scurr;
          DRT::Node* snode = idiscret_->gNode(gid);
          if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
          CoNode* csnode = static_cast<CoNode*>(snode);
          const int* sdofs = csnode->Dofs();

          double dik = (dmap[0])[sdofs[0]];
          double dikold=(dmapold[0])[sdofs[0]];

          for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = csnode->Dofs()[dimrow];
            double val = -(dik-dikold);
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }
      // in the 3D quadratic case, the values are obtained from the
      // global matrices Dmod and Doldmod
      else
      {  
        // loop over dimension of the node
        for (int dim=0;dim<cnode->NumDof();++dim)
        {  
          int NumEntries = 0;
          int NumEntriesOld = 0;
          std::vector<double> Values((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> Indices((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<double> ValuesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> IndicesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());

          // row           
          int row = cnode->Dofs()[dim];
      
          // extract entries of this row from matrix
          int err = (dmatrixmod->EpetraMatrix())->ExtractGlobalRowCopy(row,(dmatrixmod->EpetraMatrix())->MaxNumEntries(),NumEntries,&Values[0], &Indices[0]);
          if (err) dserror("ExtractMyRowView failed: err=%d", err);

          int errold = (doldmod->EpetraMatrix())->ExtractGlobalRowCopy(row,(doldmod->EpetraMatrix())->MaxNumEntries(),NumEntriesOld,&ValuesOld[0], &IndicesOld[0]);
          if (errold) dserror("ExtractMyRowView failed: err=%d", err);

          // loop over entries of this vector
          for (int j=0;j<NumEntries;++j)
          {  
            double ValueOld=0;
            bool found = false;
         
            // find value with the same index in vector of Dold
            for (int k=0;k<NumEntriesOld;++k)
            {
              if (Indices[k]==Indices[j])
              {
                ValueOld = ValuesOld[k];
                found = true;
                break;
              }  
            }

            if(found==false or abs(ValueOld) < 1e-12)
              dserror("Error in EvaluareRelMov(): No old D value exists");
            
            // write to node
            cnode->AddDerivJumpValue(dim,Indices[j],Values[j]-ValueOld);
          }  
        }
      }
 
      /*** 02  **********************************************************/
      // loop over according master nodes 
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold=(mmapold[0])[mdofs[0]];

        for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
        {
            int col = cmnode->Dofs()[dimrow];
            double val = (mik-mikold);
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
        }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
      map<int,map<int,double> >& ddmap = cnode->CoData().GetDerivD();
      map<int,map<int,double> >::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr=ddmap.begin();dscurr!=ddmap.end();++dscurr)
      {
        int gid = dscurr->first;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = static_cast<CoNode*>(snode);

        // compute entry of the current stick node / slave node pair
        map<int,double>& thisdmmap = cnode->CoData().GetDerivD(gid);

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dim=0;dim<cnode->NumDof();++dim)
          {
            int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
            double val =-colcurr->second*(*xsmod)[locid];
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dim,col,val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
      map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
      map<int,map<int,double> >::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
      {
        int gid = dmcurr->first;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = static_cast<CoNode*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
        map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);
        
        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            double val =colcurr->second*mxi[dimrow];
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }
    } // active nodes
  } // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble relative movement / jump (global)             gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRelMov(Epetra_Vector& jumpglobal)
{
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    int dim = cnode->NumDof();
    double* jump = cnode->FriData().jump();

    Epetra_SerialDenseVector jumpnode(dim);
    vector<int> jumpdof(dim);
    vector<int> jumpowner(dim);

    for( int k=0; k<dim; ++k )
    {
      jumpnode(k) = jump[k];
      jumpdof[k] = cnode->Dofs()[k];
      jumpowner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(jumpglobal, jumpnode, jumpdof, jumpowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble slave coordinates (xs)                        gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleSlaveCoord(RCP<Epetra_Vector>& xsmod)
{
   
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    int dim = cnode->NumDof();

    Epetra_SerialDenseVector xspatial(dim);
    vector<int> dof(dim);
    vector<int> owner(dim);

    for( int k=0; k<dim; ++k )
    {
      xspatial(k) = cnode->xspatial()[k];
      dof[k] = cnode->Dofs()[k];
      owner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xsmod, xspatial, dof, owner);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateTangentNorm(double& cnormtan)
{
  cnormtan=0;

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some information form node
    double* n = cnode->MoData().n();
    int dim = cnode->NumDof();

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    Epetra_SerialDenseMatrix jumpvec(dim,1);
    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // jump vector
    Epetra_SerialDenseMatrix forcevec(dim,1);
    for (int i=0;i<dim;i++)
      forcevec(i,0) = cnode->MoData().lm()[i];

    // evaluate jump in tangential direction
    Epetra_SerialDenseMatrix jumptan(dim,1);
    jumptan.Multiply('N','N',1,tanplane,jumpvec,0.0);

    // norm of tangential jumps for stick nodes
    if (cnode->Active()== true and cnode->FriData().Slip()==false)
    {
      for( int j=0;j<cnode->NumDof();++j)
        cnormtan+=jumptan(j,0)*jumptan(j,0);
    }
    else if (cnode->Active()== true and cnode->FriData().Slip()==true)
    {
      double jumptxi = 0;
      double jumpteta = 0;
      double forcen = 0;
      double forcetxi = 0;
      double forceteta = 0;

      for (int i=0;i<dim;i++)
      {
        jumptxi+=cnode->CoData().txi()[i]*cnode->FriData().jump()[i];
        jumpteta+=cnode->CoData().teta()[i]*cnode->FriData().jump()[i];

        forcen+=cnode->MoData().n()[i]*cnode->MoData().lm()[i];
        forcetxi+=cnode->CoData().txi()[i]*cnode->MoData().lm()[i];
        forceteta+=cnode->CoData().teta()[i]*cnode->MoData().lm()[i];
      }

      //cout << "FACTOR-Direction " << (jumptxi/jumpteta)/(forcetxi/forceteta) << endl;
      //cout << "FACTOR-Magnitude" << (frcoeff*forcen)/(sqrt(forcetxi*forcetxi+forceteta*forceteta))<< endl;
    }
  } // loop over slave nodes

  // get cnorm from all procs
  double sumcnormtanallprocs=0.0;
  Comm().SumAll(&cnormtan,&sumcnormtanallprocs,1);
  cnormtan=sumcnormtanallprocs;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                 popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegNormalForces(bool& localisincontact,
                                                 bool& localactivesetchange)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
  
  // penalty parameter
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->CoData().Getg();

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    //double modgap = cnode->CoData().Getg();
    //if (abs(modgap) < 1.0e-10) modgap=0.0;

    double kappa = cnode->CoData().Kappa();

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitely to constant
    // and corresponding derivatives to zero

    for( int j=0;j<dim;++j)
      cnode->MoData().lm()[j] = i*j;

    cnode->CoData().GetDerivZ().clear();

    continue;
#endif

    //********************************************************************
    // Decision on active /  inactive nodes (regularization)
    //
    // CASE 1: Penalty approach
    // A node is activated if its weighted gap is negative or deactivated
    // if its gap is equal zero or positive.
    // -> the regularization reads: lambda_n = kappa * pp * < -gap >
    //
    // CASE 2: Augmented Lagrange approach
    // A node is activated if its Lagrange multiplier, stemming from the
    // last Uzawa Lagrange multiplier AND the current regularization is
    // negative or deactivated if its LM is equal zero or positive.
    // -> the regularization reads: lambda_n = < lmuzawa_n - kappa * pp * gap >
    //
    // As the Uzawa Lagrange multipliers are zero in the penalty approach,
    // the two cases can formally be treted identically, see below.
    // We do not need an explicit separation of cases!
    //
    //********************************************************************

    // Activate/Deactivate node and notice any change
    if( (cnode->Active() == false) && (lmuzawan - kappa * pp * gap >= 0) )
    {
        cnode->Active() = true;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now active (";
        //for( int j=0; j<dim; j++)
        //  cout << " " << cnode->Dofs()[j] << " ";
        //cout << ") gap=" << gap << endl;
    }

    else if( (cnode->Active() == true) && (lmuzawan - kappa * pp * gap < 0) )
    {
        cnode->Active() = false;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now inactive, gap=" << gap << endl;
    }
    //********************************************************************

    // Compute derivZ-entries with the Macauley-Bracket
    // of course, this is only done for active constraints in order
    // for linearization and r.h.s to match!
    if( cnode->Active()==true )
    {

//      cout << "GID " << gid << endl;
//      cout << "LMUZAWAN " << lmuzawan << endl;
//      cout << "GAP " << gap << endl;

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      // contribution of derivative of normal
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();
      map<int,double>::iterator ncurr;

      for( int j=0;j<dim;++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
          cnode->AddDerivZValue(j, gcurr->first, - kappa * pp * (gcurr->second) * normal[j]);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, - kappa * pp * gap * ncurr->second);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, + lmuzawan * ncurr->second);
      }
    }

    // be sure to remove all LM-related stuff from inactive nodes
    else
    {
      // clear lagrange multipliers
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;

      // clear derivz
      cnode->CoData().GetDerivZ().clear();

    } // Macauley-Bracket
  } // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                 gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesPenalty()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);
       tanplane(0,2)=  -(n[0]*n[2]);
       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
       tanplane(1,2)=  -(n[1]*n[2]);

       tanplane(2,0)=  -(n[2]*n[0]);
       tanplane(2,1)=  -(n[2]*n[1]);
       tanplane(2,2)= 1-(n[2]*n[2]);
     }
     else if (dim==2)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);

       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
     }
     else
       dserror("Error in AssembleTangentForces: Unknown dimension.");

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // fill vector tractionold
    vector<double> tractionold(dim);
    for (int i=0;i<dim;i++)
      tractionold[i] = cnode->FriData().tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=tractionold[i]+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))
    {
      //cout << "Node " << gid << " is stick" << endl;
      cnode->FriData().Slip() = false;

      // in the stick case, traction is trailtraction
      for (int i=0;i<dim;i++)
        cnode->FriData().traction()[i]=trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
      //cout << "Node " << gid << " is slip" << endl;
      cnode->FriData().Slip() = true;

      // in the slip case, traction is evaluated with a return map algorithm
      for (int i=0;i<dim;i++)
        cnode->FriData().traction()[i]=maxtantrac/magnitude*trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap)+maxtantrac/magnitude*trailtraction[j];
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->FriData().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /**************************************** deriv(tanplane).jump  ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }
    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip()== true)
    {
      /******************** tanplane.deriv(jump).maxtantrac/magnidude ***/

      vector<map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** deriv(tanplane).jump.maxtantrac/magnitude ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->FriData().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j, gcurr->first, - frcoeff*kappa * ppnor * (gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
      for (int dim=0;dim<cnode->NumDof();++dim)
        temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

        traction+= tractionold[dimout];

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

        traction+=tractionold[dimout];

        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
        {
          int col = colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-colcurr->second*n[dim]*cnode->FriData().jump()[dim]*traction/magnitude*pptan*kappa;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

          traction += tractionold[dimout];

          for (int dim=0;dim<cnode->NumDof();++dim)
          {

          // loop over all entries of the current derivative map
          for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
          {
            int col = colcurr->first;

            double val =-colcurr->second*n[dimout]*cnode->FriData().jump()[dim]*traction/magnitude*pptan*kappa;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
               double val1 = val*temp[dimrow];
               cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces (Aug. Lagr.)    gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesAugmented()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // Evaluate trailtraction
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=lmuzawatan(i,0)+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))    {
      //cout << "Node " << gid << " is stick" << endl;
      cnode->FriData().Slip() = false;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j];
    }
    else
    {
      //cout << "Node " << gid << " is slip" << endl;
      cnode->FriData().Slip() = true;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j]*maxtantrac/magnitude;
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->FriData().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*(cnode->FriData().jump()[dim]);
            val = val - (colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*(cnode->FriData().jump()[dim]);
            val = val-(colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }

    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip()== true)
    {
      /***************************************** tanplane.deriv(jump) ***/
      vector<map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      map<int,double>::iterator colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      vector<map<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->FriData().jump()[dim];
            val = (val - (colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = colcurr->first;
            double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim];
            val = (val-(colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->CoData().GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j,gcurr->first,- frcoeff*kappa*ppnor*(gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( colcurr = (derivn[j]).begin(); colcurr != (derivn[j]).end(); ++colcurr )
        {
          for( int k=0;k<cnode->NumDof();++k)
          {
            double val = frcoeff*(colcurr->second)*lmuzawa(j,0)*trailtraction[k]/magnitude;
            cnode->AddDerivZValue(k,colcurr->first,val);
          }
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
        for (int dim=0;dim<cnode->NumDof();++dim)
          temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        // loop over all entries of the current derivative map
        for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
        {
          int col = colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-colcurr->second*n[dim]*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*pptan*kappa)*traction/magnitude;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =-colcurr->second*n[dimout]*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*pptan*kappa)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble derivatives of lagrange multipliers              popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinZ(LINALG::SparseMatrix& linzglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinZ: Node ownership inconsistency!");

    // derivz is the vector<map> we want to assemble
    vector<map<int,double> >& derivz = cnode->CoData().GetDerivZ();

    if ( (int) derivz.size()>0 )
    {
      int rowsize = cnode->NumDof();
      int colsize = (int) derivz[0].size();

      // consistency check
      for (int j=0; j<rowsize-1; ++j)
        if ((int)derivz[j].size() != (int)derivz[j+1].size())
          dserror("ERROR: AssembleLinZ: Column dim. of nodal derivz-map is inconsistent!");

      map<int,double>::iterator colcurr;

      // loop over dofs
      for ( int k=0; k<rowsize; ++k )
      {
        int row = cnode->Dofs()[k]; // row index equals global dof index of this #i node's dof k
        int l = 0;

        // loop over all directional derivative entries using the map iterator
        for( colcurr = derivz[k].begin(); colcurr != derivz[k].end(); ++colcurr )
        {
          int col = colcurr->first; // col index equals global id of directional derivative component ,l
          double val = colcurr->second;
          linzglobal.Assemble(val, row, col);
          l++;
        }

        if( l != colsize )
          dserror("ERROR: AssembleLinZ: l = %i but colsize = %i",k,colsize);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices with nodal normals / tangents           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNT(LINALG::SparseMatrix& nglobal,
                                     LINALG::SparseMatrix& tglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleNT: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(1);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(1);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(i);
      lmrowownerT[0] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->MoData().n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(1,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->CoData().txi()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }

    else if (Dim()==3)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(2);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(2);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(2*i);
      lmrowT[1] = activet_->GID(2*i+1);
      lmrowownerT[0] = cnode->Owner();
      lmrowownerT[1] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->MoData().n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(2,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->CoData().txi()[j];
        Tnode(1,j) = cnode->CoData().teta()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }
    else
      dserror("ERROR: Dim() must be either 2D or 3D");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives           popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    map<int,double>& dgmap = cnode->CoData().GetDerivG();
    map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;
      //cout << "Assemble S: " << row << " " << col << " " << val << endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix P containing tangent derivatives          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleP(LINALG::SparseMatrix& pglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleP: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      vector<map<int,double> >& dtmap = cnode->CoData().GetDerivTxi();
      map<int,double>::iterator colcurr;
      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();
      int row = activet_->GID(i);

      if (mapsize==3) mapsize=2;

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleP: Column dim. of nodal DerivT-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions)
      for (int j=0;j<mapsize;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << row << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,row,col);
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsize);
      }
    }
    else if (Dim()==3)
    {
      // prepare assembly
      vector<map<int,double> >& dtximap = cnode->CoData().GetDerivTxi();
      vector<map<int,double> >& dtetamap = cnode->CoData().GetDerivTeta();
      map<int,double>::iterator colcurr;
      int colsizexi = (int)dtximap[0].size();
      int colsizeeta = (int)dtetamap[0].size();
      int mapsizexi = (int)dtximap.size();
      int mapsizeeta = (int)dtetamap.size();
      int rowxi = activet_->GID(2*i);
      int roweta = activet_->GID(2*i+1);

      for (int j=0;j<mapsizexi-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTXi-map is inconsistent!");

      for (int j=0;j<mapsizeeta-1;++j)
        if ((int)dtetamap[j].size() != (int)dtetamap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTEta-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions) for TXi
      for (int j=0;j<mapsizexi;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << rowxi << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,rowxi,col);
          ++k;
        }

        if (k!=colsizexi)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizexi);
      }

      // loop over all derivative maps (=dimensions) for TEta
      for (int j=0;j<mapsizeeta;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->MoData().lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << roweta << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,roweta,col);
          ++k;
        }

        if (k!=colsizeeta)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizeeta);
      }
    }
    else
      dserror("ERROR: Dim() must be either 2 or 3!");

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices LinDM containing fc derivatives         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                                       LINALG::SparseMatrix& linmglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
  
  /**********************************************************************/
  // NEW VERSION (09/2010): No more communication, thanks to FE_MATRIX!
  /**********************************************************************/
  // we have: D_jk,c with j = Lagrange multiplier slave dof
  //                 with k = Displacement slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinD)_kc = D_jk,c * z_j
  /**********************************************************************/
  // we have: M_jl,c with j = Lagrange multiplier slave dof
  //                 with l = Displacement master dof
  //                 with c = Displacement slave or master dof
  // we compute (LinM)_lc = M_jl,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    int dim = cnode->NumDof();
    
    // Mortar matrix D and M derivatives
    map<int,map<int,double> >& dderiv = cnode->CoData().GetDerivD();
    map<int,map<int,double> >& mderiv = cnode->CoData().GetDerivM();

    // current Lagrange multipliers
    double* lm = cnode->MoData().lm();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    int mastersize = (int)mderiv.size();
    map<int,map<int,double> >::iterator scurr = dderiv.begin();
    map<int,map<int,double> >::iterator mcurr = mderiv.begin();
    
    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
      CoNode* csnode = static_cast<CoNode*>(snode);

      // Mortar matrix D derivatives
      map<int,double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());
 
      // inner product D_{jk,c} * z_j for index j
      for (int prodj=0;prodj<dim;++prodj)
      {
        int row = csnode->Dofs()[prodj];
        map<int,double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = scolcurr->first;
          double val = lm[prodj] * (scolcurr->second);
          ++scolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          //cout << "Assemble LinD: " << row << " " << col << " " << val << endl;
          if (abs(val)>1.0e-12) lindglobal.FEAssemble(val,row,col);
        }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr!=thisdderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
      }
    }

    // check for completeness of DerivD-Slave-iteration
    if (scurr!=dderiv.end())
      dserror("ERROR: AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    /******************************** Finished with LinDMatrix **********/
    
        
    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
      CoNode* cmnode = static_cast<CoNode*>(mnode);
      
      // Mortar matrix M derivatives
      map<int,double>&thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());
 
      // inner product M_{jl,c} * z_j for index j
      for (int prodj=0;prodj<dim;++prodj)
      {
        int row = cmnode->Dofs()[prodj];
        map<int,double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = mcolcurr->first;
          double val = lm[prodj] * (mcolcurr->second);
          ++mcolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          //cout << "Assemble LinM: " << row << " " << col << " " << val << endl;
          if (abs(val)>1.0e-12) linmglobal.FEAssemble(-val,row,col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
      }
    }

    // check for completeness of DerivM-Master-iteration
    if (mcurr!=mderiv.end())
      dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
    /******************************** Finished with LinMMatrix **********/
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble normal weighted gap                              popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleG(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDMG: Node ownership inconsistency!");

    /**************************************************** g-vector ******/
    if (cnode->CoData().Getg()!=0.0)
    {
      double gap = cnode->CoData().Getg();

      // cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj()
      //      << " IsActive: " << cnode->Active() << " Gap: " << gap << endl;

      // check if this inactive node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      //******************************************************************
      // TODO: This is only necessary for dual LM shape functions and for
      // quadratic standard LM shape functions! By the way, it makes the
      // method slightly inconsistent (e.g. patch tests with slave side
      // being wider than master side). However, we are able to solve many
      // problems with this little trick. But not all problems, e.g.
      // dropping edge problems would still fail!!! To solve this dilemma,
      // we need a clever modification of the LM shape functions such that
      // their definition is compressed to only the "projecting" element part.
      // Once we have this, the following trick can (and should) also be
      // removed in order to make the method consistent again! (10/2010)
      //******************************************************************
      if (!cnode->HasProj() && !cnode->Active())
      {
        gap = 1.0e12;
        cnode->CoData().Getg()=gap;
      }

      Epetra_SerialDenseVector gnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      gnode(0) = gap;
      lm[0] = cnode->Id();
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(gglobal,gnode,lm,lmowner);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble wear                                         gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleWear(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* frinode = static_cast<FriNode*>(node);

    if (frinode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWear: Node ownership inconsistency!");

    /**************************************************** w-vector ******/
    double wear = frinode->FriData().Wear();
      
    Epetra_SerialDenseVector wnode(1);
    vector<int> lm(1);
    vector<int> lmowner(1);

    wnode(0) = wear;
    lm[0] = frinode->Id();
    lmowner[0] = frinode->Owner();

    LINALG::Assemble(gglobal,wnode,lm,lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                          LINALG::SparseMatrix& linstickDISglobal,
                                          Epetra_Vector& linstickRHSglobal)
{
  // FIXGIT: Assemble LinStick is containing a matrix for the de-
  // rivatives of the Lagrange multipliers. This is according to Hüeber.
  // Because of worse convergence, this is not implemented, but the
  // code is commented after the algorithm.

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // create map of stick nodes
  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  // loop over all stick nodes of the interface
  for (int i=0;i<sticknodes->NumMyElements();++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    // prepare assembly, get information from node
    vector<map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
    vector<map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

    for (int j=0;j<Dim()-1;++j)
      if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

    if (Dim()==3)
    {
      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
    }

    // more information from node
    double* txi = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* jump = cnode->FriData().jump();
        
    // iterator for maps
    map<int,double>::iterator colcurr;

    // row number of entries
    vector<int> row (Dim()-1);
    if (Dim()==2)
    {
      row[0] = stickt->GID(i);
    }
    else if (Dim()==3)
    {
      row[0] = stickt->GID(2*i);
      row[1] = stickt->GID(2*i)+1;
    }
    else
      dserror("ERROR: AssemblelinStick: Dimension not correct");

    // evaluation of specific components of entries to assemble
    double jumptxi=0;
    double jumpteta=0;
    for (int dim = 0;dim < Dim();dim++)
    {
      jumptxi += txi[dim]*jump[dim];
      jumpteta += teta[dim]*jump[dim];
    }

    // check for dimensions
    if(Dim()==2 and (jumpteta != 0.0))
      dserror ("ERROR: AssembleLinStick: jumpteta must be zero in 2D");

    // Entries on right hand side
    /************************************************ (-utxi, -uteta) ***/
    Epetra_SerialDenseVector rhsnode(Dim()-1);
    vector<int> lm(Dim()-1);
    vector<int> lmowner(Dim()-1);

    rhsnode(0) = -jumptxi;
    lm[0] = cnode->Dofs()[1];
    lmowner[0] = cnode->Owner();

    if (Dim()==3)
    {
      rhsnode(1) = -jumpteta;
      lm[1] = cnode->Dofs()[2];
      lmowner[1] = cnode->Owner();
    }

    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
    
    // Entries from differentiation with respect to displacements
    /*** 1 ************************************** tangent.deriv(jump) ***/

    // get linearization of jump vector
    vector<map<int,double> > derivjump = cnode->FriData().GetDerivJump();
    
    if (derivjump.size()<1)
      dserror ("AssembleLinStick: Derivative of jump is not exiting!");

    // loop over dimensions
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
      // loop over all entries of the current derivative map (jump)
      for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
      {
        int col = colcurr->first;
        double valtxi = txi[dim]*colcurr->second;

        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

        if(Dim()==3)
        {
          double valteta = teta[dim]*colcurr->second;
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }
    }
    
    /*** 2 ************************************** deriv(tangent).jump ***/
    // loop over dimensions
    for (int j=0;j<Dim();++j)
    {
      // loop over all entries of the current derivative map (txi)
      for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
      {
        int col = colcurr->first;
        double val = jump[j]*colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
      }

      if(Dim()==3)
      {
        // loop over all entries of the current derivative map (teta)
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = jump[j]*colcurr->second;

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
        }
      }
    }
  }

//  have a look to the beginning of the function
//  // only for coulomb friction
//  string ftype   = IParams().get<string>("friction type","none");
//  double frcoeff = IParams().get<double>("friction coefficient",0.0);
//  double cn = IParams().get<double>("semismooth cn",0.0);
//  double ct = IParams().get<double>("semismooth ct",0.0);
//
//  if (ftype == "tresca")
//    dserror ("Error: AssemblelinStick: complementary function according"
//             " to Hueber only available for Coulomb friction");
//
//  // get out of here if not participating in interface
//  if (!lComm())
//    return;
//
//  // create map of stick nodes
//  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
//  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);
//
//  // nothing to do if no stick nodes
//  if (sticknodes->NumMyElements()==0)
//    return;
//
//  // not yet implemented for 3D
//    if (Dim()==3)
//      dserror("ERROR: AssembleLinStick: 3D not yet implemented");
//
//  // loop over all stick nodes of the interface
//  for (int i=0;i<sticknodes->NumMyElements();++i)
//  {
//    int gid = sticknodes->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
//    FriNode* cnode = static_cast<FriNode*>(node);
//
//    if (cnode->Owner() != Comm().MyPID())
//      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");
//
//    // prepare assembly
//    vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();
//    map<int,double>::iterator colcurr;
//
//    // calculate DerivT from DerivN
//    // only for 2D so far, in this case calculation is very easy
//    // dty =  dnx
//    // dtx = -dny
//
//    vector <map<int,double> > dtmap(Dim());
//
//    for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
//      dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));
//
//    for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
//      dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));
//
//    int colsize = (int)dtmap[0].size();
//    int mapsize = (int)dtmap.size();
//    int row = stickt->GID(i);
//    double* xi = cnode->xspatial();
//    double* txi = cnode->CoData().txi();
//    double* jump = cnode->FriData().jump();
//    double& wgap = cnode->CoData().Getg();
//
//    double utan = 0;
//    double nz = 0;
//
//    for (int dim = 0;dim < Dim();dim++)
//    {
//      utan += txi[dim]*jump[dim];
//      nz += cnode->MoData().n()[dim] * cnode->MoData().lm()[dim];
//    }
//
//    // initialization of nz if nz = wgap = 0
//    if(nz==wgap and wgap==0)
//    {
//      nz=1;
//      cout << "Warning: Initialization of nz to 1" << endl;
//    }
//
//    for (int j=0;j<mapsize-1;++j)
//      if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
//        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivT-map is inconsistent!");
//
//    // Entries on right hand side
//    /**************************************** frcoeff*cn*wgap*ct*utan ***/
//
//    Epetra_SerialDenseVector rhsnode(1);
//    vector<int> lm(1);
//    vector<int> lmowner(1);
//
//    rhsnode(0) = -frcoeff*cn*wgap*ct*utan;
//    lm[0] = cnode->Dofs()[1];
//    lmowner[0] = cnode->Owner();
//
//    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
//
//    // Entries from differentiation with respect to lagrange multipliers
//    /*******************/
//
//    // loop over the dimension
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*cnode->MoData().n()[dim]*ct*utan;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickLMglobal.Assemble(val,row,col);
//    }
//
//    // Entries from differentiation with respect to displacements
//    /************************************************** -tan.(D-Dn-1) ***/
//
//    // we need the nodal entries of the D-matrix and the old one
//    double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
//    double Dold= (cnode->FriData().GetDOld()[0])[cnode->Dofs()[0]];
//
//    // loop over all derivative maps (=dimensions)
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*txi[dim]*(D-Dold);
//     //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//     // do not assemble zeros into matrix
//     if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /*************************************************** tan.(M-Mn-1) ***/
//
//    // we need the nodal entries of the M-matrix and the old one
//    vector<map<int,double> > mmap = cnode->MoData().GetM();
//    vector<map<int,double> > mmapold = cnode->FriData().GetMOld();
//
//    // create a set of nodes including nodes according to M entries
//    // from current and previous time step
//    set <int> mnodes;
//
//    for (colcurr=mmap[0].begin(); colcurr!=mmap[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    if(mmapold.size()<1)
//    {
//      cout << "GID " << gid << endl;
//      dserror("vector too small");
//    }
//
//    for (colcurr=mmapold[0].begin(); colcurr!=mmapold[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    set<int>::iterator mcurr;
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int dim=0;dim<cnode->NumDof();++dim)
//      {
//        int col = cmnode->Dofs()[dim];
//        double val = -frcoeff*(nz-cn*wgap)*ct*txi[dim]*(mik-mikold);
//        //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//       // do not assemble zeros into matrix
//       if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /***************************************  -DerivT.(D-Dn-1).xs  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*(D-Dold)*xi[j]*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//
//    /***************************************  -DerivT.(M-Mn-1).xm  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      double* mxi = cmnode->xspatial();
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int j=0;j<mapsize;++j)
//      {
//        int k=0;
//
//        // loop over all entries of the current derivative map
//        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//        {
//          int col = colcurr->first;
//          double val = -frcoeff*(nz-cn*wgap)*ct*(mik-mikold)*mxi[j]*colcurr->second;
//          // do not assemble zeros into matrix
//          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//          ++k;
//        }
//
//        if (k!=colsize)
//          dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//      }
//    }
//
//    /**********************************************  -T.DerivD.x  *******/
//
//    // we need the dot product n*x of this node
//    double tdotx = 0.0;
//    for (int dim=0;dim<cnode->NumDof();++dim)
//      tdotx += txi[dim]*xi[dim];
//
//    // prepare assembly
//    map<int,double>& ddmap = cnode->CoData().GetDerivD();
//
//    // loop over all entries of the current derivative map
//    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*tdotx*colcurr->second;
//
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /***********************************************   -T.DerivM.x ******/
//
//    // we need the Lin(M-matrix) entries of this node
//    map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
//    map<int,map<int,double> >::iterator dmcurr;
//
//    // loop over all master nodes in the DerivM-map of the active slave node
//    for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
//    {
//      int gid = dmcurr->first;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      FriNode* cmnode = static_cast<FriNode*>(mnode);
//      double* mxi = cmnode->xspatial();
//
//      // we need the dot product ns*xm of this node pair
//      double tdotx = 0.0;
//      for (int dim=0;dim<cnode->NumDof();++dim)
//        tdotx += txi[dim]*mxi[dim];
//
//      // compute matrix entry of the current active node / master node pair
//      map<int,double>& thisdmmap = cnode->GetDerivM(gid);
//
//      // loop over all entries of the current derivative map
//      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*tdotx*colcurr->second;
//
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /**************************************-frcoeff*cn*ct*utan*DerivG ***/
//
//    // prepare assembly
//    map<int,double>& dgmap = cnode->CoData().GetDerivG();
//
//    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = +frcoeff*cn*ct*utan*colcurr->second;
//      //cout << "Assemble LinStick: " << row << " " << col << " " << val << endl;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /************************************** -frcoeff*DerivN*z+ct*utan ***/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*cnode->MoData().lm()[j]*ct*utan*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//  }
  return;
}

/*----------------------------------------------------------------------*
|  Assemble matrix LinSlip with tangential+D+M derivatives    mgit 02/09|
*----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
                                         LINALG::SparseMatrix& linslipDISglobal,
                                         Epetra_Vector& linslipRHSglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
  double frbound = IParams().get<double>("FRBOUND");
  double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  bool fulllin = DRT::INPUT::IntegralValue<int>(IParams(),"FULL_LINEARIZATION");

  // Coulomb Friction
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
#ifdef CONTACTCOMPHUEBER

    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();
      vector<map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      vector<map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       if (Dim()==3)
       {
         for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
       }

      // more information from node
      double* jump = cnode->FriData().jump();
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* z = cnode->MoData().lm();
      double& wgap = cnode->CoData().Getg();

      // iterator for maps
      map<int,double>::iterator colcurr;

      // row number of entries
      vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = slipt_->GID(2*i);
        row[1] = slipt_->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinSlip: Dimension not correct");

      // boolean variable if flag "CONTACTFRICTIONLESSFIRST" AND
      // ActiveOld = true
      bool friclessandfirst = false;

      // evaluation of specific components of entries to assemble
      double znor = 0;
      double ztxi = 0;
      double zteta = 0;
      double jumptxi = 0;
      double jumpteta = 0;
      double euclidean = 0;
      for (int i=0;i<Dim();i++)
      {
        znor += n[i]*z[i];
        ztxi += txi[i]*z[i];
        zteta += teta[i]*z[i];
        jumptxi += txi[i]*jump[i];
        jumpteta += teta[i]*jump[i];
      }

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      vector<double> sum1 (Dim()-1,0);
      sum1[0] = ztxi+ct*jumptxi;
      if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
      if (Dim()==2) euclidean = abs(sum1[0]);
      if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

      // check of dimensions
      if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      // check of euclidean norm
      if (euclidean==0.0)
        dserror ("ERROR: AssemblelinSlip: Euclidean norm is zero");

#ifdef CONTACTFRICTIONLESSFIRST

      // in the case of frictionless contact for nodes just coming into
      // contact, the frictionless contact condition is applied.
      if (cnode->FriData().ActiveOld()==false)
      {
        friclessandfirst=true;
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];
          double valteta = 0;
          if (Dim()==3) valteta = teta[dim];

          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);

        }
        if(fulllin)
        {
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (colcurr->second)*z[dim];
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            }

            if(Dim()==3)
            {
              for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
              {
                int col = colcurr->first;
                double valteta = (colcurr->second)*z[dim];
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
        }
      }
#endif

      // this is not evaluated if "FRICTIONLESSFIRST" is flaged on AND the node
      // is just coming into contact
      if(friclessandfirst==false)
      {
        /******************************************************************/
        // calculation of matrix entries of the linearized slip condition
        /******************************************************************/
        // 1) Entries from differentiation with respect to lagrange multipliers
        // 2) Entries on right hand side
        // 3) Entries from differentiation with respect to displacements

        // 1) Entries from differentiation with respect to lagrange multipliers
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          double valtxi = 0;
          int col = cnode->Dofs()[dim];
          double valtxi0 = euclidean*txi[dim];
          double valtxi1 = ((ztxi+ct*jumptxi)/euclidean*ztxi)*txi[dim];
          double valtxi3 = (zteta+ct*jumpteta)/euclidean*ztxi*teta[dim];
          double valtxi2 = -frcoeff*(znor-cn*wgap)*txi[dim]-frcoeff*(ztxi+ct*jumptxi)*n[dim];
          valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

          double valteta = 0;
          if (Dim()==3)
          {
            double valteta0 = euclidean*teta[dim];
            double valteta1 = ((ztxi+ct*jumptxi)/euclidean*zteta)*txi[dim];
            double valteta3 = (zteta+ct*jumpteta)/euclidean*zteta*teta[dim];
            double valteta2 = -frcoeff*(znor-cn*wgap)*teta[dim]-frcoeff*(zteta+ct*jumpteta)*n[dim];
            valteta = valteta0 + valteta1 + valteta2 + valteta3;
          }

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
        }

        // 2) Entries on right hand side
        /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

        double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
        double valuetxi2 = +euclidean*ztxi;
        double valuetxi3 = (ztxi+ct*jumptxi)/euclidean*ztxi*ztxi;
        double valuetxi4 = (zteta+ct*jumpteta)/euclidean*zteta*ztxi;
        double valuetxi5 = -(frcoeff*(znor-cn*wgap))*ztxi-(frcoeff*znor)*(ztxi+ct*jumptxi);

        Epetra_SerialDenseVector rhsnode(Dim()-1);
        vector<int> lm(Dim()-1);
        vector<int> lmowner(Dim()-1);

        rhsnode(0) = (valuetxi1+valuetxi2+valuetxi3+valuetxi4+valuetxi5);
        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        if(Dim()==3)
        {
          double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
          double valueteta2 = +euclidean*zteta;
          double valueteta3 = (ztxi+ct*jumptxi)/euclidean*ztxi*zteta;
          double valueteta4 = (zteta+ct*jumpteta)/euclidean*zteta*zteta;
          double valueteta5 = -(frcoeff*(znor-cn*wgap))*zteta-(frcoeff*znor)*(zteta+ct*jumpteta);

          rhsnode(1) = (valueteta1+valueteta2+valueteta3+valueteta4+valueteta5);
          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }

        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

        // 3) Entries from differentiation with respect to displacements
        /*** 01  ********* -Deriv(euclidean).ct.tangent.deriv(u)*ztan ***/

        // get linearization of jump vector
        vector<map<int,double> > derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
         {
            int col = colcurr->first;
            
            double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*colcurr->second*ztxi;
            double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*colcurr->second*zteta;
            double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*colcurr->second*ztxi;
            double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*colcurr->second*zteta;

            // do not assemble zeros into matrix
            if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
            if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
            if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
            if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
          }
        }

        /*** 02 ***************** frcoeff*znor*ct*tangent.deriv(jump) ***/
        
        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
         {
            int col = colcurr->first;
            
            //cout << "val " << colcurr->second << endl;

            double valtxi = (-1)*(frcoeff*(znor-cn*wgap))*ct*txi[dim]*colcurr->second;
            double valteta = (-1)*(frcoeff*(znor-cn*wgap))*ct*teta[dim]*colcurr->second;

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
             if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }
          
        }

        // remaining terms only in case of full linearization
        if(fulllin)
        {
          /*** 1 ********************************* euclidean.deriv(T).z ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean*(colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if (Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = euclidean*(colcurr->second)*z[j];

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*zteta;

             // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (Dim()==3)
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // 3D loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*zteta;

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/

          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

                // do not assemble zeros into matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map
            for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (-1)*(ztxi+ct*jumptxi)*frcoeff*(colcurr->second)*z[j];
              double valteta = (-1)*(zteta+ct*jumpteta)*frcoeff*(colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

          /*** 7 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
          // prepare assembly
          map<int,double>& dgmap = cnode->CoData().GetDerivG();

          // loop over all entries of the current derivative map
          for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
            double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }
        } // if fullin
      } // if (frictionlessandfirst == false)
    } // loop over all slip nodes of the interface
#else
    
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->FriData().jump();
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();
      
      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double znor    = n[0]*z[0] + n[1]*z[1];
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frcoff*znor).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frcoeff*znor)*txi[dim]-frcoeff*(ztan+ct*jumptan)*n[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1 = -(abs(ztan+ct*jumptan))*ztan+(frcoeff*znor)*(ztan+ct*jumptan);
      double value2 = -(frcoeff*znor)*ztan-(frcoeff*znor)*(ztan+ct*jumptan);
      double value3 = +sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      rhsnode(0) = (value1+value2+value3);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->FriData().GetDOld()[0])[cnode->Dofs()[0]];

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif


       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->FriData().GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->FriData().GetMNodes();
      set <int> mnodesold = cnode->FriData().GetMNodesOld();

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /********************************* frcoeff*znor*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (frcoeff*znor)*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -frcoeff*znor*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = (frcoeff*znor)*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

          // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /*********************************** -(frcoeff*znor).DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** (frcoeff*znor).ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** +(frcoeff*znor).ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /***************************** -(frcoeff*znor).ct.T.DerivM.x ***/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           FriNode* cmnode = static_cast<FriNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /***************************** -frcoeff*DerivN.z(ztan+ct*utan) ***/

       // loop over all derivative maps (=dimensions)
       for (int j=0;j<mapsize;++j)
       {
         int k=0;

         // loop over all entries of the current derivative map
         for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(ztan+ct*jumptan)*frcoeff*(colcurr->second)*z[j];
           //cout << "10 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into s matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
           ++k;
         }

         if (k!=colsize)
           dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }
      } // if fullin
    }
#endif
  } // Coulomb friction

  // Tresca Friction
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->FriData().jump();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frbound)*txi[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1= -(abs(ztan+ct*jumptan))*ztan+frbound*(ztan+ct*jumptan);
      double value2= -frbound*ztan+sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);
      rhsnode(0) = (value1+value2);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->FriData().GetDOld()[0])[cnode->Dofs()[0]];

      if (abs(Dold)<0.0001)
        dserror ("Error:No entry for Dold");

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->MoData().GetM();
      vector<map<int,double> > mmapold = cnode->FriData().GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->FriData().GetMNodes();
      set <int> mnodesold = cnode->FriData().GetMNodesOld();

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /************************************** frbound*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = frbound*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

        // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cmnode = static_cast<FriNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = frbound*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = 0;
#endif
          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
            if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*frbound*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /********************************  -frbound.ct.T.DerivM.x ******/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           FriNode* cmnode = static_cast<FriNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*frbound*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      } // if fullin
    }
  }// Tresca friction

  return;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::BuildActiveSet(bool init)
{
  // define local variables
  vector<int> mynodegids(0);
  vector<int> mydofgids(0);
  vector<int> myslipnodegids(0);
  vector<int> myslipdofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // check if node is initially active
      if (cnode->IsInitActive())
      {
        cnode->Active()=true;
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is initially in slip state
      if (friction_)
      {
        // do nothing: we always assume STICK at t=0
      }
    }

    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->Active())
      {
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is in slip state
      if (friction_)
      {
        if (static_cast<FriNode*>(cnode)->FriData().Slip())
        {
          myslipnodegids.push_back(cnode->Id());

          for (int j=0;j<numdof;++j)
            myslipdofgids.push_back(cnode->Dofs()[j]);
        }
      }
    }
  }

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(-1,(int)mynodegids.size(),&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(-1,(int)mydofgids.size(),&mydofgids[0],0,Comm()));

  if (friction_)
  {
    // create slip node map and slip dof map
    slipnodes_ = rcp(new Epetra_Map(-1,(int)myslipnodegids.size(),&myslipnodegids[0],0,Comm()));
    slipdofs_  = rcp(new Epetra_Map(-1,(int)myslipdofgids.size(),&myslipdofgids[0],0,Comm()));
  }

  // split active dofs and slip dofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs         popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::SplitActiveDofs()
{
  // get out of here if active set is empty
  if (activenodes_==null)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements()==0)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countN=0;
  int countT=0;
  vector<int> myNgids(activenodes_->NumMyElements());
  vector<int> myTgids((Dim()-1)*activenodes_->NumMyElements());

  // dimension check
  double dimcheck =(activedofs_->NumGlobalElements())/(activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j=1;j<numdof;++j)
    {
      myTgids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  Comm().SumAll(&countN,&gcountN,1);
  Comm().SumAll(&countT,&gcountT,1);

  // check global dimensions
  if ((gcountN+gcountT)!=activedofs_->NumGlobalElements())
    dserror("ERROR: SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));
  activet_ = rcp(new Epetra_Map(gcountT,countT,&myTgids[0],0,Comm()));
  
  // *******************************************************************
  // FRICTION - EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if there is no friction
  if(friction_==false)
    return true; 
  
  // get out of here if slip set is empty
  if (slipnodes_==null)
  {
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  if (slipnodes_->NumGlobalElements()==0)
  {
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countslipT=0;
  vector<int> myslipTgids((Dim()-1)*slipnodes_->NumMyElements());

  // dimension check
  dimcheck =(slipdofs_->NumGlobalElements())/(slipnodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add dofs to slipTmap
    for (int j=1;j<numdof;++j)
    {
      myslipTgids[countslipT] = cnode->Dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  Comm().SumAll(&countslipT,&gcountslipT,1);

  // create Tslipmap objects
  slipt_   = rcp(new Epetra_Map(gcountslipT,countslipT,&myslipTgids[0],0,Comm()));

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix A                                      gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleA(LINALG::SparseMatrix& aglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* frinode = static_cast<FriNode*>(node);

    if (frinode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleA: Node ownership inconsistency!");

    /**************************************************** A-matrix ******/
    if ((frinode->FriData().GetA()).size()>0)
    {
      vector<map<int,double> > amap = frinode->FriData().GetA();
      int rowsize = frinode->NumDof();
      int colsize = (int)amap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)amap[j].size() != (int)amap[j+1].size())
          dserror("ERROR: AssembleA: Column dim. of nodal A-map is inconsistent!");

      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = frinode->Dofs()[j];
        int k = 0;

        for (colcurr=amap[j].begin();colcurr!=amap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;

          // do the assembly into global A matrix
          // create the A matrix, do not assemble zeros
          aglobal.Assemble(val, row, col);
    
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleA: k = %i but colsize = %i",k,colsize);
      }
    }
  }
  return;
}


#endif  // #ifdef CCADISCRET
