/*!----------------------------------------------------------------------
\file drt_discret_evaluate.cpp
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_timecurve.H"
#include "drt_function.H"
#include "linalg_utils.H"
#include "linalg_sparsematrix.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                        Teuchos::ParameterList&              params,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                        Teuchos::RCP<Epetra_Vector>          systemvector1,
                        Teuchos::RCP<Epetra_Vector>          systemvector2,
                        Teuchos::RCP<Epetra_Vector>          systemvector3)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate");

  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblemat2 = systemmatrix2!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // call the element's register class preevaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate PreEvaluate");
    map<int,RefCountPtr<ElementRegister> >::iterator curr;
    for (curr=elementregister_.begin(); curr != elementregister_.end(); ++curr)
      curr->second->PreEvaluate(*this,params,systemmatrix1,systemmatrix2,
                                systemvector1,systemvector2,systemvector3);
  }

#ifdef THROWELEMENTERRORS
  EnterElementLoop();
#endif

  vector<int> lm;
  vector<int> lmowner;

  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate LocationVector");
    // get element location vector, dirichlet flags and ownerships
    lm.clear();
    lmowner.clear();
    actele->LocationVector(*this,lm,lmowner);
    }

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate Resize");

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    if (assemblemat1)
    {
      if (elematrix1.M()!=eledim or elematrix1.N()!=eledim)
        elematrix1.Shape(eledim,eledim);
      else
        memset(elematrix1.A(),0,eledim*eledim*sizeof(double));
    }
    if (assemblemat2)
    {
      if (elematrix2.M()!=eledim or elematrix2.N()!=eledim)
        elematrix2.Shape(eledim,eledim);
      else
        memset(elematrix2.A(),0,eledim*eledim*sizeof(double));
    }
    if (assemblevec1)
    {
      if (elevector1.Length()!=eledim)
        elevector1.Size(eledim);
      else
        memset(elevector1.Values(),0,eledim*sizeof(double));
    }
    if (assemblevec2)
    {
      if (elevector2.Length()!=eledim)
        elevector2.Size(eledim);
      else
        memset(elevector2.Values(),0,eledim*sizeof(double));
    }
    if (assemblevec3)
    {
      if (elevector3.Length()!=eledim)
        elevector3.Size(eledim);
      else
        memset(elevector3.Values(),0,eledim*sizeof(double));
    }
    }

#ifdef THROWELEMENTERRORS
    try {
#endif

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");
    // call the element evaluate method
    int err = actele->Evaluate(params,*this,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");
      int eid = actele->Id();
      if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,lm,lmowner);
      if (assemblemat2) systemmatrix2->Assemble(eid,elematrix2,lm,lmowner);
      if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
      if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
      if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);
    }

#ifdef THROWELEMENTERRORS
    }
    catch (const std::string& err)
    {
      ElementError(actele->Id(),err);
    }
#endif

  } // for (int i=0; i<numcolele; ++i)

#ifdef THROWELEMENTERRORS
  ExitElementLoop();
#endif

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        u.kue 01/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList&              params,
                                   Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                                   Teuchos::RCP<Epetra_Vector>          systemvector)
{
  Evaluate(params, systemmatrix, Teuchos::null, systemvector);
}



/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList& params, Epetra_Vector& systemvector)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"PointNeumann") continue;
    DRT::Condition& cond = *(fool->second);
    const vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    const vector<int>*    curve  = cond.Get<vector<int> >("curve");
    const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
    const vector<double>* val    = cond.Get<vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
    //cout << "Neumann load curve factor on point " << curvefac << endl;
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      vector<int> dofs = Dof(actnode);
      const unsigned numdf = dofs.size();
      for (unsigned j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= curvefac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }
  //--------------------------------------------------------
  // loop through line/surface/volume Neumann BCs and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
    if (fool->first == (string)"LineNeumann" ||
        fool->first == (string)"SurfaceNeumann" ||
        fool->first == (string)"VolumeNeumann"
       )
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*this,lm,lmowner);
        elevector.Size((int)lm.size());
        curr->second->EvaluateNeumann(params,*this,cond,lm,elevector);
        //cout << "Neumann load curve factor on element " << elevector << endl;
        LINALG::Assemble(systemvector,elevector,lm,lmowner);
      }
    }
  return;
}

/*----------------------------------------------------------------------*/
/*!
\brief Determine Dirichlet condition at given time and apply its
       values to a system vector

\param cond            The condition object
\param dis             The discretisation
\param usetime
\param time            Evaluation time
\param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
\param systemvectord   First time derivative of DBCs
\param systemvectordd  Second time derivative of DBCs
\param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this component remains untouched
\param dbcgids         Map containing DOFs subjected to Dirichlet boundary conditions
\date 02/08
*/
static void DoDirichletCondition(DRT::Condition&             cond,
                                 DRT::Discretization&        dis,
                                 const bool                  usetime,
                                 const double                time,
                                 Teuchos::RCP<Epetra_Vector> systemvector,
                                 Teuchos::RCP<Epetra_Vector> systemvectord,
                                 Teuchos::RCP<Epetra_Vector> systemvectordd,
                                 Teuchos::RCP<Epetra_Vector> toggle,
                                 Teuchos::RCP<std::set<int> > dbcgids);


//--------------------------------------------------------------------
/*!
\brief evaluate spatial function (public)
\author g.bau
\date 03/07
*/
//--------------------------------------------------------------------
static double EvaluateFunction(DRT::Node*        node,
		               int               index,
			       int		 funct_num);


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params,
                                            Teuchos::RCP<Epetra_Vector> systemvector,
                                            Teuchos::RCP<Epetra_Vector> systemvectord,
                                            Teuchos::RCP<Epetra_Vector> systemvectordd,
                                            Teuchos::RCP<Epetra_Vector> toggle,
                                            Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::null;
  if (dbcmapextractor != Teuchos::null) dbcgids = Teuchos::rcp(new std::set<int>());

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in systemvector.
  // For this reason, Dirichlet BCs are evaluated hierarchical meaning
  // in this order:
  //                VolumeDirichlet
  //                SurfaceDirichlet
  //                LineDirichlet
  //                PointDirichlet
  // This way, lower entities override higher ones which is
  // equivalent to inheritance of dirichlet BCs as done in the old
  // ccarat discretization with design          (mgee 1/07)

  // Do VolumeDirichlet first
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }

  // create DBC and free map and build their common extractor
  if (dbcmapextractor != Teuchos::null)
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids->size() > 0)
    {
      dbcgidsv.reserve(dbcgids->size());
      dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap
      = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DofRowMap()->IndexBase(), DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *dbcmapextractor = LINALG::MapExtractor(*(DofRowMap()), dbcmap);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DoDirichletCondition(DRT::Condition&             cond,
                          DRT::Discretization&        dis,
                          const bool                  usetime,
                          const double                time,
                          Teuchos::RCP<Epetra_Vector> systemvector,
                          Teuchos::RCP<Epetra_Vector> systemvectord,
                          Teuchos::RCP<Epetra_Vector> systemvectordd,
                          Teuchos::RCP<Epetra_Vector> toggle,
                          Teuchos::RCP<std::set<int> > dbcgids)
{
  const vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const vector<int>*    curve  = cond.Get<vector<int> >("curve");
  const vector<int>*    funct  = cond.Get<vector<int> >("funct");
  const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
  const vector<double>* val    = cond.Get<vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!dis.NodeRowMap()->MyGID((*nodeids)[i])) continue;
    DRT::Node* actnode = dis.gNode((*nodeids)[i]);
    if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
    vector<int> dofs = dis.Dof(actnode);
    const unsigned numdf = dofs.size();
    for (unsigned j=0; j<numdf; ++j)
    {
      if ((*onoff)[j]==0)
      {
        const int lid = (*systemvectoraux).Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (toggle!=Teuchos::null)
          (*toggle)[lid] = 0.0;
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
        continue;
      }
      const int gid = dofs[j];
      vector<double> value(deg+1, (*val)[j]);

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[j];
      if (curvenum>=0 && usetime)
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;
      //cout << "Dirichlet curve factor: ";
      //for (unsigned i=0; i<deg; ++i) cout << curvefac[i] << ", ";
      //cout << curvefac[deg] << endl;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct) funct_num = (*funct)[j];
      {
         if (funct_num>0)
           functfac = EvaluateFunction(actnode,j,funct_num);
      }
      //cout << "Dirichlet value " << value << " functfac " <<  functfac << endl;

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
      }

      // assign value
      const int lid = (*systemvectoraux).Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      if (systemvector != Teuchos::null)
        (*systemvector)[lid] = value[0];
      if (systemvectord != Teuchos::null)
        (*systemvectord)[lid] = value[1];
      if (systemvectordd != Teuchos::null)
        (*systemvectordd)[lid] = value[2];
      // set toggle vector
      if (toggle != Teuchos::null)
        (*toggle)[lid] = 1.0;
      // amend vector of DOF-IDs which are Dirichlet BCs
      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  RefCountPtr<Epetra_Vector> systemvector,
  const string& condstring,
  const int condid
)
{
  EvaluateCondition(params,Teuchos::null,Teuchos::null,systemvector,Teuchos::null,Teuchos::null,condstring,condid);
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  const string& condstring,
  const int condid
)
{
  EvaluateCondition(params,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condstring,condid);
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 02/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  RCP<LINALG::SparseOperator> systemmatrix1,
  RCP<LINALG::SparseOperator> systemmatrix2,
  RCP<Epetra_Vector> systemvector1,
  RCP<Epetra_Vector> systemvector2,
  RCP<Epetra_Vector> systemvector3,
  const string& condstring,
  const int condid
)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.Getint("ConditionID"))
      {
        map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
        // if (geom.empty()) dserror("evaluation of condition with empty geometry");
        // no check for empty geometry here since in parallel computations
        // can exist processors which do not own a portion of the elements belonging
        // to the condition geometry
        map<int,RefCountPtr<DRT::Element> >::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameterlist
        const vector<int>*    curve  = cond.Get<vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
          curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

        // Get ConditionID of current condition if defined and write value in parameterlist
        const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID",(*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
          params.set(factorname,curvefac);
        }
        else
        {
          params.set("LoadCurveFactor",curvefac);
        }
        params.set<RefCountPtr<DRT::Condition> >("condition", fool->second);

        // define element matrices and vectors
        Epetra_SerialDenseMatrix elematrix1;
        Epetra_SerialDenseMatrix elematrix2;
        Epetra_SerialDenseVector elevector1;
        Epetra_SerialDenseVector elevector2;
        Epetra_SerialDenseVector elevector3;

        for (curr=geom.begin(); curr!=geom.end(); ++curr)
        {
          // get element location vector and ownerships
          vector<int> lm;
          vector<int> lmowner;
          curr->second->LocationVector(*this,lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        if (assemblemat1) elematrix1.Shape(eledim,eledim);
        if (assemblemat2) elematrix2.Shape(eledim,eledim);
        if (assemblevec1) elevector1.Size(eledim);
        if (assemblevec2) elevector2.Size(eledim);
        if (assemblevec3) elevector3.Size(eledim);

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
                                           elevector1,elevector2,elevector3);
          if (err) dserror("error while evaluating elements");

          // assembly
          int eid = curr->second->Id();
          if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,lm,lmowner);
          if (assemblemat2) systemmatrix2->Assemble(eid,elematrix2,lm,lmowner);
          if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
          if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
          if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);
        }
      }
    }
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // end of DRT::Discretization::EvaluateCondition


/*----------------------------------------------------------------------*
 |  evaluate a condition on a surface using parent data        (public) |
 |                                                          gammi 07/08 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateConditionUsingParentData(
  ParameterList&                       params       ,
  RCP<LINALG::SparseOperator> systemmatrix1,
  RCP<LINALG::SparseOperator> systemmatrix2,
  RCP<Epetra_Vector>          systemvector1,
  RCP<Epetra_Vector>          systemvector2,
  RCP<Epetra_Vector>          systemvector3,
  const string&                        condstring   ,
  const int                            condid       )
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  multimap<string,RefCountPtr<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.Getint("ConditionID"))
      {
	map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
	// no check for empty geometry here since in parallel computations
	// can exist processors which do not own a portion of the elements belonging
	// to the condition geometry

	map<int,RefCountPtr<DRT::Element> >::iterator curr;

	// stuff the whole condition into the parameterlist
	// --- we want to be able to access boundary values
	// on the element level
	params.set<RefCountPtr<DRT::Condition> >("condition", fool->second);

	// define element matrices and vectors
	Epetra_SerialDenseMatrix elematrix1;
	Epetra_SerialDenseMatrix elematrix2;
	Epetra_SerialDenseVector elevector1;
	Epetra_SerialDenseVector elevector2;
	Epetra_SerialDenseVector elevector3;

	// element matrices and vectors will be reshaped
	// during the element call!

	for (curr=geom.begin(); curr!=geom.end(); ++curr)
	{
	  // get element location vector and ownerships
	  vector<int> lm;
	  vector<int> lmowner;
	  curr->second->LocationVector(*this,lm,lmowner);

	  // place vectors for parent lm and lmowner in
	  // the parameterlist --- the element will fill
	  // them since only the element implementation
	  // knows its parent
	  RefCountPtr<vector<int> > plm     =rcp(new vector<int>);
	  RefCountPtr<vector<int> > plmowner=rcp(new vector<int>);

	  params.set<RefCountPtr<vector<int> > >("plm",plm);
	  params.set<RefCountPtr<vector<int> > >("plmowner",plmowner);

	  // call the element specific evaluate method
	  int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
					   elevector1,elevector2,elevector3);
	  if (err) dserror("error while evaluating elements");

	  // assembly to all parent dofs even if we just integrated
	  // over a boundary element
	  int eid = curr->second->Id();

	  if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,*plm,*plmowner);
	  if (assemblemat2) systemmatrix2->Assemble(eid,elematrix2,*plm,*plmowner);
	  if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,*plm,*plmowner);
	  if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,*plm,*plmowner);
	  if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,*plm,*plmowner);
	} // end loop geometry elements of this conditions
      } // the condition number is as desired or we wanted to evaluate
        // all numbers
    } // end we have a condition of type condstring
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // end of DRT::Discretization::EvaluateConditionUsingParentData

/*----------------------------------------------------------------------*
 |  evaluate spatial function (public)                       g.bau 03/07|
 *----------------------------------------------------------------------*/
double EvaluateFunction(DRT::Node*      node,
                        int index,
                        int funct_num)
{
  return DRT::UTILS::FunctionManager::Instance().Funct(funct_num-1).Evaluate(index,node->X());
}


/*----------------------------------------------------------------------*
 |  evaluate/assemble scalars across elements (public)       bborn 08/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(
  Teuchos::ParameterList& params,
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // number of scalars
  const int numscalars = scalars->Length();
  if (numscalars <= 0) dserror("scalars vector of interest has size <=0");
  // intermediate sum of each scalar on each processor
  Epetra_SerialDenseVector cpuscalars(numscalars);

  // define element matrices and vectors
  // -- which are empty and unused, just to satisfy element Evaluate()
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over _row_ elements
  const int numrowele = NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    // pointer to current element
    DRT::Element* actele = lRowElement(i);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;  // location vector
    std::vector<int> lmowner;  // processor which owns DOFs
    actele->LocationVector(*this,lm,lmowner);

    // define element vector
    Epetra_SerialDenseVector elescalars(numscalars);

    // call the element evaluate method
    {
      int err = actele->Evaluate(params,*this,lm,
                                 elematrix1,elematrix2,elescalars,elevector2,elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    // sum up (on each processor)
    cpuscalars += elescalars;
  } // for (int i=0; i<numrowele; ++i)

  // reduce
  for (int i=0; i<numscalars; ++i) (*scalars)(i) = 0.0;
  Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

  // bye
  return;
}  // DRT::Discretization::EvaluateScalars


#endif  // #ifdef CCADISCRET
