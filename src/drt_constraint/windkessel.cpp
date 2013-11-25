/*!----------------------------------------------------------------------
\file windkessel.cpp

\brief Basic Windkessel class, dealing with Windkessel boundary conditions

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/



#include <iostream>

#include "windkessel.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Windkessel::Windkessel(RCP<DRT::Discretization> discr,
        const std::string& conditionname,
        int& offsetID,
        int& maxID,
        std::vector<int>& curID):
actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,windkesselcond_);
  if (windkesselcond_.size())
  {
  	windkesseltype_=GetWindkesselType(conditionname);
		std::vector<int> curcoupID;
    for (unsigned int i=0; i<windkesselcond_.size();i++)
    {
      //windkesselcond_[i]->Print(std::cout);
      int condID=(*(windkesselcond_[i]->Get<std::vector<int> >("id")))[0];

      //std::vector<int> curID(i);
      curID.push_back(windkesselcond_[i]->GetInt("id"));

      if (condID>maxID)
      {
        maxID=condID;
      }
      if (condID<offsetID)
      {
        offsetID=condID;
      }

      std::vector<DRT::Condition*> surfneumcond;
      std::vector<int> tmp;
      Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
      if (structdis == Teuchos::null) dserror("no structure discretization available");

      // first get all Neumann conditions on structure
      structdis->GetCondition("SurfaceNeumann",surfneumcond);
      unsigned int numneumcond = surfneumcond.size();
      if (numneumcond == 0) dserror("no Neumann conditions on structure");

      // now filter those Neumann conditions that are due to the coupling
      std::vector<DRT::Condition*> coupcond;
      for (unsigned int k = 0; k < numneumcond; ++k)
      {
        DRT::Condition* actcond = surfneumcond[k];
        if (actcond->Type() == DRT::Condition::WindkesselStructureCoupling)
    	  coupcond.push_back(actcond);
      }
      unsigned int numcond = coupcond.size();

			curcoupID.push_back((coupcond[i])->GetInt("coupling id"));
			if (curID.size() != curcoupID.size()) dserror("Coupling conditions do not match Windkessel conditions!");

      if (numcond == 0) dserror("no coupling conditions found");

      int coupcondID = 0;

			DRT::Condition* cond = coupcond[i];
			cond->Add("type","neum_orthopressure");
			std::vector<int> onoff(6,0);
			onoff[0] = 1;
			cond->Add("onoff",onoff);
			std::vector<double> val(6,0.0);
			cond->Add("val",val);
			coupcondID = (coupcond[i])->GetInt("coupling id");
			tmp.push_back(coupcondID);

      if (curID[i] != curcoupID[i]) dserror("Choose the same ids for the Windkessel and the structural coupling surface!");

    }
  }
  else
  {
	windkesseltype_=none;
  }
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Windkessel::Windkessel(RCP<DRT::Discretization> discr,
        const std::string& conditionname):
actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,windkesselcond_);
  
  if (windkesselcond_.size())
  {
	windkesseltype_=GetWindkesselType(conditionname);
    for (unsigned int i=0; i<windkesselcond_.size();i++)
    {

      int condID=windkesselcond_[i]->GetInt("id");
      inittimes_.insert(std::pair<int,double>(condID,0.0));
      activecons_.insert(std::pair<int,bool>(condID,false));

    }
  }
  else
  {
	windkesseltype_=none;
  }
}

/*-----------------------------------------------------------------------*
|(private)                                                      mhv 10/13|
*-----------------------------------------------------------------------*/
UTILS::Windkessel::WindkesselType UTILS::Windkessel::GetWindkesselType(const std::string& name)
{
  if (name=="WindkesselStructureCond")
    return rc;
  return none;
}

/*------------------------------------------------------------------------*
|(public)                                                      mhv 10/13  |
|Initialization routine computes ref base values and activates conditions |
*------------------------------------------------------------------------*/
void UTILS::Windkessel::Initialize(
    Teuchos::ParameterList&        params,
    RCP<Epetra_Vector>    systemvector3)
{

  params.set("action","calc_struct_constrvol");

  // start computing
  InitializeWindkessel(params,systemvector3);
  return;
}


/*-----------------------------------------------------------------------*
|(public)                                                       mhv 10/13|
|Evaluate Windkessel functions, choose the right action based on type             |
*-----------------------------------------------------------------------*/
void UTILS::Windkessel::Evaluate(
    Teuchos::ParameterList&        params,
    RCP<LINALG::SparseMatrix> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
    RCP<Epetra_Vector>    systemvector2,
    RCP<Epetra_Vector>    systemvector3)
{

  params.set("action","calc_struct_volconstrstiff");

  EvaluateWindkessel(params,systemmatrix1,systemmatrix2,systemvector1,systemvector2,systemvector3);
  return;
}

/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 10/13 |
 |Evaluate method, calling element evaluates of a condition and          |
 |assembing results based on this conditions                             |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateWindkessel(
    Teuchos::ParameterList&        params,
    RCP<LINALG::SparseMatrix> systemmatrix1,
    RCP<LINALG::SparseOperator> systemmatrix2,
    RCP<Epetra_Vector>    systemvector1,
		RCP<Epetra_Vector>    systemvector2,
		RCP<Epetra_Vector>    systemvector3)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_timint = params.get("scale_timint",1.0);
  double gamma = params.get("scale_gamma",1.0);
  double ts_size = params.get("time_step_size",1.0);

	std::vector<double> resistance;
	std::vector<double> invresist;
	std::vector<double> compliance;
	std::vector<double> wkstiff;
	std::vector<bool> havegid;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

    if (assemblemat1 or assemblevec1 or assemblevec2)
    {
			resistance.push_back(windkesselcond_[i]->GetDouble("resistance"));
			compliance.push_back(windkesselcond_[i]->GetDouble("compliance"));
			invresist.push_back(1./resistance[i]);
    }

		// is condition already labeled as active?
		if(activecons_.find(condID)->second==false)
		{
			const std::string action = params.get<std::string>("action");
			RCP<Epetra_Vector> displast=params.get<RCP<Epetra_Vector> >("old disp");
			actdisc_->SetState("displacement",displast);
			RCP<Epetra_Vector> disp=params.get<RCP<Epetra_Vector> >("new disp");
			actdisc_->SetState("displacement",disp);
			params.set("action",action);
		}

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = condID-offsetID;
    //std::cout << "" << offsetID << std::endl;

		// elements might need condition
		params.set<RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

		// assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
		if (assemblemat1)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			systemmatrix1->UnComplete();
			wkstiff.push_back(sc_timint*(compliance[i]/(gamma*ts_size)+(1./resistance[i])));

			havegid.push_back(systemmatrix1->RowMap().MyGID(gindex));
			if(havegid[i]) systemmatrix1->Assemble(wkstiff[i],colvec[0],colvec[0]);

		}
		// assemble the constant parts of the Windkessel rhs vector
		if (assemblevec1)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err2 = systemvector1->SumIntoGlobalValues(1,&invresist[i],&colvec[0]);
			if (err2) dserror("SumIntoGlobalValues failed!");
		}
		if (assemblevec2)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err3 = systemvector2->SumIntoGlobalValues(1,&compliance[i],&colvec[0]);
			if (err3) dserror("SumIntoGlobalValues failed!");
		}

		// define element matrices and vectors
		Epetra_SerialDenseMatrix elematrix1;
		Epetra_SerialDenseMatrix elematrix2;
		Epetra_SerialDenseVector elevector1;
		Epetra_SerialDenseVector elevector2;
		Epetra_SerialDenseVector elevector3;

		std::map<int,RCP<DRT::Element> >& geom = cond.Geometry();
		// if (geom.empty()) dserror("evaluation of condition with empty geometry");
		// no check for empty geometry here since in parallel computations
		// can exist processors which do not own a portion of the elements belonging
		// to the condition geometry
		std::map<int,RCP<DRT::Element> >::iterator curr;
		for (curr=geom.begin(); curr!=geom.end(); ++curr)
		{
			// get element location vector and ownerships
			std::vector<int> lm;
			std::vector<int> lmowner;
			std::vector<int> lmstride;
			curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

			// get dimension of element matrices and vectors
			// Reshape element matrices and vectors and init to zero
			const int eledim = (int)lm.size();

			elematrix2.Shape(eledim,eledim);
			elevector2.Size(eledim);
			elevector3.Size(1);

			// call the element specific evaluate method
			int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
					elevector1,elevector2,elevector3);
			if (err) dserror("error while evaluating elements");

			// assembly
			int eid = curr->second->Id();

			if (assemblemat2)
			{
				// assemble to rectangular matrix. The col corresponds to the Windkessel ID.
				std::vector<int> colvec(1);
				colvec[0]=gindex;
				systemmatrix2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
			}
			if (assemblevec3)
			{
				std::vector<int> windkessellm;
				std::vector<int> windkesselowner;
				windkessellm.push_back(gindex);
				windkesselowner.push_back(curr->second->Owner());
				LINALG::Assemble(*systemvector3,elevector3,windkessellm,windkesselowner);
			}

    }

  }
  return;
} // end of EvaluateCondition

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeWindkessel(
    Teuchos::ParameterList&        params,
    RCP<Epetra_Vector>    systemvector)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist

    int condID=cond.GetInt("id");
    params.set("id",condID);

		params.set<RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

		// define element matrices and vectors
		Epetra_SerialDenseMatrix elematrix1;
		Epetra_SerialDenseMatrix elematrix2;
		Epetra_SerialDenseVector elevector1;
		Epetra_SerialDenseVector elevector2;
		Epetra_SerialDenseVector elevector3;

		std::map<int,RCP<DRT::Element> >& geom = cond.Geometry();
		// no check for empty geometry here since in parallel computations
		// can exist processors which do not own a portion of the elements belonging
		// to the condition geometry
		std::map<int,RCP<DRT::Element> >::iterator curr;
		for (curr=geom.begin(); curr!=geom.end(); ++curr)
		{
			// get element location vector and ownerships
			std::vector<int> lm;
			std::vector<int> lmowner;
			std::vector<int> lmstride;
			curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

			// get dimension of element matrices and vectors
			// Reshape element matrices and vectors and init to zero
			elevector3.Size(1);

			// call the element specific evaluate method
			int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
																			 elevector1,elevector2,elevector3);
			if (err) dserror("error while evaluating elements");

			// assembly

			std::vector<int> windkessellm;
			std::vector<int> windkesselowner;
			int offsetID = params.get<int> ("OffsetID");
			windkessellm.push_back(condID-offsetID);
			windkesselowner.push_back(curr->second->Owner());
			LINALG::Assemble(*systemvector,elevector3,windkessellm,windkesselowner);
		}
		// remember next time, that this condition is already initialized, i.e. active
		activecons_.find(condID)->second=true;

		if (actdisc_->Comm().MyPID()==0)
		{
			std::cout << "===== Welcome to monolithic 3D structure 0D Windkessel coupling (coupling id = " << condID << ") ====="<< std::endl;
		}


  }
  return;
} // end of Initialize Windkessel

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::vector<int> UTILS::Windkessel::GetActiveCondID()
{
  std::vector<int> condID;
  std::map<int,bool>::const_iterator mapit;
  for(mapit = activecons_.begin();mapit!=activecons_.end();mapit++)
  {
    if (mapit->second)
      condID.push_back(mapit->first);
  }
  return condID;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::SetState
(
  const std::string& state,  ///< name of state to set
  RCP<Epetra_Vector> V  ///< values to set
)
{
  actdisc_->SetState(state,V);
}

