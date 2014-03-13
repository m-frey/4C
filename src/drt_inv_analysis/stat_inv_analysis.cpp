/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_analysis.H"
#include "invana_utils.H"
#include "invana_resulttest.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_comm/comm_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/drt_validparameters.H"

#include "timint_adjoint.H"
#include "matpar_manager.H"
#include "matpar_manager_uniform.H"
#include "objective_funct_disp.H"
#include "objective_funct_surfcurr.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 10/13  */
/*----------------------------------------------------------------------*/
STR::INVANA::StatInvAnalysis::StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis):
discret_(dis),
dofrowmap_(NULL),
regweight_(0.0)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();


  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  // this is supposed to be the number of simulation steps in the primal AND the dual problem
  msteps_ = sdyn.get<int>("NUMSTEP");
  double timestep = sdyn.get<double>("TIMESTEP");

  // initialize "state" vectors
  dis_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));

  // initialize the vector of time steps according to the structural dynamic params
  time_ = Teuchos::rcp(new std::vector<double>(msteps_,0.0));
  for (int i=0; i<=msteps_-1; i++)
  {
    (*time_)[i] = (i+1)*timestep;
  }

  // set up an objective function
  switch (DRT::INPUT::IntegralValue<INPAR::STR::StatInvObjFunctType>(statinvp,"OBJECTIVEFUNCT"))
  {
    case INPAR::STR::stat_inv_obj_disp:
    {
      objfunct_ = Teuchos::rcp(new STR::INVANA::ObjectiveFunctDisp(discret_, msteps_, time_));
    }
    break;
    case INPAR::STR::stat_inv_obj_surfcurr:
    {
      objfunct_ = Teuchos::rcp(new STR::INVANA::ObjectiveFunctSurfCurrRepresentation(discret_, msteps_, time_));
    }
    break;
    case INPAR::STR::stat_inv_obj_none:
    {
      dserror("choose some type of objective function");
    }
    break;
  }

  // do we have regularization!
  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvRegularization>(statinvp,"REGULARIZATION"))
  {
    case INPAR::STR::stat_inv_reg_none:
    {
      havereg_ = false;
    }
    break;
    case INPAR::STR::stat_inv_reg_thikonov:
    {
      // The most simple one can do. This will eventually be moved somewhere but for now we do it in here since
      // it is actually "one" line of code which will be evaluated when calling the respective
      // routines on the objective function
      havereg_ = true;
      regweight_ = statinvp.get<double>("REG_WEIGHT");
    }
      break;
  }

  // set up the material parameter handler:
  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvMatParametrization>(statinvp,"PARAMETRIZATION"))
  {
    case INPAR::STR::stat_inv_mp_smoothkernel:
    {
      dserror("no parametrization based on gaussian kernels yet!");
    }
    break;
    case INPAR::STR::stat_inv_mp_elementwise:
    {
      matman_ = Teuchos::rcp(new STR::INVANA::MatParManagerPerElement(discret_));
    }
    break;
    case INPAR::STR::stat_inv_mp_uniform:
    {
      matman_ = Teuchos::rcp(new STR::INVANA::MatParManagerUniform(discret_));
    }
      break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
    break;
  }

  //set these infeasibly high:
  objval_ = 1.0e17;
  objval_o_ = 1.0e16;
  error_incr_ = 1.0e16;

  objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumParams(),true));
  objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumParams(),true));


}

/*----------------------------------------------------------------------*/
/* MStep EpetraVector to EpetraMultiVector                   keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::MStepEpetraToEpetraMulti(Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > mstepvec,
                                                            Teuchos::RCP<Epetra_MultiVector> multivec)
{
  for (int i=0; i<msteps_; i++)
    (*multivec)(i)->Update(1.0,*(*mstepvec)(-(msteps_-i-1)),0.0);
}

/*----------------------------------------------------------------------*/
/* Mstep double to std::vector<double>                       keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::MStepDToStdVecD(Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > mstepvec,
                                                   Teuchos::RCP<std::vector<double> > stdvec)
{
  for (int i=0; i<msteps_; i++)
    (*stdvec)[i] = *(*mstepvec)(-(msteps_-i-1));

}

/*----------------------------------------------------------------------*/
/* solve primal problem                                      keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::SolveForwardProblem()
{
  // resultfiles are overwritten every run since usually only the final results are of interest
  discret_->Writer()->OverwriteResultFile();

  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
    case INPAR::STR::dyna_statics:
    {
      ADAPTER::StructureBaseAlgorithm adapterbase(sdyn,const_cast<Teuchos::ParameterList&>(sdyn), discret_);
      ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

      // do restart
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        structadaptor.ReadRestart(restart);
      }
      structadaptor.Integrate();

      // get displacement and time
      MStepEpetraToEpetraMulti(structadaptor.DispMStep(), dis_);
      MStepDToStdVecD(structadaptor.TimeMStep(), time_);

      break;
    }
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_expleuler:
    case INPAR::STR::dyna_centrdiff:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma:
    case INPAR::STR::dyna_euimsto:
      dserror("return of multistep-variables only for static analysis (so far)");
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  return;
}


/*----------------------------------------------------------------------*/
/* solve dual problem                                       keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::SolveAdjointProblem()
{
  //Setup RHS for the adjoints
  Teuchos::RCP<Epetra_MultiVector> objgrad;
  objgrad = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  objfunct_->EvaluateGradient(dis_,objgrad);

  //initialize adjoint time integration with RHS as input
  STR::TimIntAdjoint timintadj = TimIntAdjoint(discret_, *time_);
  timintadj.SetupAdjoint(objgrad, dis_);

  // adjoint time integration
  timintadj.Integrate();

  // get the solution
  disdual_ = timintadj.ExtractSolution();

}

/*----------------------------------------------------------------------*/
/* evaluate gradient of the objective function                          */
/* using the adjoint equations                              keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::EvaluateGradient()
{
  //zero out gradient vector initially
  objgrad_->Scale(0.0);

  Teuchos::RCP<Epetra_Vector> zeros;
  zeros = LINALG::CreateVector(*(dofrowmap_), true);
  zeros->Scale(0.0);

  //loop the time steps
  for (int j=0; j<msteps_; j++)
  {
    discret_->SetState(0, "displacement", Teuchos::rcp((*dis_)(j),false));
    discret_->SetState(0, "residual displacement", zeros);
    discret_->SetState(0, "dual displacement", Teuchos::rcp((*disdual_)(j),false));

    matman_->Evaluate((*time_)[j], objgrad_);
  }

  if (havereg_)
  {
    //simple tikhonov regularization on the parameter vector
    objgrad_->Update(regweight_,*(matman_->GetParams()),1.0);
  }

}

void STR::INVANA::StatInvAnalysis::ResetDiscretization()
{
  Teuchos::ParameterList p;
  p.set("action","calc_struct_reset_all");
  discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
}

/*----------------------------------------------------------------------*/
/* FD approximation of the gradient                         keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::EvaluateGradientFD()
{
  if (discret_->Comm().NumProc()>1) dserror("this does probably not run in parallel");

  objgrad_->Scale(0.0);
  EvaluateError();
  // we need to keep this!
  double objval0 = objval_;

  int numparams = matman_->NumParams();
  int numele = discret_->ElementColMap()->NumMyElements();

  double perturba = 1.0e-6;
  double perturbb = 1.0e-12;

  Teuchos::RCP<Epetra_MultiVector> perturb = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()), matman_->NumParams(),true));
  perturb->Update(1.0,*(matman_->GetParams()),0.0);

  //keep a copy of the current parameters to reset after perturbation:
  Teuchos::RCP<Epetra_MultiVector> pcurr = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementColMap()), matman_->NumParams(),true));
  pcurr->Update(1.0,*(matman_->GetParams()),0.0);

  // keep a copy of the current displacements correspnding to pcurr:
  Teuchos::RCP<Epetra_MultiVector> discurr = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  discurr->Update(1.0,*dis_,0.0);

  double pn=0.0;
  double p=0.0;
  double dp=0.0;
  for (int i=0; i<numparams; i++)
  {
    for (int j=0; j<numele; j++)
    {
      p = (*(*(matman_->GetParams()))(i))[j];
      pn = p+p*perturba+perturbb;
      perturb->ReplaceGlobalValue(j,i,pn);

      matman_->ReplaceParams(perturb);
      ResetDiscretization();
      SolveForwardProblem();
      perturb->Update(1.0,*pcurr,0.0);

      EvaluateError();
      dp=(objval0-objval_)/(p-pn);
      objgrad_->ReplaceGlobalValue(j,i,dp);
    }
  }

  //reset
  objval_ = objval0;
  matman_->ReplaceParams(pcurr);
  dis_->Update(1.0, *discurr, 0.0);

}

/*----------------------------------------------------------------------*/
/* Evaluate the objective function                          keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::StatInvAnalysis::EvaluateError()
{
  objfunct_->Evaluate(dis_,objval_);

  if (havereg_)
  {
    double val = 0.0;
    STR::INVANA::MVNorm(matman_->GetParams(),2,&val,discret_->ElementRowMap());
    objval_ += 0.5*regweight_*val*val;
  }
}


// return the value of the gradient 2-norm
double STR::INVANA::StatInvAnalysis::GetGrad2Norm()
{
  double res;
  STR::INVANA::MVNorm(objgrad_,2,&res,discret_->ElementRowMap());
  return res;
}


/*----------------------------------------------------------------------*/
/* Creates the field test                                               */
Teuchos::RCP<DRT::ResultTest> STR::INVANA::StatInvAnalysis::CreateFieldTest()
{
  return Teuchos::rcp(new InvAnaResultTest(*this));
}

void STR::INVANA::StatInvAnalysis::PrintDataToScreen(Epetra_MultiVector& vec)
{
  for (int j=0; j<vec.NumVectors(); j++)
  {
    Epetra_Vector tmp(*(vec(j)));
    for (int i=0; i<vec.MyLength(); i++)
    {
      printf("mypid: %2d %.16e \n", vec.Comm().MyPID(), tmp[i]);
    }
  }
}
