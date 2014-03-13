/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_graddesc.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "invana_utils.H"
#include "stat_inv_ana_graddesc.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_adapter/ad_str_structure.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* constructor */
STR::INVANA::StatInvAnaGradDesc::StatInvAnaGradDesc(Teuchos::RCP<DRT::Discretization> dis):
  StatInvAnalysis(dis),
stepsize_(0.0),
maxiter_(0),
runc_(0),
convcritc_(0)
{
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // max number of iterations
  maxiter_ = invap.get<int>("MAXITER");

  //set stepsize for gradient scheme
  stepsize_ = invap.get<double>("STEPSIZE");

  //get tolerance
  convtol_ = invap.get<double>("CONVTOL");

  p_= Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumParams(),true));
  step_= Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumParams(), true));

}

/*----------------------------------------------------------------------*/
/* do the optimization loop*/
void STR::INVANA::StatInvAnaGradDesc::Optimize()
{
  int success=0;

  // solve initially to get quantities:
  SolveForwardProblem();
  SolveAdjointProblem();
  EvaluateGradient();
  EvaluateError();

  //check gradient by fd:
#if 0
    std::cout << "gradient: " << *objgrad_ << std::endl;
    EvaluateGradientFD();
    std::cout << "gradient approx: " << *objgrad_ << std::endl;
#endif

  objgrad_o_->Update(1.0, *objgrad_, 0.0);

  //get search direction from gradient:
  p_->Update(-1.0, *objgrad_o_, 0.0);

  objval_o_ = objval_;

  MVNorm(objgrad_o_,2,&convcritc_,discret_->ElementRowMap());

  PrintOptStep(0,0);

  while (convcritc_ > convtol_ && runc_<maxiter_)
  {
    double tauopt=1.0;
    int numsteps=0;

    // do the line search
    success = EvaluateArmijoRule(&tauopt, &numsteps);

    if (success == 1)
    {
      std::cout << " Line Search Break Down" << std::endl;
      break;
    }

    //get the L2-norm:
    MVNorm(objgrad_,2,&convcritc_,discret_->ElementRowMap());

    //compute new direction only for runs
    p_->Update(-1.0, *objgrad_, 0.0);

    // bring quantities to the next run
    objgrad_o_->Update(1.0, *objgrad_, 0.0);
    objval_o_=objval_;
    runc_++;

    //do some on screen printing
    PrintOptStep(tauopt, numsteps);
  }

  Summarize();

  return;
}

/*----------------------------------------------------------------------*/
/* do a line search based on armijo rule */
int STR::INVANA::StatInvAnaGradDesc::EvaluateArmijoRule(double* tauopt, int* numsteps)
{
  int i=0;
  int imax=20;
  double c1=1.0e-4;
  double tau_max=1.0e10;
  double gnorm=0.0;

  // "last"/"intermediate" values for cubic model
  // these are actually safeguardly set after the first call to the quadratic model
  double tau_l=0.0;
  double e_l=0.0;

  int success=0;

  //safeguard multiplicators
  double blow=0.1;
  double bhigh=0.5;

  MVNorm(objgrad_o_,2,&gnorm,discret_->ElementRowMap());

  double tau_n=std::min(1.0, 100/(1+gnorm));
  //std::cout << "trial step size: " << tau_n << std::endl;

  while (i<imax && tau_n<tau_max)
  {
    // step based on current stepsize
    step_->Update(tau_n, *p_, 0.0);

    //make a step
    matman_->UpdateParams(step_);
    SolveForwardProblem();
    SolveAdjointProblem();
    EvaluateGradient();
    EvaluateError();

    // check sufficient decrease:
    double dfp_o=0.0;
    MVDotProduct(objgrad_o_,p_,&dfp_o,discret_->ElementRowMap());

    if ( (objval_-objval_o_) < c1*tau_n*dfp_o )
    {
      *tauopt=tau_n;
      *numsteps=i+1;
      return 0;
    }

    // do stepsize prediction based on polynomial models
    if (i==0)
      success=polymod(objval_o_, dfp_o,tau_n,objval_,blow,bhigh,tauopt);
    else
      success=polymod(objval_o_,dfp_o,tau_n,objval_,blow,bhigh,tau_l,e_l,tauopt);

    if (success==1) return 1;

    e_l=objval_;
    tau_l=tau_n;
    tau_n=*tauopt;
    matman_->ResetParams();
    i++;

#if 0
    //brute force sampling:
    if (runc_==-1)
    {
      for (int k=0; k<100; k++)
      {
        step_->Update(tau_n*k/100, *p_, 0.0);
        matman_->UpdateParams(step_);
        SolveForwardProblem();
        double ee = objfunct_->Evaluate(dis_,matman_);
        std::cout << " run " << tau_n*k/100 << " " << ee << endl;
        matman_->ResetParams();
      }
    }
#endif
  }

  return 1;
}

/*----------------------------------------------------------------------*/
/* quadratic model */
int STR::INVANA::StatInvAnaGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  *tauopt=-dfp/(2*tau_n*(e_n-e_o-dfp));
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}

/*----------------------------------------------------------------------*/
/* cubic model model */
int STR::INVANA::StatInvAnaGradDesc::polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt)
{
  double lleft=tau_n*blow;
  double lright=tau_n*bhigh;

  double a1=tau_n*tau_n;
  double a2=tau_n*tau_n*tau_n;
  double a3=tau_l*tau_l;
  double a4=tau_l*tau_l*tau_l;

  double deta=a1*a4-a2*a3;

  if (deta<1.0e-14)
    return 1;

  double b1=e_n-(e_o+dfp*tau_n);
  double b2=e_l-(e_o+dfp*tau_l);

  double c1=1/deta*(a4*b1-a2*b2);
  double c2=1/deta*(-a3*b1+a1*b2);

  *tauopt=(-c1+sqrt(c1*c1-3*c2*dfp))/(3*c2);
  if (*tauopt < lleft) *tauopt = lleft;
  if (*tauopt > lright) *tauopt = lright;

  return 0;
}


/*----------------------------------------------------------------------*/
/* print final results*/
void STR::INVANA::StatInvAnaGradDesc::PrintOptStep(double tauopt, int numsteps)
{
  printf("OPTIMIZATION STEP %3d | ", runc_);
  printf("Objective function: %10.8e | ", objval_o_);
  printf("Gradient : %10.8e | ", convcritc_);
  printf("stepsize : %10.8e | LSsteps %2d\n", tauopt, numsteps);
  fflush(stdout);

}

/*----------------------------------------------------------------------*/
/* print final results*/
void STR::INVANA::StatInvAnaGradDesc::Summarize()
{
  std::cout << "the final vector of parameters: " << std::endl;
  std::cout << *(matman_->GetParams()) << std::endl;
  return;
}

