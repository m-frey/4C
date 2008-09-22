/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_base_algorithm.cpp

\brief scalar transport field base algorithm

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_scatra_base_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  /// setup scalar transport algorithm (overriding some dynamic parameters with
  /// values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numscatra,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------

  const Teuchos::ParameterList& scatradyn = 
    DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // print out default parameters of scalar tranport parameter list
  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, scatradyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  SOLVAR *actsolv  = &solv[genprob.numscatra];

  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> scatratimeparams= rcp(new ParameterList());

  // ----problem type (type of scalar transport problem we want to solve)
  scatratimeparams->set<string>("problem type",DRT::Problem::Instance()->ProblemType());

  // --------------------type of time-integration (or stationary) scheme
  INPUTPARAMS::ScaTraTimeIntegrationScheme timintscheme =
    Teuchos::getIntegralValue<INPUTPARAMS::ScaTraTimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  scatratimeparams->set<INPUTPARAMS::ScaTraTimeIntegrationScheme>("time int algo",timintscheme);

  // --------------------------------------- time integration parameters
  // the default time step size
  scatratimeparams->set<double>   ("time step size"           ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  scatratimeparams->set<double>   ("total time"               ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  scatratimeparams->set<int>      ("max number timesteps"     ,prbdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  scatratimeparams->set           ("write restart every"       ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  scatratimeparams->set           ("write solution every"      ,prbdyn.get<int>("UPRES"));
  //scatratimeparams->set           ("write solution every"      ,prbdyn.get<int>("WRITESOLEVRY"));
  // write also flux vectors when solution is written out?
  scatratimeparams->set<string>   ("write flux"   ,scatradyn.get<string>("WRITEFLUX"));

  // ---------------------------------------------------- initial field
  scatratimeparams->set<int>("scalar initial field"     ,Teuchos::getIntegralValue<int>(scatradyn,"INITIALFIELD"));
  scatratimeparams->set<int>("scalar initial field func number",scatradyn.get<int>("INITFUNCNO"));
  
  // ----------------------------------------------------velocity field
  scatratimeparams->set<int>("velocity field"     ,Teuchos::getIntegralValue<int>(scatradyn,"VELOCITYFIELD"));
  scatratimeparams->set<int>("velocity function number",scatradyn.get<int>("VELFUNCNO"));
  
  // -------------------------------- (fine-scale) subgrid diffusivity?
  scatratimeparams->set<string>("fs subgrid diffusivity"   ,scatradyn.get<string>("FSSUGRVISC"));


  // --------------sublist for combustion-specific gfunction parameters
  /* This sublist COMBUSTION DYNAMIC/GFUNCTION contains parameters for the gfunction field
   * which are only relevant for a combustion problem.                         07/08 henke */
  if (genprob.probtyp == prb_combust)
  {
    scatratimeparams->sublist("COMBUSTION GFUNCTION")=prbdyn.sublist("COMBUSTION GFUNCTION");
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm construction depending on 
  // respective time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timintscheme == INPUTPARAMS::timeint_stationary)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for stationary scheme
    // -----------------------------------------------------------------

    // parameter theta for time-integration schemes
    scatratimeparams->set<double>           ("theta"                    ,1.0);

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntStationary::TimIntStationary(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_one_step_theta)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for one-step-theta scheme
    // -----------------------------------------------------------------

    // parameter theta for time-integration schemes
    scatratimeparams->set<double>           ("theta"                    ,scatradyn.get<double>("THETA"));

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_bdf2)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for BDF2 scheme
    // -----------------------------------------------------------------

    // parameter theta for time-integration schemes
    scatratimeparams->set<double>           ("theta"                    ,scatradyn.get<double>("THETA"));

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntBDF2::TimIntBDF2(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    scatratimeparams->set<double>           ("alpha_M"                  ,scatradyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    scatratimeparams->set<double>           ("alpha_F"                  ,scatradyn.get<double>("ALPHA_F"));

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    dserror("no adapter for generalized alpha scalar transport dynamic routine implemented.");
  }
  else
  {
    dserror("Unknown time integration scheme for scalar tranport problem");
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::~ScaTraBaseAlgorithm()
{
}


#endif  // #ifdef CCADISCRET
