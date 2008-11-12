/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "stru_genalpha_zienxie_drt.H"
#include "strugenalpha.H"
#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "stru_resulttest.H"

#include "../drt_inv_analysis/inv_analysis.H"
#include "../drt_statmech/statmech_time.H"
#include "../drt_statmech/bromotion_timeint.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
extern "C"
void caldyn_drt()
{
  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY");

  // major switch to different time integrators
  switch (Teuchos::getIntegralValue<int>(sdyn,"DYNAMICTYP"))
  {
  case STRUCT_DYNAMIC::centr_diff:
    dserror("no central differences in DRT");
    break;
  case STRUCT_DYNAMIC::gen_alfa:
  case STRUCT_DYNAMIC::gen_alfa_statics:
    switch (Teuchos::getIntegralValue<int>(tap,"KIND"))
    {
    case TIMADA_DYNAMIC::timada_kind_none:
      dyn_nlnstructural_drt();
      break;
    case TIMADA_DYNAMIC::timada_kind_zienxie:
      stru_genalpha_zienxie_drt();
      break;
    default:
      dserror("unknown time adaption scheme '%s'", sdyn.get<std::string>("TA_KIND").c_str());
    }
    break;
  case STRUCT_DYNAMIC::Gen_EMM:
    dserror("GEMM not supported");
    break;
  case STRUCT_DYNAMIC::statics:
  case STRUCT_DYNAMIC::genalpha:
  case STRUCT_DYNAMIC::onesteptheta:
  case STRUCT_DYNAMIC::gemm:
  case STRUCT_DYNAMIC::ab2:
    // direct time integration
    STR::strudyn_direct();
    break;
  default:
    dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
  }
}


/*----------------------------------------------------------------------*
  | structural nonlinear dynamics (gen-alpha)              m.gee 12/06  |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(actdis);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[0];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->StructuralContactParams();
  const Teuchos::ParameterList& statmech = DRT::Problem::Instance()->StatisticalMechanicsParams();
  const Teuchos::ParameterList& iap      = DRT::Problem::Instance()->InverseAnalysisParams();
  const Teuchos::ParameterList& bromop   = DRT::Problem::Instance()->BrownianMotionParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),DRT::Problem::Instance()->ErrorFile()->Handle());
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  switch (Teuchos::getIntegralValue<int>(sdyn,"DYNAMICTYP"))
  {
    //==================================================================
    // Generalized alpha time integration
    //==================================================================
    case STRUCT_DYNAMIC::gen_alfa :
    case STRUCT_DYNAMIC::gen_alfa_statics :
    {
      ParameterList genalphaparams;
      StruGenAlpha::SetDefaults(genalphaparams);

      genalphaparams.set<string>("DYNAMICTYP",sdyn.get<string>("DYNAMICTYP"));

      // Rayleigh damping
      genalphaparams.set<bool>  ("damping",(not (sdyn.get<std::string>("DAMPING") == "no"
                                                 or sdyn.get<std::string>("DAMPING") == "No"
                                                 or sdyn.get<std::string>("DAMPING") == "NO")));
      genalphaparams.set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
      genalphaparams.set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

      // Generalised-alpha coefficients
      genalphaparams.set<double>("beta",sdyn.get<double>("BETA"));
#ifdef STRUGENALPHA_BE
      genalphaparams.set<double>("delta",sdyn.get<double>("DELTA"));
#endif
      genalphaparams.set<double>("gamma",sdyn.get<double>("GAMMA"));
      genalphaparams.set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
      genalphaparams.set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

      genalphaparams.set<double>("total time",0.0);
      genalphaparams.set<double>("delta time",sdyn.get<double>("TIMESTEP"));
      genalphaparams.set<double>("max time",sdyn.get<double>("MAXTIME"));
      genalphaparams.set<int>   ("step",0);
      genalphaparams.set<int>   ("nstep",sdyn.get<int>("NUMSTEP"));
      genalphaparams.set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
      genalphaparams.set<int>   ("num iterations",-1);

      genalphaparams.set<string>("convcheck", sdyn.get<string>("CONV_CHECK"));
      genalphaparams.set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
      genalphaparams.set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
      genalphaparams.set<double>("tolerance constraint",sdyn.get<double>("TOLCONSTR"));

      genalphaparams.set<double>("UZAWAPARAM",sdyn.get<double>("UZAWAPARAM"));
      genalphaparams.set<double>("UZAWATOL",sdyn.get<double>("UZAWATOL"));
      genalphaparams.set<int>   ("UZAWAMAXITER",sdyn.get<int>("UZAWAMAXITER"));
      genalphaparams.set<int>   ("UZAWAALGO",getIntegralValue<int>(sdyn,"UZAWAALGO"));
      genalphaparams.set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
      genalphaparams.set<int>   ("io disp every nstep",sdyn.get<int>("RESEVRYDISP"));

      genalphaparams.set<bool>  ("ADAPTCONV",getIntegralValue<int>(sdyn,"ADAPTCONV")==1);
      genalphaparams.set<double>("ADAPTCONV_BETTER",sdyn.get<double>("ADAPTCONV_BETTER"));

      switch (Teuchos::getIntegralValue<STRUCT_STRESS_TYP>(ioflags,"STRUCT_STRESS"))
      {
      case struct_stress_none:
        genalphaparams.set<string>("io structural stress", "none");
      break;
      case struct_stress_cauchy:
        genalphaparams.set<string>("io structural stress", "cauchy");
      break;
      case struct_stress_pk:
        genalphaparams.set<string>("io structural stress", "2PK");
      break;
      default:
        genalphaparams.set<string>("io structural stress", "none");
      break;
      }

      genalphaparams.set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

      switch (Teuchos::getIntegralValue<STRUCT_STRAIN_TYP>(ioflags,"STRUCT_STRAIN"))
      {
      case struct_strain_none:
        genalphaparams.set<string>("io structural strain", "none");
      break;
      case struct_strain_ea:
        genalphaparams.set<string>("io structural strain", "euler_almansi");
      break;
      case struct_strain_gl:
        genalphaparams.set<string>("io structural strain", "green_lagrange");
      break;
      default:
        genalphaparams.set<string>("io structural strain", "none");
      break;
      }

      genalphaparams.set<int>   ("restart",probtype.get<int>("RESTART"));
      genalphaparams.set<int>   ("write restart every",sdyn.get<int>("RESTARTEVRY"));

      genalphaparams.set<bool>  ("print to screen",true);
      genalphaparams.set<bool>  ("print to err",true);
      genalphaparams.set<FILE*> ("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

      // parameters for inverse analysis
      genalphaparams.set<bool>  ("inv_analysis",Teuchos::getIntegralValue<int>(iap,"INV_ANALYSIS"));
      genalphaparams.set<double>("measured_curve0",iap.get<double>("MEASURED_CURVE0"));
      genalphaparams.set<double>("measured_curve1",iap.get<double>("MEASURED_CURVE1"));
      genalphaparams.set<double>("measured_curve2",iap.get<double>("MEASURED_CURVE2"));
      genalphaparams.set<double>("inv_ana_tol",iap.get<double>("INV_ANA_TOL"));

      // non-linear solution technique
      switch (Teuchos::getIntegralValue<int>(sdyn,"NLNSOL"))
      {
        case STRUCT_DYNAMIC::fullnewton:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
        case STRUCT_DYNAMIC::lsnewton:
          genalphaparams.set<string>("equilibrium iteration","line search newton");
        break;
        case STRUCT_DYNAMIC::modnewton:
          genalphaparams.set<string>("equilibrium iteration","modified newton");
        break;
        case STRUCT_DYNAMIC::nlncg:
          genalphaparams.set<string>("equilibrium iteration","nonlinear cg");
        break;
        case STRUCT_DYNAMIC::ptc:
          genalphaparams.set<string>("equilibrium iteration","ptc");
        break;
        case STRUCT_DYNAMIC::newtonlinuzawa:
          genalphaparams.set<string>("equilibrium iteration","newtonlinuzawa");
        break;
        case STRUCT_DYNAMIC::augmentedlagrange:
          genalphaparams.set<string>("equilibrium iteration","augmentedlagrange");              
        break;
        default:
          genalphaparams.set<string>("equilibrium iteration","full newton");
        break;
      }

      // set predictor (takes values "constant" "consistent")
      switch (Teuchos::getIntegralValue<int>(sdyn,"PREDICT"))
      {
        case STRUCT_DYNAMIC::pred_vague:
          dserror("You have to define the predictor");
          break;
        case STRUCT_DYNAMIC::pred_constdis:
          genalphaparams.set<string>("predictor","consistent");
          break;
        case STRUCT_DYNAMIC::pred_constdisvelacc:
          genalphaparams.set<string>("predictor","constant");
          break;
        default:
          dserror("Cannot cope with choice of predictor");
          break;
      }

      // detect if contact is present
      bool contact = false;
      switch (Teuchos::getIntegralValue<int>(scontact,"CONTACT"))
      {
        case INPUTPARAMS::contact_none:
          contact = false;
          break;
        case INPUTPARAMS::contact_normal:
          contact = true;
          break;
        case INPUTPARAMS::contact_frictional:
          contact = true;
          break;
        case INPUTPARAMS::contact_meshtying:
          contact = true;
          break;
        default:
          dserror("Cannot cope with choice of contact type");
          break;
      }
      
      // detect whether thermal bath is present
      bool thermalbath = false;
      switch (Teuchos::getIntegralValue<int>(statmech,"THERMALBATH"))
      {
        case INPUTPARAMS::thermalbath_none:
          thermalbath = false;
          break;
        case INPUTPARAMS::thermalbath_uniform:
          thermalbath = true;
          break;
        case INPUTPARAMS::thermalbath_shearflow:
          thermalbath = true;
          break;
        default:
          dserror("Cannot cope with choice of thermal bath");
          break;
      }
     
      // brownain motion
      genalphaparams.set<bool>  ("bro_motion",Teuchos::getIntegralValue<int>(bromop,"BROWNIAN_MOTION"));
      bool bromotion = genalphaparams.get("bro_motion",false);
      
      // create the time integrator
      bool inv_analysis = genalphaparams.get("inv_analysis",false);
      RCP<StruGenAlpha> tintegrator;
      if (!contact && !inv_analysis && !thermalbath && !bromotion)
        tintegrator = rcp(new StruGenAlpha(genalphaparams,*actdis,solver,output));
      else
      {
        if (contact)
          tintegrator = rcp(new CONTACT::ContactStruGenAlpha(genalphaparams,*actdis,solver,output));
        if (inv_analysis)
          tintegrator = rcp(new Inv_analysis(genalphaparams,*actdis,solver,output));
        if (thermalbath)
          tintegrator = rcp(new StatMechTime(genalphaparams,*actdis,solver,output));
        if (bromotion)
          tintegrator = rcp(new BroMotion_TimeInt(genalphaparams,*actdis,solver,output));
      }

      // do restart if demanded from input file
      // note that this changes time and step in genalphaparams
      if (genprob.restart)
        tintegrator->ReadRestart(genprob.restart);

      // write mesh always at beginning of calc or restart
      {
        int    step = genalphaparams.get<int>("step",0);
        double time = genalphaparams.get<double>("total time",0.0);
        output.WriteMesh(step,time);
      }

      // integrate in time and space
      tintegrator->Integrate();

      // test results
      {
        DRT::ResultTestManager testmanager(actdis->Comm());
        testmanager.AddFieldTest(rcp(new StruResultTest(*tintegrator)));
        testmanager.TestAll();
      }

    }
    break;
    //==================================================================
    // Generalized Energy Momentum Method
    //==================================================================
    case STRUCT_DYNAMIC::Gen_EMM :
    {
      dserror("Not yet impl.");
    }
    break;
    //==================================================================
    // Everything else
    //==================================================================
    default :
    {
      dserror("Time integration scheme is not available");
    }
    break;
  } // end of switch(sdyn->Typ)

  return;
} // end of dyn_nlnstructural_drt()

#endif  // #ifdef CCADISCRET
