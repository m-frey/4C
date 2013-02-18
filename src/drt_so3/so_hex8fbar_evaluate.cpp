/*!----------------------------------------------------------------------
\file so_hex8fbar_evaluate.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "Epetra_SerialDenseSolver.h"
#include "so_hex8fbar.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "prestress.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbar::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex8fbar::ActionType act = So_hex8fbar::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_hex8fbar::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_hex8fbar::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_hex8fbar::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_hex8fbar::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_hex8fbar::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_hex8fbar::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_hex8fbar::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_hex8fbar::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                         act = So_hex8fbar::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")                    act = So_hex8fbar::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")                  act = So_hex8fbar::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")                     act = So_hex8fbar::calc_struct_reset_istep;
  else if (action=="calc_struct_reset_discretization")            act = So_hex8fbar::calc_struct_reset_discretization;
  else if (action=="postprocess_stress")                          act = So_hex8fbar::postprocess_stress;
  else if (action=="multi_readrestart")                           act = So_hex8fbar::multi_readrestart;
  else if (action=="multi_calc_dens")                             act = So_hex8fbar::multi_calc_dens;
  else if (action=="calc_struct_prestress_update")                act = So_hex8fbar::prestress_update;
  else dserror("Unknown type of action for So_hex8fbar");
  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      std::vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      soh8fbar_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> myemat(true);
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        RCP<std::vector<char> > stressdata = params.get<RCP<std::vector<char> > >("stress",Teuchos::null);
        RCP<std::vector<char> > straindata = params.get<RCP<std::vector<char> > >("strain",Teuchos::null);
        if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        soh8fbar_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);
        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
        }
      }
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);

      Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
      if (poststress==Teuchos::null)
        dserror("No element stress/strain vector available");

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        soh8_expol(gpstress, *poststress);
      }
      else if (stresstype=="cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_SOH8; ++i)
          {
            double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH8;
          }
        }
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      // Update of history for plastic material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
      {
        MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
        plastic->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_growth)
      {
        MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
        grow->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
      {
        MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
        comix->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

    case calc_struct_update_imrlike:
    {
      // Update of history for plastic material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
      {
        MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
        plastic->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_growth)
      {
        MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
        grow->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
      {
        MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
        comix->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

    case calc_struct_reset_istep:
    {
      // Update of history for plastic material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
      {
	MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
	plastic->Update();
      }
    }
    break;

    //==================================================================================
    case calc_struct_reset_discretization:
    {
      // Reset of history for materials
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
      {
        MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
        comix->SetupHistory(NUMGPT_SOH8);
      }
      // Reset prestress
      if (pstype_==INPAR::STR::prestress_mulf)
      {
        time_ = 0.0;
        LINALG::Matrix<3,3> Id(true);
        Id(0,0) = Id(1,1) = Id(2,2) = 1.0;
        for (int gp=0; gp<NUMGPT_SOH8; ++gp)
        {
          prestress_->MatrixtoStorage(gp,Id,prestress_->FHistory());
          prestress_->MatrixtoStorage(gp,invJ_[gp],prestress_->JHistory());
        }
        prestress_->MatrixtoStorage(NUMGPT_SOH8,Id,prestress_->FHistory());
        LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;
        for (int i=0; i<NUMNOD_SOH8; ++i)
        {
          xrefe(i,0) = Nodes()[i]->X()[0];
          xrefe(i,1) = Nodes()[i]->X()[1];
          xrefe(i,2) = Nodes()[i]->X()[2];
        }
        LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
        DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ_0;
        invJ_0.Multiply(N_rst_0,xrefe);
        invJ_0.Invert();
        prestress_->MatrixtoStorage(NUMGPT_SOH8,invJ_0,prestress_->JHistory());
      }
      if (pstype_==INPAR::STR::prestress_id)
        dserror("Reset of Inverse Design not yet implemented");
    }
    break;

    case multi_calc_dens:
    {
      soh8_homog(params);
    }
    break;

    //==================================================================================
    case prestress_update:
    {
      time_ = params.get<double>("total time");
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // build def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOH8+1,9);
      DefGradient(mydisp,gpdefgrd,*prestress_);

      // update deformation gradient and put back to storage
      LINALG::Matrix<3,3> deltaF;
      LINALG::Matrix<3,3> Fhist;
      LINALG::Matrix<3,3> Fnew;
      for (int gp=0; gp<NUMGPT_SOH8+1; ++gp)
      {
        prestress_->StoragetoMatrix(gp,deltaF,gpdefgrd);
        prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());
        Fnew.Multiply(deltaF,Fhist);
        prestress_->MatrixtoStorage(gp,Fnew,prestress_->FHistory());
      }

      // push-forward invJ for every gaussian point
      UpdateJacobianMapping(mydisp,*prestress_);
    }
    break;

    //==================================================================================
    case inversedesign_update:
      dserror("The sohex8fbar element does not support inverse design analysis");
    break;


    // read restart of microscale
    case multi_readrestart:
    {
      RCP<MAT::Material> mat = Material();

      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
        soh8_read_restart_multi();
    }
    break;

    default:
      dserror("Unknown type of action for So_hex8fbar");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 03/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbar::InitJacobianMapping()
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    //invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(derivs[gp],xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0",detJ_[gp]);

    if (pstype_==INPAR::STR::prestress_mulf && pstime_ >= time_)
      if (!(prestress_->IsInit()))
        prestress_->MatrixtoStorage(gp,invJ_[gp],prestress_->JHistory());

  }

  // init the centroid invJ
  if (pstype_==INPAR::STR::prestress_mulf && pstime_ >= time_)
    if (!(prestress_->IsInit()))
    {
      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
      DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ_0;
      invJ_0.Multiply(N_rst_0,xrefe);
      invJ_0.Invert();
      prestress_->MatrixtoStorage(NUMGPT_SOH8,invJ_0,prestress_->JHistory());
    }


  if (pstype_==INPAR::STR::prestress_mulf && pstime_ >= time_)
    prestress_->IsInit() = true;

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)               |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbar::EvaluateNeumann(Teuchos::ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  // **

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int> >("funct");
  LINALG::Matrix<NUMDIM_SOH8,1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim=0; dim<NUMDIM_SOH8; dim++)
      if ((*funct)[dim] > 0)
	    havefunct = havefunct or true;

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
   LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct) {
      for (int dim=0; dim<NUMDIM_SOH8; dim++) {
        xrefegp(dim) = 0.0;
        for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid,dim);
      }
    }

    // integration factor
    const double fac = gpweights[gp] * curvefac * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<NUMDIM_SOH8; dim++) {
      // function evaluation
      const int functnum = (funct) ? (*funct)[dim] : -1;
      const double functfac
        = (functnum>0)
        ? DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,xrefegp.A(),time,NULL)
        : 1.0;
      const double dim_fac = (*onoff)[dim] * (*val)[dim] * fac * functfac;
      for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
        elevec1[nodid*NUMDIM_SOH8+dim] += shapefcts[gp](nodid) * dim_fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex8fbar::EvaluateNeumann

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbar::soh8fbar_nlnstiffmass(
      std::vector<int>&         lm,             // location matrix
      std::vector<double>&      disp,           // current displacements
      std::vector<double>&      residual,       // current residual displ
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix, // element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
      LINALG::Matrix<NUMDOF_SOH8,1>* force,                 // element internal force vector
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,   // stresses at GP
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,   // strains at GP
      Teuchos::ParameterList&   params,         // algorithmic parameters e.g. time
      const INPAR::STR::StressType   iostress,  // stress output option
      const INPAR::STR::StrainType   iostrain)  // strain output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];

    if (pstype_==INPAR::STR::prestress_mulf)
    {
      xdisp(i,0) = disp[i*NODDOF_SOH8+0];
      xdisp(i,1) = disp[i*NODDOF_SOH8+1];
      xdisp(i,2) = disp[i*NODDOF_SOH8+2];
    }
  }

  //****************************************************************************
	// deformation gradient at centroid of element
  //****************************************************************************
  double detF_0 = -1.0;
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd_0;
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ_0;
  //element coordinate derivatives at centroid
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
  DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);
  {
    //inverse jacobian matrix at centroid
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ_0;
    invJ_0.Multiply(N_rst_0,xrefe);
    invJ_0.Invert();
    //material derivatives at centroid
    N_XYZ_0.Multiply(invJ_0,N_rst_0);
  }

  if (pstype_==INPAR::STR::prestress_mulf)
  {
    // get Jacobian mapping wrt to the stored configuration
    // centroid is 9th Gaussian point in storage
    LINALG::Matrix<3,3> invJdef_0;
    prestress_->StoragetoMatrix(NUMGPT_SOH8,invJdef_0,prestress_->JHistory());
    // get derivatives wrt to last spatial configuration
    LINALG::Matrix<3,8> N_xyz_0;
    N_xyz_0.Multiply(invJdef_0,N_rst_0); //if (!Id()) cout << invJdef_0;

    // build multiplicative incremental defgrd
    LINALG::Matrix<3,3> defgrd_0(false);
    defgrd_0.MultiplyTT(xdisp,N_xyz_0);
    defgrd_0(0,0) += 1.0;
    defgrd_0(1,1) += 1.0;
    defgrd_0(2,2) += 1.0;

    // get stored old incremental F
    LINALG::Matrix<3,3> Fhist;
    prestress_->StoragetoMatrix(NUMGPT_SOH8,Fhist,prestress_->FHistory());

    // build total defgrd = delta F * F_old
    LINALG::Matrix<3,3> tmp;
    tmp.Multiply(defgrd_0,Fhist);
    defgrd_0 = tmp;

    // build inverse and detF
    invdefgrd_0.Invert(defgrd_0);
    detF_0=defgrd_0.Determinant();
  }
  else // no prestressing
  {
    //deformation gradient and its determinant at centroid
    LINALG::Matrix<3,3> defgrd_0(false);
    defgrd_0.MultiplyTT(xcurr,N_XYZ_0);
    invdefgrd_0.Invert(defgrd_0);
    detF_0=defgrd_0.Determinant();
  }
  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(false);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]);
    double detJ = detJ_[gp];

    if (pstype_==INPAR::STR::prestress_mulf)
    {
      // get Jacobian mapping wrt to the stored configuration
      LINALG::Matrix<3,3> invJdef;
      prestress_->StoragetoMatrix(gp,invJdef,prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<3,8> N_xyz;
      N_xyz.Multiply(invJdef,derivs[gp]);

      // build multiplicative incremental defgrd
      defgrd.MultiplyTT(xdisp,N_xyz);
      defgrd(0,0) += 1.0;
      defgrd(1,1) += 1.0;
      defgrd(2,2) += 1.0;

      // get stored old incremental F
      LINALG::Matrix<3,3> Fhist;
      prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());

      // build total defgrd = delta F * F_old
      LINALG::Matrix<3,3> Fnew;
      Fnew.Multiply(defgrd,Fhist);
      defgrd = Fnew;
    }
    else // no prestressing
    {
      // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
      defgrd.MultiplyTT(xcurr,N_XYZ);
    }
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd;
    invdefgrd.Invert(defgrd);
    double detF=defgrd.Determinant();

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // F_bar deformation gradient =(detF_0/detF)^1/3*F
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd_bar(defgrd);
    double f_bar_factor=pow(detF_0/detF,1.0/3.0);
    defgrd_bar.Scale(f_bar_factor);

    // Right Cauchy-Green tensor(Fbar) = F_bar^T * F_bar
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen_bar;
    cauchygreen_bar.MultiplyTN(defgrd_bar,defgrd_bar);

    // Green-Lagrange strains(F_bar) matrix E = 0.5 * (Cauchygreen(F_bar) - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_bar_epetra(NUMSTR_SOH8);
    LINALG::Matrix<NUMSTR_SOH8,1> glstrain_bar(glstrain_bar_epetra.A(),true);
    glstrain_bar(0) = 0.5 * (cauchygreen_bar(0,0) - 1.0);
    glstrain_bar(1) = 0.5 * (cauchygreen_bar(1,1) - 1.0);
    glstrain_bar(2) = 0.5 * (cauchygreen_bar(2,2) - 1.0);
    glstrain_bar(3) = cauchygreen_bar(0,1);
    glstrain_bar(4) = cauchygreen_bar(1,2);
    glstrain_bar(5) = cauchygreen_bar(2,0);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      for (int i = 0; i < 3; ++i)
        (*elestrain)(gp,i) = glstrain_bar(i);
      for (int i = 3; i < 6; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain_bar(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl_bar;
      gl_bar(0,0) = glstrain_bar(0);
      gl_bar(0,1) = 0.5*glstrain_bar(3);
      gl_bar(0,2) = 0.5*glstrain_bar(5);
      gl_bar(1,0) = gl_bar(0,1);
      gl_bar(1,1) = glstrain_bar(1);
      gl_bar(1,2) = 0.5*glstrain_bar(4);
      gl_bar(2,0) = gl_bar(0,2);
      gl_bar(2,1) = gl_bar(1,2);
      gl_bar(2,2) = glstrain_bar(2);

      // inverse of fbar deformation gradient
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd_bar;
      invdefgrd_bar.Invert(defgrd_bar);

      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi_bar;
      temp.Multiply(gl_bar,invdefgrd_bar);
      euler_almansi_bar.MultiplyTN(invdefgrd_bar,temp);

      (*elestrain)(gp,0) = euler_almansi_bar(0,0);
      (*elestrain)(gp,1) = euler_almansi_bar(1,1);
      (*elestrain)(gp,2) = euler_almansi_bar(2,2);
      (*elestrain)(gp,3) = euler_almansi_bar(0,1);
      (*elestrain)(gp,4) = euler_almansi_bar(1,2);
      (*elestrain)(gp,5) = euler_almansi_bar(0,2);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not available");
    }

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      bop(0,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stress_bar(true);
    LINALG::Matrix<NUMSTR_SOH8,1> plglstrain(true);
    soh8_mat_sel(&stress_bar,&cmat,&density,&glstrain_bar,&plglstrain,&defgrd_bar,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH8; ++i)
        (*elestress)(gp,i) = stress_bar(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF_bar = defgrd_bar.Determinant();

      LINALG::Matrix<3,3> pkstress_bar;
      pkstress_bar(0,0) = stress_bar(0);
      pkstress_bar(0,1) = stress_bar(3);
      pkstress_bar(0,2) = stress_bar(5);
      pkstress_bar(1,0) = pkstress_bar(0,1);
      pkstress_bar(1,1) = stress_bar(1);
      pkstress_bar(1,2) = stress_bar(4);
      pkstress_bar(2,0) = pkstress_bar(0,2);
      pkstress_bar(2,1) = pkstress_bar(1,2);
      pkstress_bar(2,2) = stress_bar(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystress_bar;
      temp.Multiply(1.0/detF_bar,defgrd_bar,pkstress_bar);
      cauchystress_bar.MultiplyNT(temp,defgrd_bar);

      (*elestress)(gp,0) = cauchystress_bar(0,0);
      (*elestress)(gp,1) = cauchystress_bar(1,1);
      (*elestress)(gp,2) = cauchystress_bar(2,2);
      (*elestress)(gp,3) = cauchystress_bar(0,1);
      (*elestress)(gp,4) = cauchystress_bar(1,2);
      (*elestress)(gp,5) = cauchystress_bar(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
    }

    double detJ_w = detJ*gpweights[gp];

    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w/f_bar_factor, bop, stress_bar, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH8> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w*f_bar_factor,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(stress_bar); // auxiliary integrated stress
      sfac.Scale(detJ_w/f_bar_factor); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************

      // integrate additional fbar matrix
     LINALG::Matrix<NUMSTR_SOH8,1> cauchygreenvector;
     cauchygreenvector(0) = cauchygreen(0,0);
     cauchygreenvector(1) = cauchygreen(1,1);
     cauchygreenvector(2) = cauchygreen(2,2);
     cauchygreenvector(3) = 2*cauchygreen(0,1);
     cauchygreenvector(4) = 2*cauchygreen(1,2);
     cauchygreenvector(5) = 2*cauchygreen(2,0);

     LINALG::Matrix<NUMSTR_SOH8,1> ccg;
     ccg.Multiply(cmat,cauchygreenvector);

     LINALG::Matrix<NUMDOF_SOH8,1> bopccg(false); // auxiliary integrated stress
     bopccg.MultiplyTN(detJ_w*f_bar_factor/3.0,bop,ccg);

     double htensor[NUMDOF_SOH8];
     for(int n=0;n<NUMDOF_SOH8;n++)
     {
       htensor[n]=0;
       for(int i=0;i<NUMDIM_SOH8;i++)
       {
         htensor[n] += invdefgrd_0(i,n%3)*N_XYZ_0(i,n/3)-invdefgrd(i,n%3)*N_XYZ(i,n/3);
       }
     }

     LINALG::Matrix<NUMDOF_SOH8,1> bops(false); // auxiliary integrated stress
     bops.MultiplyTN(-detJ_w/f_bar_factor/3.0,bop,stress_bar);
     for(int i=0;i<NUMDOF_SOH8;i++)
     {
       for (int j=0;j<NUMDOF_SOH8;j++)
       {
         (*stiffmatrix)(i,j) += htensor[j]*(bops(i,0)+bopccg(i,0));
       }
     } // end of integrate additional `fbar' stiffness**********************
    }  // if (stiffmatrix != NULL)

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  }/* ==================================================== end of Loop over GP */

  return;
} // DRT::ELEMENTS::So_hex8fbar::soh8fbar_nlnstiffmass

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex8fbar* actele = dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8fbar* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbar::DefGradient(const std::vector<double>& disp,
                                         Epetra_SerialDenseMatrix& gpdefgrd,
                                         DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  // derivatives at centroid point
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
  DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOH8+0];
    xdisp(i,1) = disp[i*NODDOF_SOH8+1];
    xdisp(i,2) = disp[i*NODDOF_SOH8+2];
  }

  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::Matrix<3,3> invJdef;
    prestress.StoragetoMatrix(gp,invJdef,prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_xyz;
    N_xyz.Multiply(invJdef,derivs[gp]);

    // build defgrd (independent of xrefe!)
    LINALG::Matrix<3,3> defgrd;
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;

    prestress.MatrixtoStorage(gp,defgrd,gpdefgrd);
  }

  {
    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::Matrix<3,3> invJdef;
    prestress.StoragetoMatrix(NUMGPT_SOH8,invJdef,prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_xyz;
    N_xyz.Multiply(invJdef,N_rst_0);

    // build defgrd (independent of xrefe!)
    LINALG::Matrix<3,3> defgrd;
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;

    prestress.MatrixtoStorage(NUMGPT_SOH8,defgrd,gpdefgrd);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbar::UpdateJacobianMapping(
                                            const std::vector<double>& disp,
                                            DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  // derivatives at centroid
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
  DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);

  // get incremental disp
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOH8+0];
    xdisp(i,1) = disp[i*NODDOF_SOH8+1];
    xdisp(i,2) = disp[i*NODDOF_SOH8+2];
  }

  LINALG::Matrix<3,3> invJhist;
  LINALG::Matrix<3,3> invJ;
  LINALG::Matrix<3,3> defgrd;
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_xyz;
  LINALG::Matrix<3,3> invJnew;
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp,invJhist,prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist,derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd,invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp,invJnew,prestress.JHistory());
  } // for (int gp=0; gp<NUMGPT_SOH8; ++gp)

  {
    // get the invJ old state
    prestress.StoragetoMatrix(NUMGPT_SOH8,invJhist,prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist,N_rst_0);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.MultiplyTT(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd,invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(NUMGPT_SOH8,invJnew,prestress.JHistory());
  }

  return;
}


