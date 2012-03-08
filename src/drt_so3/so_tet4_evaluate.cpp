/*!----------------------------------------------------------------------*
\file so_tet4_evaluate.cpp
\brief quadratic nonlinear tetrahedron

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet4.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_patspec/patspec.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mortar/mortar_analytical.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"

//#define PRINT_DEBUG
#ifdef PRINT_DEBUG
#include <string>
#include <sstream>
#include <cstring>
template <class T>
void writeArray(const T& mat, std::string name = "unnamed")
{
  std::stringstream header;
  header << 'M' << name << ':' << mat.M() << 'x' << mat.N() << ':';
  unsigned int s = header.str().size() + mat.M()*mat.N()*sizeof(double);
  std::cerr.write(reinterpret_cast<const char*>(&s),sizeof(unsigned int));
  std::cerr << header.str();
  for (int i = 0; i < mat.M()*mat.N(); ++i) {
    std::cerr.write(reinterpret_cast<const char*>(&(mat.A()[i])),sizeof(double));
  }
}

void writeComment(const std::string v)
{
  unsigned int s = v.size()+1;
  std::cerr.write(reinterpret_cast<const char*>(&s),sizeof(unsigned int));
  std::cerr << 'C' << v;
}
#endif // PRINT_DEBUG


using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4::Evaluate(ParameterList&           params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET4,1>             elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET4,1>             elevec2(elevec2_epetra.A(),true);

  // start with "none"
  DRT::ELEMENTS::So_tet4::ActionType act = So_tet4::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")             act = So_tet4::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")             act = So_tet4::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")        act = So_tet4::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")         act = So_tet4::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")         act = So_tet4::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")        act = So_tet4::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")               act = So_tet4::calc_struct_stress;
  else if (action=="postprocess_stress")               act = So_tet4::postprocess_stress;
  else if (action=="calc_struct_eleload")              act = So_tet4::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")              act = So_tet4::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")         act = So_tet4::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")       act = So_tet4::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")          act = So_tet4::calc_struct_reset_istep;
  else if (action=="calc_struct_reset_discretization") act = So_tet4::calc_struct_reset_discretization;
  else if (action=="calc_struct_errornorms")           act = So_tet4::calc_struct_errornorms;
  else if (action=="calc_struct_prestress_update")     act = So_tet4::prestress_update;
  else if (action=="calc_struct_energy")	       act = So_tet4::calc_struct_energy;
  else if (action=="calc_struct_output_E")	       act = So_tet4::calc_struct_output_E; 
  else if (action=="calc_struct_inversedesign_update") act = So_tet4::inversedesign_update;
  else if (action=="calc_struct_inversedesign_switch") act = So_tet4::inversedesign_switch;
  else if (action=="multi_calc_dens")                  act = So_tet4::multi_calc_dens;
  else if (action=="multi_readrestart")                act = So_tet4::multi_readrestart;
  else dserror("Unknown type of action for So_tet4");

  // check for patient specific data
  PATSPEC::GetILTDistance(Id(),params,discretization);
  PATSPEC::GetLocalRadius(Id(),params,discretization);

  // get the material law
  Teuchos::RCP<MAT::Material> actmat = Material();

  // what should the element do
  switch(act)
  {
    //==================================================================================
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      so_tet4_nlnstiffmass(params,lm,mydisp,myres, &elemat1, NULL, &elevec1, NULL,NULL,actmat,
                           INPAR::STR::stress_none,INPAR::STR::strain_none);
   }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      if (pstype_==INPAR::STR::prestress_id && time_ <= pstime_) // inverse design analysis
        invdesign_->so_tet4_nlnstiffmass(params,this,lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat,
                                         INPAR::STR::stress_none,INPAR::STR::strain_none);
      else
        so_tet4_nlnstiffmass(params,lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat,
                             INPAR::STR::stress_none,INPAR::STR::strain_none);

    }
    break;

    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4> myemat(true); // to zero
      so_tet4_nlnstiffmass(params,lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,actmat,
                           INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      if (pstype_==INPAR::STR::prestress_id && time_ <= pstime_) // inverse design analysis
        invdesign_->so_tet4_nlnstiffmass(params,this,lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,
                                         INPAR::STR::stress_none,INPAR::STR::strain_none);
      else
        so_tet4_nlnstiffmass(params,lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,
                             INPAR::STR::stress_none,INPAR::STR::strain_none);

      if (act==calc_struct_nlnstifflmass) so_tet4_lumpmass(&elemat2);
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
        RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==null) dserror("Cannot get 'stress' data");
        if (straindata==null) dserror("Cannot get 'strain' data");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4> stress(true); // set to zero
        LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4> strain(true);
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        if (pstype_==INPAR::STR::prestress_id && time_ <= pstime_) // inverse design analysis
          invdesign_->so_tet4_nlnstiffmass(params,this,lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,iostress,iostrain);
        else
          so_tet4_nlnstiffmass(params,lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,iostress,iostrain);

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

    //==================================================================================
    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
          params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
        if (gpstressmap==null) dserror("no gp stress/strain map available for postprocessing");

        string stresstype = params.get<string>("stresstype","ndxyz");
        int gid = Id();
        LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4> gpstress(((*gpstressmap)[gid])->A(),true);

        RCP<Epetra_MultiVector> poststress=params.get<RCP<Epetra_MultiVector> >("poststress",null);
        if (poststress==null) dserror("No element stress/strain vector available");

        if (stresstype=="ndxyz")
        {
          // extrapolate stresses/strains at Gauss points to nodes
          so_tet4_expol(gpstress, *poststress);
        }
        else if (stresstype=="cxyz")
        {
          const Epetra_BlockMap elemap = poststress->Map();
          int lid = elemap.LID(Id());
          if (lid!=-1)
          {
            for (int i = 0; i < NUMSTR_SOTET4; ++i)
            {
              double& s = (*((*poststress)(i)))[lid];
              s = 0.;
              for (int j = 0; j < NUMGPT_SOTET4; ++j)
              {
                s += gpstress(j,i);
              }
              s /= NUMGPT_SOTET4;
            }
          }
        }
        else
          dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    //==================================================================================
    case prestress_update:
    {
      time_ = params.get<double>("total time");
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // build incremental def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET4,9);
      DefGradient(mydisp,gpdefgrd,*prestress_);

      // update deformation gradient and put back to storage
      LINALG::Matrix<3,3> deltaF;
      LINALG::Matrix<3,3> Fhist;
      LINALG::Matrix<3,3> Fnew;
      for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
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
    // this is a dummy output for strain energy
    case calc_struct_energy:
    {
    // check length of elevec1
    if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");
    RCP<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_stvenant)
    {
    	double intenergy = -1.0;
        elevec1_epetra(0) = intenergy;
    }  else dserror("No dummy strain energy for material other than INPAR::MAT::m_stvenant");
    }
    break;
    // this is needed by bone topology optimization
    case calc_struct_output_E:
    {
     RCP<MAT::Material> mat = Material();
      // check length of elevec1
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");
      MAT::Material* rawmat = mat.get();
      MAT::StVenantKirchhoff* stvk = dynamic_cast<MAT::StVenantKirchhoff*>(rawmat);
      if (!stvk) dserror("dynamic cast to stvenant failed");
      double E = stvk->Youngs();
      elevec1_epetra(0) = E;
    
    }
    break;

    //==================================================================================
    case inversedesign_update:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      invdesign_->sot4_StoreMaterialConfiguration(this,mydisp);
      invdesign_->IsInit() = true; // this is to make the restart work
    }
    break;

    //==================================================================================
    case inversedesign_switch:
    {
      time_ = params.get<double>("total time");
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    //==================================================================================
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      // determine new fiber directions
      RCP<MAT::Material> mat = Material();
      bool remodel;
      const Teuchos::ParameterList& patspec = DRT::Problem::Instance()->PatSpecParams();
      remodel = DRT::INPUT::IntegralValue<int>(patspec,"REMODEL");
      if (remodel &&
          ((mat->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular) ||
           (mat->MaterialType() == INPAR::MAT::m_humphreycardiovascular) ||
           (mat->MaterialType() == INPAR::MAT::m_constraintmixture)))// && timen_ <= timemax_ && stepn_ <= stepmax_)
      {
        RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        so_tet4_remodel(lm,mydisp,params,mat);
      }
      // Update of history for visco material
      if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
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

    //==================================================================================
    case calc_struct_update_imrlike:
    {
      // determine new fiber directions
      RCP<MAT::Material> mat = Material();
      bool remodel;
      const Teuchos::ParameterList& patspec = DRT::Problem::Instance()->PatSpecParams();
      remodel = DRT::INPUT::IntegralValue<int>(patspec,"REMODEL");
      if (remodel &&
          ((mat->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular) ||
           (mat->MaterialType() == INPAR::MAT::m_humphreycardiovascular) ||
           (mat->MaterialType() == INPAR::MAT::m_constraintmixture)))// && timen_ <= timemax_ && stepn_ <= stepmax_)
      {
        RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        so_tet4_remodel(lm,mydisp,params,mat);
      }
      // Update of history for visco material
      if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
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

    //==================================================================================
    case calc_struct_reset_istep:
    {
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
      {
        MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
        comix->Reset();
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
        comix->SetupHistory(NUMGPT_SOTET4);
      }
      // Reset prestress
      if (pstype_==INPAR::STR::prestress_mulf)
      {
        time_ = 0.0;
        LINALG::Matrix<3,3> Id(true);
        Id(0,0) = Id(1,1) = Id(2,2) = 1.0;
        for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
        {
          prestress_->MatrixtoStorage(gp,Id,prestress_->FHistory());
          prestress_->MatrixtoStorage(gp,nxyz_,prestress_->JHistory());
        }
      }
      if (pstype_==INPAR::STR::prestress_id)
        dserror("Reset of Inverse Design not yet implemented");
    }
    break;

    //==================================================================================
    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    //==================================================================================
    case calc_struct_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - 4 Gauss point rule is used for integration of error norms
      // - only implemented for SVK material (relevant for energy norm only, L2 and
      //   H1 norms are of course valid for arbitrary materials)
      // - analytical solutions are currently stored in a repository in the MORTAR
      //   namespace, however they could (should?) be moved to a more general location

      // check length of elevec1
      if (elevec1_epetra.Length() < 3) dserror("The given result vector is too short.");

      // check material law
      RCP<MAT::Material> mat = Material();

      //******************************************************************
      // only for St.Venant Kirchhoff material
      //******************************************************************
      if (mat->MaterialType() == INPAR::MAT::m_stvenant)
      {
        // declaration of variables
        double l2norm = 0.0;
        double h1norm = 0.0;
        double energynorm = 0.0;

        // use 4 Gauss points for integration, not only 1
        int ngp = 4;

        // shape functions, derivatives and integration weights
        const static vector<LINALG::Matrix<NUMNOD_SOTET4,1> > vals = so_tet4_4gp_shapefcts();
        const static std::vector<double> weights = so_tet4_4gp_weights();

        // get displacements and extract values of this element
        RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state displacement vector");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

        // nodal displacement vector
        LINALG::Matrix<NUMDOF_SOTET4,1> nodaldisp;
        for (int i=0; i<NUMDOF_SOTET4; ++i) nodaldisp(i,0) = mydisp[i];

        // reference geometry (nodal positions)
        LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xrefe;
        DRT::Node** nodes = Nodes();
        for (int i=0; i<NUMNOD_SOTET4; ++i)
        {
          xrefe(i,0) = nodes[i]->X()[0];
          xrefe(i,1) = nodes[i]->X()[1];
          xrefe(i,2) = nodes[i]->X()[2];
        }

        // deformation gradient = identity tensor (geometrically linear case!)
        LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> defgrd(true);
        for (int i=0;i<NUMDIM_SOTET4;++i) defgrd(i,i) = 1;

        //----------------------------------------------------------------
        // loop over all Gauss points
        //----------------------------------------------------------------
        for (int gp=0; gp<ngp; gp++)
        {
          // Gauss weights and Jacobian determinant
          double fac = V_ * weights[gp];

          // Gauss point in reference configuration
          LINALG::Matrix<NUMDIM_SOTET4,1> xgp(true);
          for (int k=0;k<NUMDIM_SOTET4;++k)
            for (int n=0;n<NUMNOD_SOTET4;++n)
              xgp(k,0) += (vals[gp])(n) * xrefe(n,k);

          //**************************************************************
          // get analytical solution
          LINALG::Matrix<NUMDIM_SOTET4,1> uanalyt(true);
          LINALG::Matrix<NUMSTR_SOTET4,1> strainanalyt(true);
          LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> derivanalyt(true);

          MORTAR::AnalyticalSolutions3D(xgp,uanalyt,strainanalyt,derivanalyt);
          //**************************************************************

          //--------------------------------------------------------------
          // (1) L2 norm
          //--------------------------------------------------------------

          // compute displacements at GP
          LINALG::Matrix<NUMDIM_SOTET4,1> ugp(true);
          for (int k=0;k<NUMDIM_SOTET4;++k)
            for (int n=0;n<NUMNOD_SOTET4;++n)
              ugp(k,0) += (vals[gp])(n) * nodaldisp(NODDOF_SOTET4*n+k,0);

          // displacement error
          LINALG::Matrix<NUMDIM_SOTET4,1> uerror(true);
          for (int k=0;k<NUMDIM_SOTET4;++k)
            uerror(k,0) = uanalyt(k,0) - ugp(k,0);

          // compute GP contribution to L2 error norm
          l2norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (2) H1 norm
          //--------------------------------------------------------------

          // compute derivatives N_XYZ at GP w.r.t. material coordinates
          // (we have the transposed stored in nxyz_ and then switch indices)
          LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> tr_N_XYZ(nxyz_); // copy!
          LINALG::Matrix<NUMDIM_SOTET4,NUMNOD_SOTET4> N_XYZ(true);
          for (int k=0;k<NUMNOD_SOTET4;++k)
            for (int m=0;m<NUMDIM_SOTET4;++m)
              N_XYZ(m,k) = tr_N_XYZ(k,m);

          // compute partial derivatives at GP
          LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> derivgp(true);
          for (int l=0;l<NUMDIM_SOTET4;++l)
            for (int m=0;m<NUMDIM_SOTET4;++m)
              for (int k=0;k<NUMNOD_SOTET4;++k)
                derivgp(l,m) += N_XYZ(m,k) * nodaldisp(NODDOF_SOTET4*k+l,0);

          // derivative error
          LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> deriverror(true);
          for (int k=0;k<NUMDIM_SOTET4;++k)
            for (int m=0;m<NUMDIM_SOTET4;++m)
              deriverror(k,m) = derivanalyt(k,m) - derivgp(k,m);

          // compute GP contribution to H1 error norm
          h1norm += fac * deriverror.Dot(deriverror);
          h1norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (3) Energy norm
          //--------------------------------------------------------------

          // compute linear B-operator
          LINALG::Matrix<NUMSTR_SOTET4,NUMDOF_SOTET4> bop;
          for (int i=0; i<NUMNOD_SOTET4; ++i)
          {
            bop(0,NODDOF_SOTET4*i+0) = N_XYZ(0,i);
            bop(0,NODDOF_SOTET4*i+1) = 0.0;
            bop(0,NODDOF_SOTET4*i+2) = 0.0;
            bop(1,NODDOF_SOTET4*i+0) = 0.0;
            bop(1,NODDOF_SOTET4*i+1) = N_XYZ(1,i);
            bop(1,NODDOF_SOTET4*i+2) = 0.0;
            bop(2,NODDOF_SOTET4*i+0) = 0.0;
            bop(2,NODDOF_SOTET4*i+1) = 0.0;
            bop(2,NODDOF_SOTET4*i+2) = N_XYZ(2,i);

            bop(3,NODDOF_SOTET4*i+0) = N_XYZ(1,i);
            bop(3,NODDOF_SOTET4*i+1) = N_XYZ(0,i);
            bop(3,NODDOF_SOTET4*i+2) = 0.0;
            bop(4,NODDOF_SOTET4*i+0) = 0.0;
            bop(4,NODDOF_SOTET4*i+1) = N_XYZ(2,i);
            bop(4,NODDOF_SOTET4*i+2) = N_XYZ(1,i);
            bop(5,NODDOF_SOTET4*i+0) = N_XYZ(2,i);
            bop(5,NODDOF_SOTET4*i+1) = 0.0;
            bop(5,NODDOF_SOTET4*i+2) = N_XYZ(0,i);
          }

          // compute linear strain at GP
          LINALG::Matrix<NUMSTR_SOTET4,1> straingp(true);
          straingp.Multiply(bop,nodaldisp);

          // strain error
          LINALG::Matrix<NUMSTR_SOTET4,1> strainerror(true);
          for (int k=0;k<NUMSTR_SOTET4;++k)
            strainerror(k,0) = strainanalyt(k,0) - straingp(k,0);

          // compute stress vector and constitutive matrix
          double density = 0.0;
          LINALG::Matrix<NUMSTR_SOTET4,NUMSTR_SOTET4> cmat(true);
          LINALG::Matrix<NUMSTR_SOTET4,1> stress(true);
          so_tet4_mat_sel(&stress,&cmat,&density,&strainerror,&defgrd,gp,params);

          // compute GP contribution to energy error norm
          energynorm += fac * stress.Dot(strainerror);

          //cout << "UAnalytical:      " << ugp << endl;
          //cout << "UDiscrete:        " << uanalyt << endl;
          //cout << "StrainAnalytical: " << strainanalyt << endl;
          //cout << "StrainDiscrete:   " << straingp << endl;
          //cout << "DerivAnalytical:  " << derivanalyt << endl;
          //cout << "DerivDiscrete:    " << derivgp << endl;
        }
        //----------------------------------------------------------------

        // return results
        elevec1_epetra(0) = l2norm;
        elevec1_epetra(1) = h1norm;
        elevec1_epetra(2) = energynorm;
      }
      else
        dserror("ERROR: Error norms only implemented for SVK material");
    }
    break;

    //==================================================================================
    case multi_calc_dens:
    {
      sotet4_homog(params);
    }
    break;

    //==================================================================================
    // read restart of microscale
    case multi_readrestart:
    {
      sotet4_read_restart_multi();
    }
    break;

    default:
      dserror("Unknown type of action for so_tet4");
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  dserror("DRT::ELEMENTS::So_tet4::EvaluateNeumann not implemented");
#if 0
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  // **

/* =============================================================================*
 * CONST SHAPE FUNCTIONS and WEIGHTS for TET_4 with 1 GAUSS POINTS              *
 * =============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOTET4,1> > shapefcts = so_tet4_1gp_shapefcts();
  const static vector<double> gpweights = so_tet4_1gp_weights();
/* ============================================================================*/

/* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  {
    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1    1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **		   [ z_1  z_2  z_3  z_4 ]
    */
    LINALG::Matrix<NUMCOORD_SOTET4,NUMCOORD_SOTET4> jac_coord;
    for (int i=0; i<4; i++) jac_coord(0,i)=1;
    DRT::Node** nodes = Nodes();
    for (int col=0;col<4;col++)
    {
      const double* x = nodes[col]->X();
      for (int row=0;row<3;row++)
        jac_coord(row+1,col) = x[row];
    }

    // compute determinant of Jacobian with own algorithm
    // !!warning detJ is not the actual determinant of the jacobian (here needed for the quadrature rule)
    // but rather the volume of the tetrahedara
    double detJ=jac_coord.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = gpweights[gp] * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for(int dim=0; dim<NUMDIM_SOTET4; dim++) {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<NUMNOD_SOTET4; ++nodid) {
        elevec1(nodid*NUMDIM_SOTET4+dim) += shapefcts[gp](nodid) * dim_fac;
      }
    }
  }/* ==================================================== end of Loop over GP */
#endif
  return 0;
} // DRT::ELEMENTS::So_tet4::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::InitJacobianMapping()
{
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **		 [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::Matrix<NUMCOORD_SOTET4,NUMCOORD_SOTET4> jac;
  for (int i=0; i<4; i++)  jac(0,i)=1;
  for (int row=0;row<3;row++)
    for (int col=0;col<4;col++)
      jac(row+1,col)= xrefe(col,row);
  // volume of the element
  V_ = jac.Determinant()/6.0;
  if (V_ <= 0.0) dserror("Element volume %10.5e <= 0.0",V_);

  //nxyz_.resize(NUMGPT_SOTET4);
  const static vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > derivs = so_tet4_1gp_derivs();
  LINALG::Matrix<NUMCOORD_SOTET4-1,NUMCOORD_SOTET4> tmp;
  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    tmp.MultiplyTN(xrefe,derivs[gp]);
    for (int i=0; i<4; i++) jac(0,i)=1;
    for (int row=0;row<3;row++)
      for (int col=0;col<4;col++)
        jac(row+1,col)=tmp(row,col);
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4,NUMDIM_SOTET4> I_aug(true);
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4,NUMDIM_SOTET4> partials(true);
    I_aug(1,0)=1;
    I_aug(2,1)=1;
    I_aug(3,2)=1;

    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_SOTET4,NUMCOORD_SOTET4,NUMDIM_SOTET4> solve_for_inverseJac;
    solve_for_inverseJac.SetMatrix(jac);  // set A=jac
    solve_for_inverseJac.SetVectors(partials, I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) || (err2!=0))
    	dserror("Inversion of Jacobian failed");

    //nxyz_[gp] = N_xsi_k*partials
    nxyz_.Multiply(derivs[gp],partials);
    /* structure of N_XYZ:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    N_XYZ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_4     dN_4     dN_4   ]
    **             [  -------  -------  ------- ]
    **             [    dX       dY       dZ    ]
    */

    if (pstype_==INPAR::STR::prestress_mulf && pstime_ >= time_)
      if (!(prestress_->IsInit()))
        prestress_->MatrixtoStorage(gp,nxyz_,prestress_->JHistory());

    if (pstype_==INPAR::STR::prestress_id && pstime_ < time_)
      if (!(invdesign_->IsInit()))
      {
        invdesign_->MatrixtoStorage(gp,nxyz_,invdesign_->JHistory());
        invdesign_->DetJHistory()[gp] = V_;
      }

  } // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  if (pstype_==INPAR::STR::prestress_mulf && pstime_ >= time_)
    prestress_->IsInit() = true;

  if (pstype_==INPAR::STR::prestress_id && pstime_ < time_)
    invdesign_->IsInit() = true;

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 08/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_nlnstiffmass(
      ParameterList&            params,
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residual displ
      LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* stiffmatrix,    // element stiffness matrix
      LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* massmatrix,     // element mass matrix
      LINALG::Matrix<NUMDOF_SOTET4,1>* force,          // element internal force vector
      LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4>* elestress,      // stresses at GP
      LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4>* elestrain,      // strains at GP
      Teuchos::RCP<const MAT::Material>            material,       // element material data
      const INPAR::STR::StressType                 iostress,         // stress output options
      const INPAR::STR::StrainType                 iostrain)         // strain output options
{
/* =============================================================================*
** CONST DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
** =============================================================================*/
  const static vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > derivs = so_tet4_1gp_derivs();
  const static vector<double> gpweights = so_tet4_1gp_weights();
/* ============================================================================*/
  double density;
  // element geometry
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_4   Y_4   Z_4  ]
    */
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_4   y_4   z_4  ]
    */
  // current  displacements of element
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xdisp;
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }


  //volume of a tetrahedra
  double detJ = V_;

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  {
    LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> nxyz(nxyz_); // copy!

    //                                      d xcurr
    // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    // size is 3x3
    LINALG::Matrix<3,3> defgrd(false);

    if (pstype_==INPAR::STR::prestress_mulf)
    {
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> N_xyz;
      prestress_->StoragetoMatrix(gp,N_xyz,prestress_->JHistory());

      // build multiplicative incremental defgrd
      //defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
      defgrd.MultiplyTN(xdisp,N_xyz);
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
    else
    {
      defgrd.MultiplyTN(xdisp,nxyz);
      defgrd(0,0)+=1;
      defgrd(1,1)+=1;
      defgrd(2,2)+=1;
    }

    if (pstype_==INPAR::STR::prestress_id && pstime_ < time_)
    {
      // make the multiplicative update so that defgrd refers to
      // the reference configuration that resulted from the inverse
      // design analysis
      LINALG::Matrix<3,3> Fhist;
      invdesign_->StoragetoMatrix(gp,Fhist,invdesign_->FHistory());
      LINALG::Matrix<3,3> tmp3x3;
      tmp3x3.Multiply(defgrd,Fhist);
      defgrd = tmp3x3;

      // make detJ and nxyzmat refer to the ref. configuration that resulted from
      // the inverse design analysis
      detJ = invdesign_->DetJHistory()[gp];
      invdesign_->StoragetoMatrix(gp,nxyz,invdesign_->JHistory());
    }

    // Right Cauchy-Green tensor = F^T * F
    // size is 3x3
    LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::Matrix<6,1> glstrain(false);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("no strain data available");
      for (int i = 0; i < 3; ++i)
        (*elestrain)(gp,i) = glstrain(i);
      for (int i = 3; i < 6; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("no strain data available");

      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> gl;
      gl(0,0) = glstrain(0);
      gl(0,1) = 0.5*glstrain(3);
      gl(0,2) = 0.5*glstrain(5);
      gl(1,0) = gl(0,1);
      gl(1,1) = glstrain(1);
      gl(1,2) = 0.5*glstrain(4);
      gl(2,0) = gl(0,2);
      gl(2,1) = gl(1,2);
      gl(2,2) = glstrain(2);

      // inverse of deformation gradient
      //Epetra_SerialDenseMatrix invdefgrd(defgrd); // make a copy here otherwise defgrd is destroyed!
      //LINALG::NonsymInverse3x3(invdefgrd);
      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> invdefgrd;
      invdefgrd.Invert(defgrd);

      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> temp;
      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> euler_almansi;
      temp.Multiply(gl,invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd,temp);

      (*elestrain)(gp,0) = euler_almansi(0,0);
      (*elestrain)(gp,1) = euler_almansi(1,1);
      (*elestrain)(gp,2) = euler_almansi(2,2);
      (*elestrain)(gp,3) = euler_almansi(0,1);
      (*elestrain)(gp,4) = euler_almansi(1,2);
      (*elestrain)(gp,5) = euler_almansi(0,2);
    }
    break;
    case INPAR::STR::strain_none:
    break;
    default:
      dserror("requested strain option not available");
    }

    /*----------------------------------------------------------------------*
      the B-operator used is equivalent to the one used in hex8, this needs
      to be checked if it is ok, but from the mathematics point of view, the only
      thing that needed to be changed is the NUMDOF
      ----------------------------------------------------------------------*/
    /*
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
    // size is 6x12
    LINALG::Matrix<NUMSTR_SOTET4,NUMDOF_SOTET4> bop;
    for (int i=0; i<NUMNOD_SOTET4; i++)
    {
      bop(0,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,0);
      bop(1,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,1);
      bop(2,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,1) + defgrd(0,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,1) + defgrd(1,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,1) + defgrd(2,1)*nxyz(i,0);
      bop(4,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,2) + defgrd(0,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,2) + defgrd(1,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,2) + defgrd(2,2)*nxyz(i,1);
      bop(5,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,0) + defgrd(0,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,0) + defgrd(1,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,0) + defgrd(2,0)*nxyz(i,2);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    LINALG::Matrix<NUMSTR_SOTET4,NUMSTR_SOTET4> cmat(true);
    LINALG::Matrix<NUMSTR_SOTET4,1> stress(true);
    so_tet4_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("no stress data available");
      for (int i = 0; i < NUMSTR_SOTET4; ++i)
        (*elestress)(gp,i) = stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("no stress data available");
      double detF = defgrd.Determinant();

      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> pkstress;
      pkstress(0,0) = stress(0);
      pkstress(0,1) = stress(3);
      pkstress(0,2) = stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = stress(1);
      pkstress(1,2) = stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = stress(2);

      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> temp;
      LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> cauchystress;
      temp.Multiply(1.0/detF,defgrd,pkstress,0.);
      cauchystress.MultiplyNT(temp,defgrd);

      (*elestress)(gp,0) = cauchystress(0,0);
      (*elestress)(gp,1) = cauchystress(1,1);
      (*elestress)(gp,2) = cauchystress(2,2);
      (*elestress)(gp,3) = cauchystress(0,1);
      (*elestress)(gp,4) = cauchystress(1,2);
      (*elestress)(gp,5) = cauchystress(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
    }

    double detJ_w = detJ * (gpweights)[gp];

    // update of internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w,bop,stress,1.0);
    }

    // update of stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      // size is 6x12
      LINALG::Matrix<NUMSTR_SOTET4,NUMDOF_SOTET4> cb;
      cb.Multiply(cmat,bop);          // temporary C . B
      // size is 12x12
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu
      // auxiliary integrated stress
      LINALG::Matrix<NUMSTR_SOTET4,1> sfac(stress);
      // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      sfac.Scale(detJ_w);
      // intermediate Sm.B_L
      double SmB_L[NUMDIM_SOTET4];
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)
      // with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod)
      {
        SmB_L[0] = sfac(0) * nxyz(inod,0) + sfac(3) * nxyz(inod,1) + sfac(5) * nxyz(inod,2);
        SmB_L[1] = sfac(3) * nxyz(inod,0) + sfac(1) * nxyz(inod,1) + sfac(4) * nxyz(inod,2);
        SmB_L[2] = sfac(5) * nxyz(inod,0) + sfac(4) * nxyz(inod,1) + sfac(2) * nxyz(inod,2);
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod)
        {
          double bopstrbop = 0.0;            // intermediate value
          for (int idim=0; idim<NUMDIM_SOTET4; ++idim)
            bopstrbop += nxyz(jnod,idim)* SmB_L[idim];
          (*stiffmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += bopstrbop;
        }
      }
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/


  // static integrator created in any case to safe "if-case"
  const static vector<LINALG::Matrix<NUMNOD_SOTET4,1> > shapefcts4gp = so_tet4_4gp_shapefcts();
  const static vector<double> gpweights4gp = so_tet4_4gp_weights();
  // evaluate mass matrix
  if (massmatrix != NULL)
  {
    //consistent mass matrix evaluated using a 4-point rule
    for (int gp=0; gp<4; gp++)
    {
      double factor = density * detJ * gpweights4gp[gp];
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod)
      {
        ifactor = (shapefcts4gp[gp])(inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod)
        {
          massfactor = (shapefcts4gp[gp])(jnod) * ifactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += massfactor;
        }
      }
    }
  }// end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  return;
} // DRT::ELEMENTS::So_tet4::so_tet4_nlnstiffmass

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_lumpmass(LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4Type::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_tet4* actele = dynamic_cast<DRT::ELEMENTS::So_tet4*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 1 Gauss Point                           |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOTET4,1> > DRT::ELEMENTS::So_tet4::so_tet4_1gp_shapefcts()
{
  vector<LINALG::Matrix<NUMNOD_SOTET4,1> > shapefcts(NUMGPT_SOTET4);

  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int gp=0; gp<NUMGPT_SOTET4; gp++) {
    (shapefcts[gp])(0) = 0.25;
    (shapefcts[gp])(1) = 0.25;
    (shapefcts[gp])(2) = 0.25;
    (shapefcts[gp])(3) = 0.25;
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fct derivs at 1 Gauss Point                     |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > DRT::ELEMENTS::So_tet4::so_tet4_1gp_derivs()
{
  vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > derivs(NUMGPT_SOTET4);
  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int gp=0; gp<NUMGPT_SOTET4; gp++) {
    (derivs[gp])(0,0) = 1;
    (derivs[gp])(1,0) = 0;
    (derivs[gp])(2,0) = 0;
    (derivs[gp])(3,0) = 0;

    (derivs[gp])(0,1) = 0;
    (derivs[gp])(1,1) = 1;
    (derivs[gp])(2,1) = 0;
    (derivs[gp])(3,1) = 0;

    (derivs[gp])(0,2) = 0;
    (derivs[gp])(1,2) = 0;
    (derivs[gp])(2,2) = 1;
    (derivs[gp])(3,2) = 0;

    (derivs[gp])(0,3) = 0;
    (derivs[gp])(1,3) = 0;
    (derivs[gp])(2,3) = 0;
    (derivs[gp])(3,3) = 1;
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Weights at 1 Gauss Point                              |
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_tet4::so_tet4_1gp_weights()
{
  vector<double> weights(NUMGPT_SOTET4);
  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int i = 0; i < NUMGPT_SOTET4; ++i)
    weights[i] = 1.0;
  return weights;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 4 Gauss Points                          |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOTET4,1> > DRT::ELEMENTS::So_tet4::so_tet4_4gp_shapefcts()
{
  vector<LINALG::Matrix<NUMNOD_SOTET4,1> > shapefcts(4);

  const double gploc_alpha    = (5.0 + 3.0*sqrt(5.0))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 - sqrt(5.0))/20.0;

  const double xsi1[4] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
  const double xsi2[4] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
  const double xsi3[4] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
  const double xsi4[4] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};

  for (int gp=0; gp<4; gp++) {
    (shapefcts[gp])(0) = xsi1[gp];
    (shapefcts[gp])(1) = xsi2[gp];
    (shapefcts[gp])(2) = xsi3[gp];
    (shapefcts[gp])(3) = xsi4[gp];
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fct derivs at 4 Gauss Points                    |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > DRT::ELEMENTS::So_tet4::so_tet4_4gp_derivs()
{
  vector<LINALG::Matrix<NUMDIM_SOTET4+1,NUMNOD_SOTET4> > derivs(4);

  for (int gp=0; gp<4; gp++) {
    (derivs[gp])(0,0) = 1;
    (derivs[gp])(1,0) = 0;
    (derivs[gp])(2,0) = 0;
    (derivs[gp])(3,0) = 0;

    (derivs[gp])(0,1) = 0;
    (derivs[gp])(1,1) = 1;
    (derivs[gp])(2,1) = 0;
    (derivs[gp])(3,1) = 0;

    (derivs[gp])(0,2) = 0;
    (derivs[gp])(1,2) = 0;
    (derivs[gp])(2,2) = 1;
    (derivs[gp])(3,2) = 0;

    (derivs[gp])(0,3) = 0;
    (derivs[gp])(1,3) = 0;
    (derivs[gp])(2,3) = 0;
    (derivs[gp])(3,3) = 1;
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Weights at 4 Gauss Points                             |
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_tet4::so_tet4_4gp_weights()
{
  vector<double> weights(4);
  for (int i = 0; i < 4; ++i) {
    weights[i] = 0.25;
  }
  return weights;
}


/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::DefGradient(const vector<double>& disp,
                                         Epetra_SerialDenseMatrix& gpdefgrd,
                                         DRT::ELEMENTS::PreStress& prestress)
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xdisp;
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    // get derivatives wrt to last spatial configuration
    LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> N_xyz;
    prestress_->StoragetoMatrix(gp,N_xyz,prestress_->JHistory());

    // build multiplicative incremental defgrd
    LINALG::Matrix<3,3> defgrd;
    defgrd.MultiplyTN(xdisp,N_xyz);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;

    prestress.MatrixtoStorage(gp,defgrd,gpdefgrd);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::UpdateJacobianMapping(
                                            const vector<double>& disp,
                                            DRT::ELEMENTS::PreStress& prestress)
{
  // get incremental disp
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xdisp;
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> nxyzhist;
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> nxyznew;
  LINALG::Matrix<3,3>                         defgrd;

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    // get the nxyz old state
    prestress.StoragetoMatrix(gp,nxyzhist,prestress.JHistory());
    // build multiplicative incremental defgrd
    defgrd.MultiplyTN(xdisp,nxyzhist);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();

    // push-forward of nxyz
    nxyznew.Multiply(nxyzhist,defgrd);
    // store new reference configuration
    prestress.MatrixtoStorage(gp,nxyznew,prestress.JHistory());

  } // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  return;
}

/*----------------------------------------------------------------------*
 |  remodeling of fiber directions (protected)               tinkl 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_remodel(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      ParameterList&            params,         // algorithmic parameters e.g. time
      RCP<MAT::Material>        mat)            // material
{
// in a first step ommit everything with prestress

  // current  displacements of element
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xdisp;
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  {
    const LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4>& nxyz = nxyz_;

    //                                      d xcurr
    // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    // size is 3x3
    LINALG::Matrix<3,3> defgrd(false);

    if (pstype_==INPAR::STR::prestress_mulf)
    {
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> N_xyz;
      prestress_->StoragetoMatrix(gp,N_xyz,prestress_->JHistory());

      // build multiplicative incremental defgrd
      //defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
      defgrd.MultiplyTN(xdisp,N_xyz);
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
    else
    {
      defgrd.MultiplyTN(xdisp,nxyz);
      defgrd(0,0)+=1;
      defgrd(1,1)+=1;
      defgrd(2,2)+=1;
    }

    // Right Cauchy-Green tensor = F^T * F
    // size is 3x3
    LINALG::Matrix<NUMDIM_SOTET4,NUMDIM_SOTET4> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::Matrix<6,1> glstrain(false);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);


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
    LINALG::Matrix<NUMSTR_SOTET4,NUMDOF_SOTET4> bop;
    for (int i=0; i<NUMNOD_SOTET4; i++)
    {
      bop(0,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,0);
      bop(1,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,1);
      bop(2,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,1) + defgrd(0,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,1) + defgrd(1,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,1) + defgrd(2,1)*nxyz(i,0);
      bop(4,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,2) + defgrd(0,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,2) + defgrd(1,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,2) + defgrd(2,2)*nxyz(i,1);
      bop(5,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,0) + defgrd(0,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,0) + defgrd(1,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,0) + defgrd(2,0)*nxyz(i,2);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOTET4,NUMSTR_SOTET4> cmat(true);
    LINALG::Matrix<NUMSTR_SOTET4,1> stress(true);
    so_tet4_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // Cauchy stress
    const double detF = defgrd.Determinant();

    LINALG::Matrix<3,3> pkstress;
    pkstress(0,0) = stress(0);
    pkstress(0,1) = stress(3);
    pkstress(0,2) = stress(5);
    pkstress(1,0) = pkstress(0,1);
    pkstress(1,1) = stress(1);
    pkstress(1,2) = stress(4);
    pkstress(2,0) = pkstress(0,2);
    pkstress(2,1) = pkstress(1,2);
    pkstress(2,2) = stress(2);

    LINALG::Matrix<3,3> temp(true);
    LINALG::Matrix<3,3> cauchystress(true);
    temp.Multiply(1.0/detF,defgrd,pkstress);
    cauchystress.MultiplyNT(temp,defgrd);

    // evaluate eigenproblem based on stress of previous step
    LINALG::Matrix<3,3> lambda(true);
    LINALG::Matrix<3,3> locsys(true);
    LINALG::SYEV(cauchystress,lambda,locsys);

    // modulation function acc. Hariton: tan g = 2nd max lambda / max lambda
    double newgamma = atan(lambda(1,1)/lambda(2,2));
    //compression in 2nd max direction, thus fibers are alligned to max principal direction
    if (lambda(1,1) < 0) newgamma = 0.0;

    if (mat->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular) {
      MAT::HolzapfelCardio* holz = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holz->EvaluateFiberVecs(gp,newgamma,locsys,defgrd);
    } else if (mat->MaterialType() == INPAR::MAT::m_humphreycardiovascular) {
      MAT::HumphreyCardio* hum = static_cast <MAT::HumphreyCardio*>(mat.get());
      hum->EvaluateFiberVecs(gp,locsys,defgrd);
    } else if (mat->MaterialType() == INPAR::MAT::m_constraintmixture) {
      MAT::ConstraintMixture* comi = static_cast <MAT::ConstraintMixture*>(mat.get());
      comi->EvaluateFiberVecs(gp,locsys,defgrd);
    } else dserror("material not implemented for remodeling");

  } // end loop over gauss points
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
