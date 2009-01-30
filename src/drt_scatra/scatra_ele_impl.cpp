/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl.cpp

\brief Internal implementation of scalar transport elements

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_ele_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/linalg_serialdensematrix.H"

//#define VISUALIZE_ELEMENT_DATA
#include "scatra_element.H" // only for visualization of element data

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraImplInterface* DRT::ELEMENTS::ScaTraImplInterface::Impl(DRT::Element* ele)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  int numscal = numdofpernode;
  if (DRT::Problem::Instance()->ProblemType() == "elch")
    numscal -= 1;

  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    static ScaTraImpl<DRT::Element::hex8>* ch8;
    if (ch8==NULL)
      ch8 = new ScaTraImpl<DRT::Element::hex8>(numdofpernode,numscal);
    return ch8;
  }
/*  case DRT::Element::hex20:
  {
    static ScaTraImpl<DRT::Element::hex20>* ch20;
    if (ch20==NULL)
      ch20 = new ScaTraImpl<DRT::Element::hex20>(numdofpernode,numscal);
    return ch20;
  }
  case DRT::Element::hex27:
  {
    static ScaTraImpl<DRT::Element::hex27>* ch27;
    if (ch27==NULL)
      ch27 = new ScaTraImpl<DRT::Element::hex27>(numdofpernode,numscal);
    return ch27;
  }*/
  case DRT::Element::tet4:
  {
    static ScaTraImpl<DRT::Element::tet4>* ct4;
    if (ct4==NULL)
      ct4 = new ScaTraImpl<DRT::Element::tet4>(numdofpernode,numscal);
    return ct4;
  }
 /* case DRT::Element::tet10:
  {
    static ScaTraImpl<DRT::Element::tet10>* ct10;
    if (ct10==NULL)
      ct10 = new ScaTraImpl<DRT::Element::tet10>(numdofpernode,numscal);
    return ct10;
  } */
  case DRT::Element::wedge6:
  {
    static ScaTraImpl<DRT::Element::wedge6>* cw6;
    if (cw6==NULL)
      cw6 = new ScaTraImpl<DRT::Element::wedge6>(numdofpernode,numscal);
    return cw6;
  }
/*  case DRT::Element::wedge15:
  {
    static ScaTraImpl<DRT::Element::wedge15>* cw15;
    if (cw15==NULL)
      cw15 = new ScaTraImpl<DRT::Element::wedge15>(numdofpernode,numscal);
    return cw15;
  } */
  case DRT::Element::pyramid5:
  {
    static ScaTraImpl<DRT::Element::pyramid5>* cp5;
    if (cp5==NULL)
      cp5 = new ScaTraImpl<DRT::Element::pyramid5>(numdofpernode,numscal);
    return cp5;
  }
  case DRT::Element::quad4:
  {
    static ScaTraImpl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new ScaTraImpl<DRT::Element::quad4>(numdofpernode,numscal);
    return cp4;
  }
/*  case DRT::Element::quad8:
  {
    static ScaTraImpl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new ScaTraImpl<DRT::Element::quad8>(numdofpernode,numscal);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static ScaTraImpl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new ScaTraImpl<DRT::Element::quad9>(numdofpernode,numscal);
    return cp9;
  }*/
  case DRT::Element::tri3:
  {
    static ScaTraImpl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new ScaTraImpl<DRT::Element::tri3>(numdofpernode,numscal);
    return cp3;
  }
/*  case DRT::Element::tri6:
  {
    static ScaTraImpl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new ScaTraImpl<DRT::Element::tri6>(numdofpernode,numscal);
    return cp6;
  }*/
  case DRT::Element::line2:
  {
    static ScaTraImpl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new ScaTraImpl<DRT::Element::line2>(numdofpernode,numscal);
    return cl2;
  }/*
  case DRT::Element::line3:
  {
    static ScaTraImpl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new ScaTraImpl<DRT::Element::line3>(numdofpernode,numscal);
    return cl3;
  }*/
  default:
    dserror("shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraImpl<distype>::ScaTraImpl(int numdofpernode, int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    xyze_(true),  // initialize to zero
    bodyforce_(numdofpernode_),
    diffus_(numscal_),
    valence_(numscal_),
    diffusvalence_(numscal_),
    shcacp_(0.0),
    xsi_(true),
    funct_(true),
    densfunct_(true),
    deriv_(true),
    deriv2_(true),
    xjm_(true),
    xij_(true),
    derxy_(true),
    derxy2_(true),
    rhs_(numdofpernode_),
    hist_(numdofpernode_),
    velint_(true),
    migvelint_(true),
    mderxy_(true),
    mdiv_(0.0),
    tau_(numscal_),
    kart_(numscal_),
    xder2_(true),
    fac_(0.0),
    conv_(true),
    diff_(true),
    migconv_(true),
    migrea_(true),
    gradpot_(true),
    conint_(numscal_),
    gradphi_(true),
    lapphi_(true),
    laplace_(true)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraImpl<distype>::Evaluate(
    DRT::Element*              ele,
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1_epetra,
    Epetra_SerialDenseMatrix&  elemat2_epetra,
    Epetra_SerialDenseVector&  elevec1_epetra,
    Epetra_SerialDenseVector&  elevec2_epetra,
    Epetra_SerialDenseVector&  elevec3_epetra
    )
{
  // get the material
  RefCountPtr<MAT::Material> mat = ele->Material();
  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_condif)
    actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();
  else if (mat->MaterialType()== m_matlist)
    actmat = static_cast<MAT::MatList*>(mat.get())->MaterialData();
  else
    dserror("condif or matlist material expected but got type %d", mat->MaterialType());

  // check for the action parameter
  const string action = params.get<string>("action","none");
  if (action=="calc_condif_systemmat_and_residual")
  {
    // need current history vector and density vector
    RefCountPtr<const Epetra_Vector> hist = discretization.GetState("hist");
    RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (hist==null || densnp==null || phinp==null)
      dserror("Cannot get state vector 'hist', 'densnp' and/or 'phinp'");

    // extract local values from the global vector
    vector<double> myhist(lm.size());
    vector<double> mydensnp(lm.size());
    vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
    DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // get control parameter
    const bool is_stationary = params.get<bool>("using stationary formulation");
    const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");
    const double time = params.get<double>("total time");
    const bool islinear = params.get<bool>("is linear problem");

    // get time-step length
    const double dt = params.get<double>("time-step length");

    // One-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // set parameters for stabilization
    ParameterList& stablist = params.sublist("STABILIZATION");

    // select tau definition
    SCATRA::TauType whichtau = SCATRA::tau_franca_valentin; //default
    {
      const string taudef = stablist.get<string>("DEFINITION_TAU");

      if(taudef == "Zero") whichtau = SCATRA::tau_zero;
      else if(taudef == "Bazilevs")   whichtau = SCATRA::tau_bazilevs;
    }

    // get (weighted) velocity at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    Epetra_SerialDenseVector evel(nsd_*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_);

    // get flag for fine-scale subgrid diffusivity
    string fssgd = params.get<string>("fs subgrid diffusivity","No");

    // check for non-existing subgrid-diffusivity models
    if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" || fssgd == "Smagorinsky_small")
      dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

    // set flag for type of scalar whether it is temperature or not
    string scaltypestr=params.get<string>("problem type");
    bool temperature = false;
    double thermpressdt = 0.0;
    if(scaltypestr =="loma")
    {
      temperature = true;
      thermpressdt = params.get<double>("time derivative of thermodynamic pressure");
    }

    // set flag for conservative form
    string convform = params.get<string>("form of convective term");
    bool conservative = false;
    if(convform =="conservative") conservative = true;

    // get parameter F/RT needed for ELCH ;-)
    double frt(0.0);
    if(scaltypestr =="elch") frt = params.get<double>("frt");

    // create objects for element arrays
    vector<LINALG::Matrix<iel,1> > ephinp(numscal_);
    vector<LINALG::Matrix<iel,1> > ehist(numdofpernode_);
    LINALG::Matrix<iel,1>    edensnp;
    LINALG::Matrix<nsd_,iel> evelnp;
    LINALG::Matrix<iel,1>    epotnp;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      }
      for (int k = 0; k< numdofpernode_; ++k)
      {
        // the history vectors contains information of time step t_n
        ehist[k](i,0) = myhist[k+(i*numdofpernode_)];
      }

      // insert velocity field into element array
      for (int idim=0 ; idim < nsd_ ; idim++)
      {
        evelnp(idim,i) = evel[idim + (i*nsd_) ];
      }

      // insert density vector into element array
      // (only take values belonging to the first transported scalar!)
      edensnp(i,0) = mydensnp[0+(i*numdofpernode_)];

      if(scaltypestr =="elch")
      {
        // get values for el. potential at element nodes
        epotnp(i) = myphinp[i*numdofpernode_+numscal_];
      }
      else
        epotnp(i) = 0.0;
    } // for i

    // calculate element coefficient matrix and rhs
    Sysmat(
        ele,
        ephinp,
        ehist,
        edensnp,
        epotnp,
        elemat1_epetra,
        elevec1_epetra,
        elevec2_epetra,
        actmat,
        time,
        dt,
        timefac,
        alphaF,
        evelnp,
        temperature,
        conservative,
        whichtau,
        fssgd,
        is_stationary,
        is_genalpha,
        islinear,
        frt,
        thermpressdt);
  }
  else if (action =="calc_initial_time_deriv")
    // calculate time derivative for time value t_0
  {
    const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");

    const double time = params.get<double>("total time");
    const double dt = params.get<double>("time-step length");

    // One-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    timefac = params.get<double>("time factor");
    if (is_genalpha)
    {
      alphaF = params.get<double>("alpha_F");
      timefac *= alphaF;
    }
    if (timefac < 0.0) dserror("time factor is negative.");

    // set parameters for stabilization
    ParameterList& stablist = params.sublist("STABILIZATION");

    // select tau definition
    SCATRA::TauType whichtau = SCATRA::tau_franca_valentin; //default
    {
      const string taudef = stablist.get<string>("DEFINITION_TAU");

      if(taudef == "Zero") whichtau = SCATRA::tau_zero;
      else if(taudef == "Bazilevs")   whichtau = SCATRA::tau_bazilevs;
    }

    // need initial field
    RefCountPtr<const Epetra_Vector> phi0 = discretization.GetState("phi0");
    RefCountPtr<const Epetra_Vector> dens0 = discretization.GetState("dens0");
    if (phi0==null || dens0==null)
      dserror("Cannot get state vector 'phi0' and/or 'densnp'");

    // extract local values from the global vector
    vector<double> myphi0(lm.size());
    vector<double> mydens0(lm.size());
    DRT::UTILS::ExtractMyValues(*phi0,myphi0,lm);
    DRT::UTILS::ExtractMyValues(*dens0,mydens0,lm);

    // get initial velocity values at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    Epetra_SerialDenseVector evel(nsd_*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_);

    // get flag for fine-scale subgrid diffusivity
    string fssgd = params.get<string>("fs subgrid diffusivity","No");

    // check for non-existing subgrid-diffusivity models
    if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" || fssgd == "Smagorinsky_small")
      dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

    // set flag for type of scalar whether it is temperature or not
    string scaltypestr=params.get<string>("problem type");
    bool temperature = false;
    double thermpressdt = 0.0;
    if(scaltypestr =="loma")
    {
      temperature = true;
      thermpressdt = params.get<double>("time derivative of thermodynamic pressure");
    }

    // set flag for conservative form
    string convform = params.get<string>("form of convective term");
    bool conservative = false;
    if(convform =="conservative") conservative = true;

    // create objects for element arrays
    vector<LINALG::Matrix<iel,1> > ephi0(numscal_);
    LINALG::Matrix<iel,1> edens0;
    LINALG::Matrix<nsd_,iel> evel0;
    LINALG::Matrix<iel,1> epot0;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephi0[k](i,0) = myphi0[k+(i*numdofpernode_)];
      }

      // split for each tranported scalar, insert into element arrays
      for (int idim = 0; idim< nsd_; idim++)
      {
        evel0(idim,i) = evel[idim + (i*nsd_)];
      }

      // insert density vector into element array
      // (only take values belonging to the first transported scalar!)
      edens0(i,0) = mydens0[0+(i*numdofpernode_)];

      if(scaltypestr =="elch")
      {
        // get values for el. potential at element nodes
        epot0(i) = myphi0[i*numdofpernode_+numscal_];
      }
      else
        epot0(i) = 0.0;
    } // for i

    double frt(0.0);
    if(scaltypestr =="elch")
    {
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // calculate mass matrix and rhs
    InitialTimeDerivative(
        ele,
        ephi0,
        edens0,
        epot0,
        elemat1_epetra,
        elevec1_epetra,
        elevec2_epetra,
        actmat,
        time,
        dt,
        timefac,
        evel0,
        temperature,
        conservative,
        whichtau,
        fssgd,
        frt,
        thermpressdt);
  }
  else if (action =="calc_subgrid_diffusivity_matrix")
  // calculate normalized subgrid-diffusivity matrix
  {
    // get control parameter
    const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");
    const bool is_stationary = params.get<bool>("using stationary formulation");

    // One-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // calculate mass matrix and rhs
    CalcSubgridDiffMatrix(
        ele,
        elemat1_epetra,
        timefac,
        is_stationary);
  }
  else if (action=="calc_condif_flux")
  {
    // get velocity values at the nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    Epetra_SerialDenseVector evel(nsd_*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_);

    // need current values of transported scalar
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vectors
    vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // assure, that the values are in the same order as the element nodes
    for(int k=0;k<iel;++k)
    {
      Node* node = (ele->Nodes())[k];
      vector<int> dof = discretization.Dof(node);
      int numdof = dof.size();
      for (int i=0;i<numdof;++i)
      {
        if (dof[i]!=lm[k*numdof+i])
        {
          cout<<"dof[i]= "<<dof[i]<<"  lm[k*numdof+i]="<<lm[k*numdof+i]<<endl;
          dserror("Dofs are not in the same order as the element nodes. Implement some resorting!");
        }
      }
    }

    // access control parameter
    SCATRA::FluxType fluxtype;
    string fluxtypestring = params.get<string>("fluxtype","noflux");
    if (fluxtypestring == "totalflux")
      fluxtype = SCATRA::totalflux;
    else if (fluxtypestring == "diffusiveflux")
      fluxtype = SCATRA::diffusiveflux;
    else
      fluxtype=SCATRA::noflux;  //default value

    // set flag for type of scalar
    string scaltypestr=params.get<string>("problem type");
    int numscal = numdofpernode_;
    bool temperature = false;
    if (scaltypestr =="loma") temperature = true;

    double frt(0.0);
    if (scaltypestr =="elch")
    {
      numscal -= 1; // ELCH case: last dof is for el. potential
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // we always get an 3D flux vector for each node
    LINALG::Matrix<3,iel> eflux(true);

    // do a loop for systems of transported scalars
    for (int i = 0; i<numscal; ++i)
    {
      // calculate flux vectors for actual scalar
      eflux.Clear();
      CalculateFlux(eflux,ele,myphinp,actmat,temperature,frt,evel,fluxtype,i);

      // assembly
      for (int k=0;k<iel;k++)
      { // form arithmetic mean of assembled nodal flux vectors
        // => factor is the number of adjacent elements for each node
        double factor = ((ele->Nodes())[k])->NumElement();
        elevec1_epetra[k*numdofpernode_+i]+=eflux(0,k)/factor;
        elevec2_epetra[k*numdofpernode_+i]+=eflux(1,k)/factor;
        elevec3_epetra[k*numdofpernode_+i]+=eflux(2,k)/factor;
      }
    } // loop over numscal
  }
  else if (action=="calc_temp_and_dens")
  {
    // calculate mean temperature and density

    // NOTE: add integral values only for elements which are NOT ghosted!
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      // need current scalar and density vector
      RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
      RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
      if (phinp==null || densnp==null)
        dserror("Cannot get state vector 'phinp' and/or 'densnp'");

      // extract local values from the global vectors
      vector<double> myphinp(lm.size());
      vector<double> mydensnp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
      DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);

      // calculate temperature/concentrations, density and domain integral
      CalculateTempAndDens(ele,myphinp,mydensnp,elevec1_epetra);
    }
  }
  else if (action=="calc_domain_and_bodyforce")
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      const double time = params.get<double>("total time");

      // calculate domain and bodyforce integral
      CalculateDomainAndBodyforce(elevec1_epetra,ele,time);
    }
  }
  else if (action=="calc_elch_kwok_error")
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    // need current solution
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // create objects for element arrays
    LINALG::Matrix<iel,2> ephinp;
    LINALG::Matrix<iel,1> epotnp;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< 2; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephinp(i,k) = myphinp[k+(i*numdofpernode_)];
      }

      // get values for el. potential at element nodes
      epotnp(i) = myphinp[i*numdofpernode_+numscal_];
    } // for i

    CalErrorComparedToAnalytSolution(
        ele,
        params,
        ephinp,
        epotnp,
        elevec1_epetra,
        actmat);
  }
  else
    dserror("Unknown type of action for Scatra Implementation: %s",action.c_str());

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate mass flux                              (private) gjb 06/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateFluxSerialDense(
    LINALG::SerialDenseMatrix&      flux,
    DRT::Element*&            ele,
    vector<double>&           ephinp,
    struct _MATERIAL*         material,
    bool&                     temperature,
    double&                   frt,
    Epetra_SerialDenseVector& evel,
    std::string&              fluxtypestring,
    int&                      dofindex
)
{
  // access control parameter
  SCATRA::FluxType fluxtype;
  if (fluxtypestring == "totalflux")
    fluxtype = SCATRA::totalflux;
  else if (fluxtypestring == "diffusiveflux")
    fluxtype = SCATRA::diffusiveflux;
  else
    fluxtype=SCATRA::noflux;  //default value

  // we always get an 3D flux vector for each node
  LINALG::Matrix<3,iel> eflux(true);
  eflux.Clear();
  CalculateFlux(eflux,ele,ephinp,material,temperature,frt,evel,fluxtype,dofindex);
  for (int j = 0; j< iel; j++)
  {
    flux(0,j) = eflux(0,j);
    flux(1,j) = eflux(1,j);
    flux(2,j) = eflux(2,j);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (public)                 g.bau 08/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::Sysmat(
    DRT::Element*                         ele, ///< the element those matrix is calculated
    const vector<LINALG::Matrix<iel,1> >& ephinp,///< current scalar field
    const vector<LINALG::Matrix<iel,1> >& ehist, ///< rhs from beginning of time step
    const LINALG::Matrix<iel,1>&          edensnp, ///< current density field
    const LINALG::Matrix<iel,1>&          epotnp, ///< el. potential at element nodes
    Epetra_SerialDenseMatrix&             sys_mat,///< element matrix to calculate
    Epetra_SerialDenseVector&             residual, ///< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff, ///< subgrid-diff.-scaling vector
    struct _MATERIAL*                     material, ///< material pointer
    const double                          time, ///< current simulation time
    const double                          dt, ///< current time-step length
    const double                          timefac, ///< time discretization factor
    const double                          alphaF, ///< factor for generalized-alpha time integration
    const LINALG::Matrix<nsd_,iel>&       evelnp,///< nodal velocities at t_{n+1}
    const bool                            temperature, ///< temperature flag
    const bool                            conservative, ///< flag for conservative form
    const enum SCATRA::TauType            whichtau, ///< flag for stabilization parameter definition
    const string                          fssgd, ///< subgrid-diff. flag
    const bool                            is_stationary, ///< stationary flag
    const bool                            is_genalpha, ///< generalized-alpha flag
    const bool                            islinear, ///< flag for linear/nonlinear problem
    const double                          frt, ///< factor F/RT needed for ELCH calculations
    const double                          thermpressdt ///< time deriv. of therm. press. (for temp. eq.)
)
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele,time);

  // get material constants
  GetMaterialParams(material,temperature);

  /*----------------------------------------------------------------------*/
  // calculation of stabilization parameter(s) tau
  /*----------------------------------------------------------------------*/
  CalTau(ele,subgrdiff,evelnp,edensnp,epotnp,dt,timefac,whichtau,fssgd,is_stationary,false,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool use2ndderiv = SCATRA::useSecondDerivatives<distype>();

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,use2ndderiv,ele->Id());

    // density-weighted shape functions
    densfunct_.EMultiply(funct_,edensnp);

    // get (density-weighted) velocity at integration point
    velint_.Multiply(evelnp,funct_);

    if (conservative)
    {
      // get momentum (density-weighted velocity) derivatives at integration point
      mderxy_.MultiplyNT(evelnp,derxy_);

      // momentum divergence
      GetMomentumDivergence(mdiv_,mderxy_);
    }

    //------------ get values of variables at integration point
    for (int k = 0;k<numdofpernode_;++k)     // loop of each transported sclar
    {
      // get history data at integration point
      if (is_genalpha and not conservative)
        hist_[k] = densfunct_.Dot(ehist[k]);
      else
        hist_[k] = funct_.Dot(ehist[k]);

      // get bodyforce in gausspoint (divided by shcacp)
      // (For temperature equation, time derivative of thermodynamic pressure
      //  is added, if not constant.)
      rhs_[k] = bodyforce_[k].Dot(funct_) / shcacp_;
      rhs_[k] += thermpressdt / shcacp_;
    }

    //----------- perform integration for entire matrix and rhs
    if (numdofpernode_-numscal_== 0) // 'standard' scalar transport
    {
      if (islinear)
      {
        for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
        {
          if (not is_stationary)
            CalMat(sys_mat,residual,ephinp,use2ndderiv,conservative,is_genalpha,timefac,alphaF,k);
          else
            CalMatStationary(sys_mat,residual,use2ndderiv,conservative,k);
        } // loop over each scalar
      }
      else
        CalMatInc(sys_mat,residual,ephinp,use2ndderiv,is_stationary,timefac);
    }
    else  // ELCH problems
    {
     if(is_genalpha) dserror("GenAlpha is not supported by ELCH!");
     CalMatElch(sys_mat,residual,ephinp,epotnp,use2ndderiv,frt,is_stationary,timefac);
    }

  } // integration loop

  return;
}


/*----------------------------------------------------------------------*
 |  get the body force  (private)                              gjb 06/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::BodyForce(
    const DRT::Element*    ele,
    const double           time
)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  switch(nsd_)
  {
  case 3:
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
  break;
  case 2:
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
  break;
  case 1:
    DRT::UTILS::FindElementConditions(ele, "LineNeumann", myneumcond);
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
  }

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac(0.0);

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
        // A negative time value indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

    // set this condition to the bodyforce array
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      for (int jnode=0; jnode<iel; jnode++)
      {
        (bodyforce_[idof])(jnode) = (*onoff)[idof]*(*val)[idof]*curvefac;
      }
    }
  }
  else
  {
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      // we have no dead load
      bodyforce_[idof].Clear();
    }
  }

  return;

} //ScaTraImpl::BodyForce


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                      gjb 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetMaterialParams(
    struct _MATERIAL*  material,
    const bool&        temperature
)
{
// get diffusivity / diffusivities
if (material->mattyp == m_matlist)
{
  for (int k = 0;k<numscal_;++k)
  {
    const int matid = material->m.matlist->matids[k];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_ion)
    {
      valence_[k]= singlemat.m.ion->valence;
      diffus_[k]= singlemat.m.ion->diffusivity;
      diffusvalence_[k] = valence_[k]*diffus_[k];
    }
    else if (singlemat.mattyp == m_condif)
      diffus_[k]= singlemat.m.condif->diffusivity;
    else
      dserror("material type is not allowed");
  }
  // set specific heat capacity at constant pressure to 1.0
  shcacp_ = 1.0;
}
else if (material->mattyp == m_condif)
{
  dsassert(numdofpernode_==1,"more than 1 dof per node for condif material");

  // in case of a temperature equation, we get thermal conductivity instead of
  // diffusivity and have to divide by the specific heat capacity at constant
  // pressure; otherwise, it is the "usual" diffusivity
  if (temperature)
  {
    shcacp_ = material->m.condif->shc;
    diffus_[0] = material->m.condif->diffusivity/shcacp_;
  }
  else
  {
    // set specific heat capacity at constant pressure to 1.0, get diffusivity
    shcacp_ = 1.0;
    diffus_[0] = material->m.condif->diffusivity;
  }
}
else
  dserror("Material type is not supported");

return;
} //ScaTraImpl::GetMaterialParams


/*----------------------------------------------------------------------*
 |  calculate stabilization parameter  (private)              gjb 06/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalTau(
    DRT::Element*                   ele,
    Epetra_SerialDenseVector&       subgrdiff,
    const LINALG::Matrix<nsd_,iel>& evel,
    const LINALG::Matrix<iel,1>&    edens,
    const LINALG::Matrix<iel,1>&    epot,
    const double                    dt,
    const double&                   timefac,
    const enum SCATRA::TauType      whichtau,
    const string                    fssgd,
    const bool&                     is_stationary,
    const bool                      initial,
    const double&                   frt
  )
{
  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // use one-point Gauss rule to calculate tau at element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // coordinates of the integration point
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim = 0; idim < nsd_; idim++)
    {xsi_(idim) = gpcoord[idim];}

  // integration weight
  const double wquad = intpoints_tau.IP().qwgt[0];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(), det);

  // get (density-weighted) velocity at element center
  velint_.Multiply(evel,funct_);

  // get "migration velocity" divided by D_k*z_k at element center
  if (numdofpernode_-numscal_== 1) // ELCH
  {
    // compute global derivatives
    derxy_.Multiply(xij_,deriv_);

    migvelint_.Multiply(-frt,derxy_,epot);
  } // if ELCH

  // stabilization parameter definition according to Bazilevs et al. (2007)
  switch (whichtau){
  case SCATRA::tau_bazilevs:
  {
    for(int k = 0;k<numscal_;++k) // loop over all transported scalars
    {
      // effective velocity at element center:
      // (weighted) convective velocity + individual migration velocity
      LINALG::Matrix<nsd_,1> veleff(velint_,false);
      if (numdofpernode_ - numscal_ == 1) // ELCH
      {
        const double Dkzk = diffus_[k]*valence_[k];
        veleff.Update(Dkzk,migvelint_,1.0);
      }

      /*
                                                                1.0
               +-                                          -+ - ---
               |                                            |   2.0
               | 4.0    n+1       n+1             2         |
        tau  = | --- + u     * G u     + C * kappa  * G : G |
               |   2           -          I           -   - |
               | dt            -                      -   - |
               +-                                          -+

      */

      /*            +-           -+   +-           -+   +-           -+
                    |             |   |             |   |             |
                    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
              G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
               ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                    |    i     j  |   |    i     j  |   |    i     j  |
                    +-           -+   +-           -+   +-           -+
      */
      /*            +----
                     \
            G : G =   +   G   * G
            -   -    /     ij    ij
            -   -   +----
                     i,j
      */
      /*                      +----
             n+1       n+1     \     n+1          n+1
            u     * G u     =   +   u    * G   * u
                    -          /     i     -ij    j
                    -         +----        -
                               i,j
      */
      double G;
      double normG(0.0);
      double Gnormu(0.0);
      for (int nn=0;nn<nsd_;++nn)
      {
        for (int rr=0;rr<nsd_;++rr)
        {
          G = xij_(nn,0)*xij_(rr,0);
          for(int tt=1;tt<nsd_;tt++)
          {
            G += xij_(nn,tt)*xij_(rr,tt);
          }
          normG+=G*G;
          Gnormu+=veleff(nn,0)*G*veleff(rr,0);
        }
      }

      // definition of constant
      // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
      //  brought 144.0 from Austin...)
      const double CI = 12.0/mk;

      // stabilization parameters for instationary and stationary case, respectively
      if (is_stationary == false)
      {
        // get density at element center
        const double dens = funct_.Dot(edens);

        tau_[k] = 1.0/(sqrt((4.0*dens*dens)/(dt*dt)+Gnormu+CI*diffus_[k]*diffus_[k]*normG));
      }
      else
        tau_[k] = 1.0/(sqrt(Gnormu+CI*diffus_[k]*diffus_[k]*normG));

      // compute artificial diffusivity kappa_art_[k] if required
      if (fssgd == "artificial_all" and (not initial))
      {
        // get Euclidean norm of (weighted) velocity at element center
        const double vel_norm = velint_.Norm2();

        kart_[k] = DSQR(vel_norm)/(sqrt(Gnormu+CI*diffus_[k]*diffus_[k]*normG));

        for (int vi=0; vi<iel; ++vi)
        {
          subgrdiff(vi) = kart_[k]/ele->Nodes()[vi]->NumElement();
        }
      } // for artificial diffusivity
    } // for k
  }
  break;
  case SCATRA::tau_franca_valentin:
  // stabilization parameter definition according to Franca and Valentin (2000)
  {
    // volume of the element (2D: element surface area; 1D: element length)
    // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
    const double vol = wquad*det;

    // There exist different definitions for 'the' characteristic element length hk:
    // 1) get element length for tau_Mp/tau_C: volume-equival. diameter -> not default
    // const double hk = pow((6.*vol/PI),(1.0/3.0));

    // 2) streamlength (based on velocity vector at element centre) -> not default

    // 3) use cubic root of the element volume as characteristic length -> default
    //    2D case: characterisitc length is the square root of the element area
    //    1D case: characteristic length is the element length
    const double dim = (double) nsd_;
    const double hk = pow(vol,(1.0/dim));

    // some necessary parameter definitions
    double vel_norm, epe1, epe2, xi1, xi2;

    for(int k = 0;k<numscal_;++k) // loop over all transported scalars
    {
      if (numdofpernode_ - numscal_ == 1) // ELCH
      {
        const double Dkzk = diffus_[k]*valence_[k];
        // get Euclidean norm of effective velocity at element center:
        // (weighted) convective velocity + individual migration velocity
        LINALG::Matrix<nsd_,1> veleff(velint_,false);
        veleff.Update(Dkzk,migvelint_,1.0);
        vel_norm = veleff.Norm2();

#ifdef VISUALIZE_ELEMENT_DATA
        veleff.Update(Dkzk,migvelint_,0.0);
        double vel_norm_mig = veleff.Norm2();
        double migepe2 = mk * vel_norm_mig * hk / diffus_[k];

        DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
        if (!actele) dserror("cast to Transport* failed");
        vector<double> v(1,migepe2);
        ostringstream temp;
        temp << k;
        string name = "Pe_mig_"+temp.str();
        actele->AddToData(name,v);
#endif
      }
      else
      {
        // get Euclidean norm of (weighted) velocity at element center
        vel_norm = velint_.Norm2();
      }

      // check whether there is zero diffusivity
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");

      // stabilization parameter for instationary case
      if (is_stationary == false)
      {
        /* parameter relating diffusive : reactive forces */
        epe1 = 2.0 * timefac * diffus_[k] / (mk * DSQR(hk));
        /* parameter relating convective : diffusive forces */
        epe2 = mk * vel_norm * hk / diffus_[k];
        xi1 = DMAX(epe1,1.0);
        xi2 = DMAX(epe2,1.0);

        tau_[k] = DSQR(hk)/((DSQR(hk)*xi1)/timefac + (2.0*diffus_[k]/mk)*xi2);
      }
      // stabilization parameter for stationary case
      else
      {
        /* parameter relating convective : diffusive forces */
        epe2 = mk * vel_norm * hk / diffus_[k];
        xi2 = DMAX(epe2,1.0);

        tau_[k] = (DSQR(hk)*mk)/(2.0*diffus_[k]*xi2);
      }

      // compute artificial diffusivity kappa_art_[k]
      if (fssgd == "artificial_all" and (not initial))
      {
        /* parameter relating convective : diffusive forces */
        epe2 = mk * vel_norm * hk / diffus_[k];
        xi2 = DMAX(epe2,1.0);

        kart_[k] = (DSQR(hk)*mk*DSQR(vel_norm))/(2.0*diffus_[k]*xi2);

        for (int vi=0; vi<iel; ++vi)
        {
          subgrdiff(vi) = kart_[k]/ele->Nodes()[vi]->NumElement();
        }
      }

#ifdef VISUALIZE_ELEMENT_DATA 
      // visualize resultant Pe number and stabilization parameter
      DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
      if (!actele) dserror("cast to Transport* failed");
      vector<double> v(1,epe2);
      ostringstream temp;
      temp << k;
      string name = "Pe_"+temp.str();
      actele->AddToData(name,v);
      v[0] = tau_[k];
      name = "tau_"+ temp.str();
      actele->AddToData(name,v);
#endif

    } // loop over scalars
  }
  break;
  case SCATRA::tau_zero:
  {
    // set tau's to zero (-> no stabilization effect)
    for(int k = 0;k<numscal_;++k)
      tau_[k] = 0.0;
  }
  break;
  default: dserror("Unknown definition of tau\n");
  } //switch (whichtau)

  return;
} //ScaTraImpl::Caltau


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     gjb 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int&                                   iquad,      ///< id of current Gauss point
    const bool&                                  use2ndderiv,///< are second derivatives needed?
    const int&                                   eleid       ///< the element id
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
    {xsi_(idim) = gpcoord[idim];}

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // set integration factor: fac = Gauss weight * det(J)
  fac_ = intpoints.IP().qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderiv)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);

    // get global second derivatives
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else
    derxy2_.Clear();

  // say goodbye
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate instationary convection-diffusion matrix (private)gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized scalar transport element are calculated for the 
instationary case. The procedure is based on the Rothe method of first 
discretizing in time. Hence the resulting terms include coefficients 
containing time integration variables such as theta or delta t which are 
represented by 'timefac'.

The stabilization is based on the residuum:

R = rho * c_p * phi + timefac * rho * c_p * u * grad(phi)
                    - timefac * diffus * laplace(phi) - rhsint

The corresponding weighting operators are
L = timefac * rho * c_p * u * grad(w) +/- timefac * diffus * laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMat(
    Epetra_SerialDenseMatrix&             estif,
    Epetra_SerialDenseVector&             eforce,
    const vector<LINALG::Matrix<iel,1> >&  ephinp,
    const bool                            use2ndderiv,
    const bool                            conservative,
    const bool                            is_genalpha,
    const double                          timefac,
    const double                          alphaF,
    const int                             dofindex
    )
{
// number of degrees of freedom per node
const int numdof = numdofpernode_;

// stabilization parameter and integration factors
const double taufac     = tau_[dofindex]*fac_;
const double timefacfac = timefac*fac_;
const double timetaufac = timefac*taufac;
const double fac_diffus = timefacfac*diffus_[dofindex];

// evaluate rhs at integration point
double rhsint = hist_[dofindex] + rhs_[dofindex]*(timefac/alphaF);

// convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
conv_.MultiplyTN(derxy_,velint_);

// diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
if (use2ndderiv)
{
  GetLaplacianStrongForm(diff_, derxy2_);
  diff_.Scale(diffus_[dofindex]);
}

//----------------------------------------------------------------
// element matrix: standard Galerkin terms
//----------------------------------------------------------------
// transient term
for (int vi=0; vi<iel; ++vi)
{
  const double v = fac_*funct_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*densfunct_(ui);
  }
}

// convective term in convective form
for (int vi=0; vi<iel; ++vi)
{
  const double v = timefacfac*funct_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*conv_(ui);
  }
}

// addition to convective term for conservative form
if (conservative)
{
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = timefacfac*funct_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) += v*funct_(ui)*mdiv_;
    }
  }
}

// diffusive term
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;
    double laplawf(0.0);
    GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
    estif(fvi,fui) += fac_diffus*laplawf;
  }
}

//----------------------------------------------------------------
// element matrix: stabilization terms
//----------------------------------------------------------------
// convective stabilization of transient term (in convective form)
for (int vi=0; vi<iel; ++vi)
{
  const double v = taufac*conv_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*densfunct_(ui);
  }
}

// convective stabilization of convective term (in convective form)
for (int vi=0; vi<iel; ++vi)
{
  const double v = timetaufac*conv_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*conv_(ui);
  }
}

if (use2ndderiv)
{
  // convective stabilization of diffusive term (in convective form)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = timetaufac*conv_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) -= v*diff_(ui);
    }
  }

  // diffusive stabilization of transient term
  // (USFEM assumed here, sign change necessary for GLS)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = taufac*diff_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) += v*densfunct_(ui);
    }
  }

  // diffusive stabilization of convective term (in convective form)
  // (USFEM assumed here, sign change necessary for GLS)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = timetaufac*diff_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) += v*conv_(ui);
    }
  }

  // diffusive stabilization of diffusive term
  // (USFEM assumed here, sign change necessary for GLS)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = timetaufac*diff_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) -= v*diff_(ui);
    }
  }
}

//----------------------------------------------------------------
// element right hand side: standard Galerkin bodyforce term
//----------------------------------------------------------------
double vrhs = fac_*rhsint;
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  eforce[fvi] += vrhs*funct_(vi);
}

//----------------------------------------------------------------
// element right hand side: stabilization terms
//----------------------------------------------------------------
// convective stabilization of bodyforce term
vrhs = taufac*rhsint;
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  eforce[fvi] += vrhs*conv_(vi);
}

// diffusive stabilization of bodyforce term (only for higher-order elements)
// (USFEM assumed here, sign change necessary for GLS)
if (use2ndderiv)
{
  for (int vi=0; vi<iel; ++vi)
  {
    const int fvi = vi*numdof+dofindex;

    eforce[fvi] += vrhs*diff_(vi);
  }
}

//----------------------------------------------------------------
// part of element right hand side only required for
// generalized-alpha time integration: temporal terms
//----------------------------------------------------------------
if (is_genalpha)
{
  // integration factors for temporal rhs
  const double rhstimefacfac = timefacfac*(1.0-alphaF)/alphaF;
  const double rhstimetaufac = timetaufac*(1.0-alphaF)/alphaF;

  // gradient of scalar at time step n
  gradphi_.Multiply(derxy_,ephinp[dofindex]);

  // convective part in convective form at time step n
  const double convn = velint_.Dot(gradphi_);

  // convective temporal rhs term in convective form
  vrhs = rhstimefacfac*convn;
  // addition to convective temporal rhs term for conservative form
  if (conservative)
  {
    // scalar at integration point at time step n
    const double phi = funct_.Dot(ephinp[dofindex]);

    // convective temporal rhs term in conservative form
    vrhs += rhstimefacfac*phi*mdiv_;
  }
  for (int vi=0; vi<iel; ++vi)
  {
    const int fvi = vi*numdof+dofindex;

    eforce[fvi] -= vrhs*funct_(vi);
  }

  // diffusive temporal rhs term
  vrhs = rhstimefacfac*diffus_[dofindex];
  for (int vi=0; vi<iel; ++vi)
  {
    const int fvi = vi*numdof+dofindex;

    eforce[fvi] -= vrhs*(derxy_(0,vi)*gradphi_(0)+derxy_(1,vi)*gradphi_(1)+derxy_(2,vi)*gradphi_(2));
  }

  // convective stabilization of convective temporal rhs term (in convective form)
  vrhs = rhstimetaufac*convn;
  for (int vi=0; vi<iel; ++vi)
  {
    const int fvi = vi*numdof+dofindex;

    eforce[fvi] -= vrhs*conv_(vi);
  }

  double diffn = 0.0;
  if (use2ndderiv)
  {
    // second gradient (Laplacian) of scalar at time step n
    lapphi_.Multiply(derxy2_,ephinp[dofindex]);

    // diffusive part at time step n
    diffn = diffus_[dofindex] * (lapphi_(0) + lapphi_(1) + lapphi_(2));

    // diffusive stabilization of convective temporal rhs term (in convective form)
    vrhs = rhstimetaufac*convn;
    for (int vi=0; vi<iel; ++vi)
    {
      const int fvi = vi*numdof+dofindex;

      eforce[fvi] -= vrhs*diff_(vi);
    }

    // convective stabilization of diffusive temporal rhs term
    vrhs = rhstimetaufac*diffn;
    for (int vi=0; vi<iel; ++vi)
    {
      const int fvi = vi*numdof+dofindex;

      eforce[fvi] -= vrhs*conv_(vi);
    }

    // diffusive stabilization of diffusive temporal rhs term
    vrhs = rhstimetaufac*diffn;
    for (int vi=0; vi<iel; ++vi)
    {
      const int fvi = vi*numdof+dofindex;

      eforce[fvi] -= vrhs*diff_(vi);
    }
  }
}

return;
} //ScaTraImpl::CalMat


/*----------------------------------------------------------------------*
 |  evaluate stationary convection-diffusion matrix (private)  gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized scalar transport element are calculated for the 
stationary case.

The stabilization is based on the residuum:

R = rho * c_p * u * grad(phi) - diffus *  laplace(phi) - rhsint

The corresponding weighting operators are
L = rho * c_p * u * grad(w) +/- diffus *  laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatStationary(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    const bool                use2ndderiv,
    const bool                conservative,
    const int                 dofindex
    )
{
// number of degrees of freedom per node
const int numdof = numdofpernode_;

// stabilization parameter and integration factor
const double taufac     = tau_[dofindex]*fac_;
const double fac_diffus = fac_*diffus_[dofindex];

// evaluate rhs at integration point
double rhsint = rhs_[dofindex];

// convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
conv_.MultiplyTN(derxy_,velint_);

// diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
if (use2ndderiv)
{
  GetLaplacianStrongForm(diff_, derxy2_);
  diff_.Scale(diffus_[dofindex]);
}

//----------------------------------------------------------------
// element matrix: standard Galerkin terms
//----------------------------------------------------------------
// convective term in convective form
for (int vi=0; vi<iel; ++vi)
{
  const double v = fac_*funct_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*conv_(ui);
  }
}

// addition to convective term for conservative form
if (conservative)
{
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = fac_*funct_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) += v*funct_(ui)*mdiv_;
    }
  }
}

// diffusive term
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;
    double laplawf(0.0);
    GetLaplacianWeakForm(laplawf,derxy_,vi,ui);
    estif(fvi,fui) += fac_diffus*laplawf;
  }
}

//----------------------------------------------------------------
// element matrix: stabilization terms
//----------------------------------------------------------------
// convective stabilization of convective term (in convective form)
for (int vi=0; vi<iel; ++vi)
{
  const double v = taufac*conv_(vi);
  const int fvi = vi*numdof+dofindex;

  for (int ui=0; ui<iel; ++ui)
  {
    const int fui = ui*numdof+dofindex;

    estif(fvi,fui) += v*conv_(ui);
  }
}

if (use2ndderiv)
{
  // convective stabilization of diffusive term (in convective form)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = taufac*conv_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) -= v*diff_(ui);
    }
  }

  // diffusive stabilization of convective term (in convective form)
  // (USFEM assumed here, sign change necessary for GLS)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = taufac*diff_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) += v*conv_(ui);
    }
  }

  // diffusive stabilization of diffusive term
  // (USFEM assumed here, sign change necessary for GLS)
  for (int vi=0; vi<iel; ++vi)
  {
    const double v = taufac*diff_(vi);
    const int fvi = vi*numdof+dofindex;

    for (int ui=0; ui<iel; ++ui)
    {
      const int fui = ui*numdof+dofindex;

      estif(fvi,fui) -= v*diff_(ui);
    }
  }
}

//----------------------------------------------------------------
// element right hand side: standard Galerkin bodyforce term
//----------------------------------------------------------------
double vrhs = fac_*rhsint;
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  eforce[fvi] += vrhs*funct_(vi);
}

//----------------------------------------------------------------
// element right hand side: stabilization terms
//----------------------------------------------------------------
// convective stabilization of bodyforce term
vrhs = taufac*rhsint;
for (int vi=0; vi<iel; ++vi)
{
  const int fvi = vi*numdof+dofindex;

  eforce[fvi] += vrhs*conv_(vi);
}

// diffusive stabilization of bodyforce term
// (USFEM assumed here, sign change necessary for GLS)
if (use2ndderiv)
{
  for (int vi=0; vi<iel; ++vi)
  {
    const int fvi = vi*numdof+dofindex;

    eforce[fvi] += vrhs*diff_(vi);
  }
}

return;
} //ScaTraImpl::CalMatStationary


/*----------------------------------------------------------------------*
 | calculate mass matrix + rhs for determ. initial time deriv. gjb 08/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::InitialTimeDerivative(
    DRT::Element*                         ele,
    const vector<LINALG::Matrix<iel,1> >& ephi0,
    const LINALG::Matrix<iel,1>&          edens0,
    const LINALG::Matrix<iel,1>&          epot0,
    Epetra_SerialDenseMatrix&             massmat,
    Epetra_SerialDenseVector&             rhs,
    Epetra_SerialDenseVector&             subgrdiff,
    struct _MATERIAL*                     material,
    const double                          time,
    const double                          dt,
    const double                          timefac,
    const LINALG::Matrix<nsd_,iel>&       evel0,
    const bool                            temperature,
    const bool                            conservative,
    const enum SCATRA::TauType            whichtau,
    const string                          fssgd,
    const double                          frt,
    const double                          thermpressdt
)
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // dead load in element nodes
  BodyForce(ele,time);

  // get material constants
  GetMaterialParams(material,temperature);

  /*----------------------------------------------------------------------*/
  // calculation of instationary(!) stabilization parameter(s)
  /*----------------------------------------------------------------------*/
  CalTau(ele,subgrdiff,evel0,edens0,epot0,dt,timefac,whichtau,fssgd,false,true,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool use2ndderiv = SCATRA::useSecondDerivatives<distype>();

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,use2ndderiv,ele->Id());

    // density-weighted shape functions
    densfunct_.EMultiply(funct_,edens0);

    // get (density-weighted) velocity at element center
    velint_.Multiply(evel0,funct_);

    if (conservative)
    {
      // get momentum (density-weighted velocity) derivatives at integration point
      mderxy_.MultiplyNT(evel0,derxy_);

      // momentum divergence
      GetMomentumDivergence(mdiv_,mderxy_);
    }

    //------------ get values of variables at integration point
    for (int k = 0;k<numdofpernode_;++k)     // loop of each transported sclar
    {
      // get bodyforce in gausspoint (divided by shcacp)
      // (For temperature equation, time derivative of thermodynamic pressure
      //  is added, if not constant.)
      rhs_[k] = bodyforce_[k].Dot(funct_) / shcacp_;
      rhs_[k] += thermpressdt / shcacp_;
    }

    // get values of all transported scalars at integration point
    for (int k=0; k<numscal_; ++k)
    {
      conint_[k] = funct_.Dot(ephi0[k]);
    }

    // get gradient of el. potential at integration point
    gradpot_.Multiply(derxy_,epot0);

    // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
    conv_.MultiplyTN(derxy_,velint_);

    // migration part
    migconv_.MultiplyTN(-frt,derxy_,gradpot_);

    /*-------------- perform integration for entire matrix and rhs ---*/
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // stabilization parameter  and integration factor
      const double taufac     = tau_[k]*fac_;
      const double fac_diffus = fac_*diffus_[k];

      // evaluate rhs at integration point
      double rhsint = rhs_[k];

      if (use2ndderiv)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff_, derxy2_);
        diff_.Scale(diffus_[k]);
      }
      else
        diff_.Clear(); // for security reasons

      // convective and diffusive (if required) part times initial scalar field
      const double conv_ephi0_k = conv_.Dot(ephi0[k]);
      double diff_ephi0_k(0.0);
      if (use2ndderiv) diff_ephi0_k = diff_.Dot(ephi0[k]);

      //----------------------------------------------------------------
      // element matrix: standard Galerkin terms
      //----------------------------------------------------------------
      // transient term
      for (int vi=0; vi<iel; ++vi)
      {
        const double v = fac_*funct_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          massmat(fvi,fui) += v*densfunct_(ui);
        }
      }

      // convective term in convective form
      double vrhs = fac_*conv_ephi0_k;
      // addition to convective temporal rhs term for conservative form
      if (conservative) vrhs += fac_*conint_[k]*mdiv_;
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        rhs[fvi] -= vrhs*funct_(vi);
      }

      // diffusive term
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<iel; ++ui)
        {
          double laplawf(0.0);
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
          rhs[fvi] -= fac_diffus*laplawf*(ephi0[k])(ui);
        }
      }

      // nonlinear migration term
      vrhs = fac_diffus*conint_[k]*valence_[k];
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        rhs[fvi] += vrhs*migconv_(vi);
      }

      //----------------------------------------------------------------
      // element matrix: stabilization terms
      //----------------------------------------------------------------
      // convective stabilization of transient term (in convective form)
      for (int vi=0; vi<iel; ++vi)
      {
        const double v = taufac*conv_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          massmat(fvi,fui) += v*densfunct_(ui);
        }
      }

      // convective stabilization of convective term (in convective form)
      vrhs = taufac*conv_ephi0_k;
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        rhs[fvi] -= vrhs*conv_(vi);
      }

      if (use2ndderiv)
      {
        // convective stabilization of diffusive term (in convective form)
        vrhs = taufac*diff_ephi0_k;
        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          rhs[fvi] += vrhs*conv_(vi);
        }

        // diffusive stabilization of transient term
        // (USFEM assumed here, sign change necessary for GLS)
        for (int vi=0; vi<iel; ++vi)
        {
          const double v = taufac*diff_(vi);
          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<iel; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            massmat(fvi,fui) += v*densfunct_(ui);
          }
        }

        // diffusive stabilization of convective term (in convective form)
        // (USFEM assumed here, sign change necessary for GLS)
        vrhs = taufac*conv_ephi0_k;
        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          rhs[fvi] -= vrhs*diff_(vi);
        }

        // diffusive stabilization of diffusive term
        // (USFEM assumed here, sign change necessary for GLS)
        vrhs = taufac*diff_ephi0_k;
        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          rhs[fvi] += vrhs*diff_(vi);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: standard Galerkin bodyforce term
      //----------------------------------------------------------------
      vrhs = fac_*rhsint;
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        rhs[fvi] += vrhs*funct_(vi);
      }

      //----------------------------------------------------------------
      // element right hand side: stabilization terms
      //----------------------------------------------------------------
      // convective stabilization of bodyforce term
      vrhs = taufac*rhsint;
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        rhs[fvi] += vrhs*conv_(vi);
      }

      // diffusive stabilization of bodyforce term
      // (USFEM assumed here, sign change necessary for GLS)
      if (use2ndderiv)
      {
        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          rhs[fvi] += vrhs*diff_(vi);
        }
      }
    } // loop over each scalar k

    if (numdofpernode_-numscal_== 1) // ELCH
    {
      // we put a dummy mass matrix here in order to have a regular 
      // matrix in the lower right block of the whole system-matrix
      // A identity matrix would cause problems with ML solver in the SIMPLE 
      // schemes since ML needs to have off-diagonal entries for the aggregation!
      for (int vi=0; vi<iel; ++vi)
      {
        const double v = fac_*funct_(vi);
        const int fvi = vi*numdofpernode_+numscal_;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui = ui*numdofpernode_+numscal_;

          massmat(fvi,fui) += v*densfunct_(ui);
        }
      }
      // dof for el. potential have no 'velocity' -> rhs is zero!
    }

  } // integration loop

  return;
} // ScaTraImpl::InitialTimeDerivative


/*----------------------------------------------------------------------*
 | calculate normalized subgrid-diffusivity matrix              vg 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcSubgridDiffMatrix(
    const DRT::Element*           ele,
    Epetra_SerialDenseMatrix&     sys_mat_sd,
    const double                  timefac,
    const bool                    is_stationary
    )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

/*----------------------------------------------------------------------*/
// integration loop for one element
/*----------------------------------------------------------------------*/
// integrations points and weights
DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

// integration loop
for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
{
  EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele->Id());

  for (int k=0;k<numscal_;++k)
  {
    // parameter for artificial diffusivity (scaled to one here)
    double kartfac = fac_;
    if (not is_stationary) kartfac *= timefac;

    for (int vi=0; vi<iel; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<iel; ++ui)
      {
        const int fui = ui*numdofpernode_+k;
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
        sys_mat_sd(fvi,fui) += kartfac*laplawf;

        /*subtract SUPG term */
        //sys_mat_sd(fvi,fui) -= taufac*conv(vi)*conv(ui);
      }
    }
  }
} // integration loop

return;
} // ScaTraImpl::CalcSubgridDiffMatrix


/*----------------------------------------------------------------------*
 | calculate matrix and rhs vector (incremental condif form)  gjb 11/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatInc(
    Epetra_SerialDenseMatrix&             emat,
    Epetra_SerialDenseVector&             erhs,
    const vector<LINALG::Matrix<iel,1> >& ephinp,
    const bool                            use2ndderiv,
    const bool                            is_stationary,
    const double                          timefac
)
{
  // get values of all transported scalars at integration point
  for (int k=0; k<numscal_; ++k)
  {
    conint_[k] = funct_.Dot(ephinp[k]);
  }

  // convective part
  /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
  conv_.MultiplyTN(derxy_,velint_);

#if 0
  // DEBUG output
  cout<<endl<<"values at GP:"<<endl;
  cout<<"factor F/RT = "<<frt<<endl;
  for (int k=0;k<numscal_;++k)
  {cout<<"conint_["<<k<<"] = "<<conint_[k]<<endl;}
#endif
  // some 'working doubles'
  double rhsint(0.0);  // rhs at int. point
  // integration factors and coefficients of single terms
  double timefacfac(0.0);
  double timetaufac(0.0);
  double taufac(0.0);

  for (int k = 0; k < numscal_;++k) // loop over all transported sclars
  {
    // stabilization parameters
    taufac = tau_[k]*fac_;

    if (is_stationary)
    {
      timefacfac  = fac_;
      timetaufac  = taufac;
      rhsint = rhs_[k];     // set rhs at integration point
    }
    else
    {
      timefacfac  = timefac * fac_;
      timetaufac  = timefac * taufac;
      rhsint = hist_[k] + rhs_[k]*timefac;     // set rhs at integration point
    }

    // compute gradient of scalar k at integration point
    gradphi_.Multiply(derxy_,ephinp[k]);

    // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
    if (use2ndderiv)
    {
      GetLaplacianStrongForm(diff_, derxy2_);
      diff_.Scale(diffus_[k]);
    }

    // ----------------------------------------matrix entries
    for (int vi=0; vi<iel; ++vi)
    {
      const double timetaufac_conv_vi = timetaufac*conv_(vi);
      const double timefacfac_funct_vi = timefacfac*funct_(vi);

      for (int ui=0; ui<iel; ++ui)
      {
        /* Standard Galerkin terms: */
        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac_funct_vi*conv_(ui) ;

        /* diffusive term */
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac*diffus_[k]*laplawf;

        /* Stabilization term: */
        /* 0) transient stabilization */
        // not implemented

        /* 1) convective stabilization */

        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac_conv_vi*conv_(ui);

      } // for ui
    } // for vi

    if (use2ndderiv)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        for (int ui=0; ui<iel; ++ui)
        {
          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += -timetaufac*conv_(vi)*diff_(ui) ;

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

          /* convective term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac*diff_(vi)*conv_(ui);

          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) -= timetaufac*diff_(vi)*diff_(ui) ;
        } // for ui
      } // for vi

    } // use2ndderiv

    // ----------------------------------------------RHS
    const double conv_ephinp_k = conv_.Dot(ephinp[k]);
    const double densfunct_ephinp_k = densfunct_.Dot(ephinp[k]);
    double diff_ephinp_k(0.0);
    if (use2ndderiv) diff_ephinp_k = diff_.Dot(ephinp[k]); // only necessary for higher order ele!

    // compute residual of strong form for stabilization
    double taufacresidual = taufac*rhsint - timetaufac*(conv_ephinp_k + diff_ephinp_k);
    if (!is_stationary) // add transient term to the residual
      taufacresidual -= taufac*densfunct_ephinp_k;

    //------------residual formulation (Newton iteration)
    for (int vi=0; vi<iel; ++vi)
    {
      // RHS source term
      erhs[vi*numdofpernode_+k] += fac_*funct_(vi)*rhsint ;

      // convective term
      erhs[vi*numdofpernode_+k] -= timefacfac*funct_(vi)*conv_ephinp_k;

      // diffusive term
      erhs[vi*numdofpernode_+k] -= timefacfac*diffus_[k]*(gradphi_(0)*derxy_(0, vi) + gradphi_(1)*derxy_(1, vi)+ gradphi_(2)*derxy_(2, vi));

      // Stabilization terms:

      // 0) transient stabilization
      // not implemented

      // 1) convective stabilization
      erhs[vi*numdofpernode_+k] += conv_(vi) * taufacresidual;

    } // for vi

    if (use2ndderiv)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        // 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS)
        erhs[vi*numdofpernode_+k] += diff_(vi)*taufacresidual ;

      } // for vi
    } // use2ndderiv

    // -----------------------------------INSTATIONARY TERMS
    if (!is_stationary)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const double fac_funct_vi = fac_*funct_(vi);
        for (int ui=0; ui<iel; ++ui)
        {
          /* Standard Galerkin terms: */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += fac_funct_vi*densfunct_(ui) ;

          /* 1) convective stabilization */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*conv_(vi)*densfunct_(ui);

          if (use2ndderiv)
          {
            /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
            /* transient term */
            emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*diff_(vi)*densfunct_(ui) ;
          }
        } // for ui

        // residuum on RHS:

        /* Standard Galerkin terms: */
        /* transient term */
        erhs[vi*numdofpernode_+k] -= fac_funct_vi*densfunct_ephinp_k;

      } // for vi
    } // instationary case

  } // loop over scalars

  return;
} // ScaTraImpl::CalMatInc


/*----------------------------------------------------------------------*
 | calculate matrix and rhs for electrochemistry problem      gjb 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatElch(
    Epetra_SerialDenseMatrix&             emat,
    Epetra_SerialDenseVector&             erhs,
    const vector<LINALG::Matrix<iel,1> >& ephinp,
    const LINALG::Matrix<iel,1>&          epotnp,
    const bool                            use2ndderiv,
    const double                          frt,
    const bool                            is_stationary,
    const double                          timefac
)
{
  // get values of all transported scalars at integration point
  for (int k=0; k<numscal_; ++k)
  {
    conint_[k] = funct_.Dot(ephinp[k]);

    // when concentration becomes zero, the coupling terms in the system matrix get lost!
    if (conint_[k] < 1e-18) 
      dserror("species concentration %d at GP is zero or negative: %g",k,conint_[k]);
  }

  // get gradient of el. potential at integration point
  gradpot_.Multiply(derxy_,epotnp);

  // convective part
  /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
  conv_.MultiplyTN(derxy_,velint_);

  // migration term (convective part)
  migconv_.MultiplyTN(-frt,derxy_,gradpot_);

  // Laplacian of shape functions at integration point
  if (use2ndderiv)
  {
    GetLaplacianStrongForm(laplace_, derxy2_);
  }

#if 0
  // DEBUG output
  cout<<endl<<"values at GP:"<<endl;
  cout<<"factor F/RT = "<<frt<<endl;
  for (int k=0;k<numscal_;++k)
  {cout<<"conint_["<<k<<"] = "<<conint_[k]<<endl;}
  for (int k=0;k<3;++k)
  {cout<<"gradpot_["<<k<<"] = "<<gradpot_[k]<<endl;}
#endif
  // some 'working doubles'
  double diffus_valence_k;
  double rhsint;  // rhs at int. point
  // integration factors and coefficients of single terms
  double timefacfac;
  double timetaufac;
  double taufac;

  for (int k = 0; k < numscal_;++k) // loop over all transported sclars
  {
    // stabilization parameters
    taufac = tau_[k]*fac_;

    if (is_stationary)
    {
      timefacfac  = fac_;
      timetaufac  = taufac;
      rhsint = rhs_[k];     // set rhs at integration point
    }
    else
    {
      timefacfac  = timefac * fac_;
      timetaufac  = timefac * taufac;
      rhsint = hist_[k] + rhs_[k]*timefac;     // set rhs at integration point
    }

    // compute gradient of scalar k at integration point
    gradphi_.Multiply(derxy_,ephinp[k]);

    // factor D_k * z_k
    diffus_valence_k = diffus_[k]*valence_[k];

    if (use2ndderiv)
    {
      // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
      diff_.Update(diffus_[k],laplace_);

      // get Laplacian of el. potential at integration point
      double lappot = laplace_.Dot(epotnp);
      // reactive part of migration term
      migrea_.Update(-frt*diffus_valence_k*lappot,funct_);
    }

    const double frt_timefacfac_diffus_valence_k_conint_k = frt*timefacfac*diffus_valence_k*conint_[k];

    // ----------------------------------------matrix entries
    for (int vi=0; vi<iel; ++vi)
    {
      const int    fvi = vi*numdofpernode_+k;
      const double timetaufac_conv_eff_vi = timetaufac*(conv_(vi)+diffus_valence_k*migconv_(vi));
      const double timefacfac_funct_vi = timefacfac*funct_(vi);
      const double timefacfac_diffus_valence_k_mig_vi = timefacfac*diffus_valence_k*migconv_(vi);
      const double valence_k_fac_funct_vi = valence_[k]*fac_*funct_(vi);

      for (int ui=0; ui<iel; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        /* Standard Galerkin terms: */
        /* convective term */
        emat(fvi, fui) += timefacfac_funct_vi*conv_(ui) ;

        /* diffusive term */
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
        emat(fvi, fui) += timefacfac*diffus_[k]*laplawf;

        /* migration term (directional derivatives) */
        emat(fvi, fui) -= timefacfac_diffus_valence_k_mig_vi*funct_(ui);
        emat(fvi,ui*numdofpernode_+numscal_) += frt_timefacfac_diffus_valence_k_conint_k*laplawf;

        /* electroneutrality condition */
        emat(vi*numdofpernode_+numscal_, fui) += valence_k_fac_funct_vi*densfunct_(ui);

        /* Stabilization term: */
        /* 0) transient stabilization */
        // not implemented

        /* 1) convective stabilization */

        /* convective term */
        emat(fvi, fui) += timetaufac_conv_eff_vi*(conv_(ui)+diffus_valence_k*migconv_(ui));

      } // for ui
    } // for vi

    if (use2ndderiv)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        const double timetaufac_conv_eff_vi = timetaufac*(conv_(vi)+diffus_valence_k*migconv_(vi));
        for (int ui=0; ui<iel; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          /* 1) convective stabilization */

          // diffusive term 
          emat(fvi, fui) += -timetaufac_conv_eff_vi*diff_(ui) ;

          // migration term (reactive part)
          emat(fvi, fui) += -timetaufac_conv_eff_vi*migrea_(ui) ;

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

          // convective term
          emat(fvi, fui) += timetaufac*diff_(vi)*(conv_(ui)+diffus_valence_k*migconv_(ui));

          // diffusive term
          emat(fvi, fui) -= timetaufac*diff_(vi)*diff_(ui) ;

          // migration term (reactive part)
          emat(fvi, fui) += timetaufac*diff_(vi)*migrea_(ui) ;

          /* 3) reactive stabilization (reactive part of migration term) */

          // convective terms
          //emat(fvi, ui*numdofpernode_+k) += timetaufac*migrea_(vi)*(conv_(ui)+diffus_valence_k*migconv_(ui));

          // diffusive term
          //emat(fvi, ui*numdofpernode_+k) -= timetaufac*migrea_(vi)*diff_(ui) ;

          // migration term (reactive part)
          //emat(fvi, ui*numdofpernode_+k) += timetaufac*migrea_(vi)*migrea_(ui) ;

        } // for ui
      } // for vi

    } // use2ndderiv

    // ----------------------------------------------RHS
    const double conv_ephinp_k = conv_.Dot(ephinp[k]);
    const double Dkzk_mig_ephinp_k = diffus_valence_k*(migconv_.Dot(ephinp[k]));
    const double conv_eff_k = conv_ephinp_k + Dkzk_mig_ephinp_k;
    const double densfunct_ephinp_k = densfunct_.Dot(ephinp[k]);
    double diff_ephinp_k(0.0);
    double migrea_k(0.0);
    if (use2ndderiv) 
    { // only necessary for higher order elements
      diff_ephinp_k = diff_.Dot(ephinp[k]);   // diffusion
      migrea_k      = migrea_.Dot(ephinp[k]); // reactive part of migration term
    }

    // compute residual of strong form for residual-based stabilization
    double taufacresidual = taufac*rhsint - timetaufac*(conv_eff_k + diff_ephinp_k + migrea_k);
    if (!is_stationary) // add transient term to the residual
      taufacresidual -= taufac*densfunct_ephinp_k;

    //------------residual formulation (Newton iteration)
    for (int vi=0; vi<iel; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      // RHS source term
      erhs[fvi] += fac_*funct_(vi)*rhsint ;

      // nonlinear migration term
      erhs[fvi] += conint_[k]*timefacfac*diffus_valence_k*migconv_(vi);

      // convective term
      erhs[fvi] -= timefacfac*funct_(vi)*conv_ephinp_k;

      // diffusive term
      erhs[fvi] -= timefacfac*diffus_[k]*(gradphi_(0)*derxy_(0, vi) + gradphi_(1)*derxy_(1, vi)+ gradphi_(2)*derxy_(2, vi));

      // electroneutrality condition
      // for incremental formulation, there is the residuum on the rhs! : 0-ENC*phi_i
      erhs[vi*numdofpernode_+numscal_] -= valence_[k]*fac_*funct_(vi)*densfunct_ephinp_k;

      // Stabilization terms:

      // 0) transient stabilization
      // not implemented

      // 1) convective stabilization
      erhs[fvi] += (conv_(vi)+diffus_valence_k*migconv_(vi)) * taufacresidual;

    } // for vi

    if (use2ndderiv)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        // 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS)
        erhs[fvi] += diff_(vi)*taufacresidual ;

        /* 3) reactive stabilization (reactive part of migration term) */

      } // for vi
    } // use2ndderiv

    // -----------------------------------INSTATIONARY TERMS
    if (!is_stationary)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        const double fac_funct_vi = fac_*funct_(vi);
        for (int ui=0; ui<iel; ++ui)
        {
         const int fui = ui*numdofpernode_+k;

          /* Standard Galerkin terms: */
          /* transient term */
          emat(fvi, fui) += fac_funct_vi*densfunct_(ui) ;

          /* 1) convective stabilization */
          /* transient term */
          emat(fvi, fui) += taufac*(conv_(vi)+diffus_valence_k*migconv_(vi))*densfunct_(ui);

          if (use2ndderiv)
          {
            /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
            /* transient term */
            emat(fvi, fui) += taufac*diff_(vi)*densfunct_(ui);
          }
        } // for ui

        // residuum on RHS:

        /* Standard Galerkin terms: */
        /* transient term */
        erhs[fvi] -= fac_funct_vi*densfunct_ephinp_k;

      } // for vi
    } // instationary case

  } // loop over scalars

  return;
} // ScaTraImpl::CalMatElch


/*---------------------------------------------------------------------*
 |  calculate error compared to analytical solution           gjb 10/08|
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element*            ele,
    ParameterList&                 params,
    const LINALG::Matrix<iel,2>&   ephinp,
    const LINALG::Matrix<iel,1>&   epotnp,
    Epetra_SerialDenseVector&      errors,
    struct _MATERIAL*              material
)
{
  //at the moment, there is only one analytical test problem available!
  if (params.get<string>("action") != "calc_elch_kwok_error")
    dserror("Unknown analytical solution");

  //------------------------------------------------- Kwok et Wu,1995
  //   Reference:
  //   Kwok, Yue-Kuen and Wu, Charles C. K.
  //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
  //   Numerical Methods for Partial Differential Equations
  //   1995, Vol 11, 389-397

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // set constants for analytical solution
  const double t = params.get<double>("total time");
  const double frt = params.get<double>("frt");

  // get material constants
  GetMaterialParams(material,false);

  // working arrays
  double                  potint;
  LINALG::Matrix<2,1>     conint;
  LINALG::Matrix<nsd_,1>  xint;
  LINALG::Matrix<2,1>     c;
  double                  deltapot;
  LINALG::Matrix<2,1>     deltacon(true);

  // integrations points and weights
  // more GP than usual due to cos/exp fcts in analytical solution
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele->Id());

    // get both concentration solutions at integration point
    conint.MultiplyTN(ephinp,funct_);

    // get el. potential solution at integration point
    potint = funct_.Dot(epotnp);

    // get global coordinate of integration point
    xint.Multiply(xyze_,funct_);

    // compute various constants
    const double d = frt*((diffus_[0]*valence_[0]) - (diffus_[1]*valence_[1]));
    if (abs(d) == 0.0) dserror("division by zero");
    const double D = frt*((valence_[0]*diffus_[0]*diffus_[1]) - (valence_[1]*diffus_[1]*diffus_[0]))/d;

    // compute analytical solution for cation and anion concentrations
    const double A0 = 5.0;
    const double m = 2.0;
    const double n = 2.0;
    const double k = 2.0;
    const double A_mnk = 1.0;
    double expterm;
    double c_0_0_0_t;

    if (nsd_==3)
    {
      expterm = exp((-D)*(m*m + n*n + k*k)*t*PI*PI);
      c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1))*cos(k*PI*xint(2)))*expterm);
      c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n + k*k)*t*PI*PI));
    }
    else if (nsd_==2)
    {
      expterm = exp((-D)*(m*m + n*n)*t*PI*PI);
      c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1)))*expterm);
      c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n)*t*PI*PI));
    }
    else if (nsd_==1)
    {
      expterm = exp((-D)*(m*m)*t*PI*PI);
      c(0) = A0 + (A_mnk*(cos(m*PI*xint(0)))*expterm);
      c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m)*t*PI*PI));
    }
    else
      dserror("Illegal number of space dimenions for analyt. solution: %d",nsd_);

    // compute analytical solution for anion concentration
    c(1) = (-valence_[0]/valence_[1])* c(0);
    // compute analytical solution for el. potential
    const double pot = ((diffus_[1]-diffus_[0])/d) * log(c(0)/c_0_0_0_t);

    // compute differences between analytical solution and numerical solution
    deltapot = potint - pot;
    deltacon.Update(1.0,conint,-1.0,c);

    // add square to L2 error
    errors[0] += deltacon(0)*deltacon(0)*fac_; // cation concentration
    errors[1] += deltacon(1)*deltacon(1)*fac_; // anion concentration
    errors[2] += deltapot*deltapot*fac_; // electric potential in electrolyte solution

  } // end of loop over integration points

  return;
} // ScaTraImpl::CalErrorComparedToAnalytSolution


/*----------------------------------------------------------------------*
 |  calculate mass flux                              (private) gjb 06/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateFlux(
    LINALG::Matrix<3,iel>&          flux,
    const DRT::Element*             ele,
    const vector<double>&           ephinp,
    struct _MATERIAL*               material,
    const bool                      temperature,
    const double                    frt,
    const Epetra_SerialDenseVector& evel,
    const SCATRA::FluxType          fluxtype,
    const int                       dofindex
)
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // get material constants
  double diffus(0.0);
  double valence(0.0);
  double diffus_valence_frt(0.0);

  //GetMaterialParams(material,temperature);

  if (material->mattyp == m_matlist)
  {
    const int matid = material->m.matlist->matids[dofindex];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_condif)
      diffus = singlemat.m.condif->diffusivity;
    else if (singlemat.mattyp == m_ion)
    {
      diffus = singlemat.m.ion->diffusivity;
      valence = singlemat.m.ion->valence;
      diffus_valence_frt = diffus*valence*frt;
    }
    else
      dserror("type of material found in material list is not supported.");
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material"); // paranoia?
    diffus = material->m.condif->diffusivity;
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------- declaration of variables ---*/
  LINALG::SerialDenseMatrix nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  if ((int) nodecoords.N() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // reference coordinates of the current node
    for (int idim=0;idim<nsd_;idim++)
      {xsi_(idim) = nodecoords(idim, iquad);}

    // first derivatives
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

    // compute Jacobian matrix and determinant
    // actually compute its transpose....
    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xij_.Invert(xjm_);

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(), det);

    // compute global derivatives
    derxy_.Multiply(xij_,deriv_);

    // gradient of electric potential
    gradpot_.Clear();
    if (frt > 0.0) // ELCH
    {
      for (int k=0;k<iel;k++)
      {
        for (int idim=0; idim<nsd_ ;idim++)
        {
          gradpot_(idim) += derxy_(idim,k)*ephinp[k*numdofpernode_+numscal_];
        }
      }
    }

    const double ephinpatnode = ephinp[iquad*numdofpernode_+dofindex];
    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case SCATRA::totalflux:
        if (frt > 0.0) // ELCH
        {
          //migration flux terms
          for (int idim=0; idim<nsd_ ;idim++)
          {
            flux(idim,iquad) += diffus_valence_frt*gradpot_(idim)*ephinpatnode;
          }
        }
        //convective flux terms
        for (int idim=0; idim<nsd_ ;idim++)
        {
          flux(idim,iquad) -= evel[idim+iquad*nsd_]*ephinpatnode;
        }
        // no break statement here!
      case SCATRA::diffusiveflux:
        //diffusive flux terms
        for (int k=0;k<iel;k++)
        {
          for (int idim=0; idim<nsd_ ;idim++)
          {
            flux(idim,iquad) += diffus*derxy_(idim,k)*ephinp[k*numdofpernode_+dofindex];
          }
        }
        break;
      case SCATRA::noflux:
        dserror("received noflux flag inside flux evaluation");
    };

    //set zeros for unused space dimenions
    for (int idim=nsd_; idim<3; idim++)
    {
      flux(idim,iquad) = 0.0;
    }
  } // loop over nodes

  return;
} // ScaTraImpl::CalculateFlux


/*----------------------------------------------------------------------*
 |  calculate temperature, density and domain integral          vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateTempAndDens(
    const DRT::Element*             ele,
    const vector<double>&           ephinp,
    const vector<double>&           edensnp,
    Epetra_SerialDenseVector&       scalars
)
{
  /*------------------------------------------------- set element data */
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // integrations points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele->Id());

    // calculate integrals of temperature or concentrentations, 
    // then of density and domain
    for (int i=0; i<iel; i++)
    {
      const double fac_funct_i = fac_*funct_(i);
      for (int k = 0; k < numscal_; k++)
      {
        scalars[k] += fac_funct_i*ephinp[i*numdofpernode_+k];
      }
      scalars[numscal_]    += fac_funct_i*edensnp[i];
      scalars[numscal_+1]  += fac_funct_i;
    }
  } // loop over integration points

  return;
} // ScaTraImpl::CalculateTempAndDens


/*----------------------------------------------------------------------*
 |  calculate domain integral                                   vg 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateDomainAndBodyforce(
    Epetra_SerialDenseVector&  scalars,
    const DRT::Element*        ele,
    const double               time
)
{
  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele,time);

  /*------------------------------------------------- set element data */
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // integrations points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele->Id());

    // get bodyforce in gausspoint
    rhs_[0] = bodyforce_[0].Dot(funct_);

    // calculate integrals of domain and bodyforce
    for (int i=0; i<iel; i++)
    {
      scalars[0] += fac_*funct_(i);
    }
    scalars[1] += fac_*rhs_[0];

  } // loop over integration points

  return;
} // ScaTraImpl::CalculateDomain


#endif // CCADISCRET
#endif // D_FLUID3 or D_FLUID2
