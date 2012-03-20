/*----------------------------------------------------------------------*/
/*!
\file elasthyper.cpp
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 0 MATIDS  DENS 0 GAMMA 0 INIT_MODE -1

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "elasthyper.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElastHyper::ElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  density_(matdata->GetDouble("DENS")),
  gamma_(matdata->GetDouble("GAMMA")),
  init_mode_(matdata->GetInt("INIT_MODE"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=matids_->begin(); m!=matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> potsum = MAT::ELASTIC::Summand::Factory(matid);
    if (potsum == Teuchos::null) dserror("Failed to allocate");
    potsum_.insert(std::pair<int,Teuchos::RCP<MAT::ELASTIC::Summand> >(matid,potsum));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElastHyper(this));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const MAT::ELASTIC::Summand> MAT::PAR::ElastHyper::MaterialById(const int id) const
{
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::const_iterator m = potsum_.find(id);
  if (m == potsum_.end())
  {
    dserror("Material %d could not be found", id);
    return Teuchos::null;
  }
  else
  {
    return m->second;
  }
}


MAT::ElastHyperType MAT::ElastHyperType::instance_;


DRT::ParObject* MAT::ElastHyperType::Create( const std::vector<char> & data )
{
  MAT::ElastHyper* elhy = new MAT::ElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 08/09|
 *----------------------------------------------------------------------*/
// 6-Voigt C-index                              0 1 2  3 4 5
const int MAT::ElastHyper::VOIGT6ROW_[6] = {0,1,2, 0,1,2};
const int MAT::ElastHyper::VOIGT6COL_[6] = {0,1,2, 1,2,0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
const int MAT::ElastHyper::VOIGT3X3SYM_[9] = {0,3,5, 3,1,4, 5,4,2};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
  AddtoPack(data,isoprinc_);
  AddtoPack(data,isomod_);
  AddtoPack(data,anisoprinc_);
  AddtoPack(data,anisomod_);
  AddtoPack(data,a1_);
  AddtoPack(data,a2_);
  AddtoPack(data,A1_);
  AddtoPack(data,A2_);
  AddtoPack(data,A1A2_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->PackSummand(data);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Unpack(const std::vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }
  }

  isoprinc_ = false;
  isomod_ = false;
  anisoprinc_ = false;
  anisomod_ = false;

  int isoprinc;
  int isomod;
  int anisoprinc;
  int anisomod;

  ExtractfromPack(position,data,isoprinc);
  ExtractfromPack(position,data,isomod);
  ExtractfromPack(position,data,anisoprinc);
  ExtractfromPack(position,data,anisomod);

  if (isoprinc != 0) isoprinc_ = true;
  if (isomod != 0) isomod_ = true;
  if (anisoprinc != 0) anisoprinc_ = true;
  if (anisomod != 0) anisomod_ = true;

  ExtractfromPack(position,data,a1_);
  ExtractfromPack(position,data,a2_);
  ExtractfromPack(position,data,A1_);
  ExtractfromPack(position,data,A2_);
  ExtractfromPack(position,data,A1A2_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->UnpackSummand(data,position);
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ElastHyper::MatID(
  const unsigned index
) const
{
  if ((int)index < params_->nummat_)
    return params_->matids_->at(index);
  else
  {
    dserror("Index too large");
    return -1;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElastHyper::ShearMod() const
{
  // principal coefficients
  bool haveshearmod = false;
  double shearmod = 0.0;
  {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddShearMod(haveshearmod,shearmod);
    }
  }

  if (haveshearmod)
  {
    return shearmod;
  }
  else
  {
    dserror("Cannot provide shear modulus equivalent");
    return -1.0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::SetupAAA(Teuchos::ParameterList& params)
{
  // loop map of associated potential summands
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
  for (p=pot.begin(); p!=pot.end(); ++p)
  {
    p->second->SetupAAA(params);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // Setup summands
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
  for (p=pot.begin(); p!=pot.end(); ++p)
  {
    p->second->Setup(linedef);
  }

  // find out which formulations are used

  isoprinc_ = false ;
  isomod_ = false ;
  anisoprinc_ = false ;
  anisomod_ = false;

  for (p=pot.begin(); p!=pot.end(); ++p)
  {
    p->second->SpecifyFormulation(isoprinc_,isomod_,anisoprinc_,anisomod_);
  }


  if (anisoprinc_ or anisomod_)
  {
    // fibers aligned in local element cosy with gamma_i around circumferential direction
    // -> check whether element supports local element cosy
    vector<double> rad;
    vector<double> axi;
    vector<double> cir;
    if (linedef->HaveNamed("RAD") and
        linedef->HaveNamed("AXI") and
        linedef->HaveNamed("CIR"))
    {
      // read local (cylindrical) cosy-directions at current element
      // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
      LINALG::Matrix<3,3> locsys(true);
      linedef->ExtractDoubleVector("RAD",rad);
      linedef->ExtractDoubleVector("AXI",axi);
      linedef->ExtractDoubleVector("CIR",cir);
      double radnorm=0.; double axinorm=0.; double cirnorm=0.;

      for (int i = 0; i < 3; ++i)
      {
        radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
      }
      radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);

      for (int i=0; i<3; ++i)
      {
        locsys(i,0) = rad[i]/radnorm;
        locsys(i,1) = axi[i]/axinorm;
        locsys(i,2) = cir[i]/cirnorm;
      }
      // INIT_MODE = 0 : Fiber direction derived from local cosy
      if( 0 == params_->init_mode_)
      {
        // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
        if ((params_->gamma_<0) || (params_->gamma_ >90)) dserror("Fiber angle not in [0,90]");
        //convert
        const double gamma = (params_->gamma_*PI)/180.;

        for (int i = 0; i < 3; ++i)
        {
          // a1 = cos gamma e3 + sin gamma e2
          a1_(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
          // a2 = cos gamma e3 - sin gamma e2
          a2_(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
        }
      }
      // INIT_MODE = 1 : Fiber direction aligned to local cosy
      else if (1 == params_->init_mode_)
      {
        for (int i = 0; i < 3; ++i)
        {
          a1_(i) = locsys(i,0);
          a2_(i) = locsys(i,1);
        }
      }
      // INIT_MODE = -1 = default value; usage of fiber direction without initialization mode
      else if (-1 == params_->init_mode_)
      {
        dserror("Forgotten to give INIT_MODE in .dat-file");
      }
      else
      {
        dserror("Problem with fiber initialization");
      }
      for (int i = 0; i < 3; ++i) {
        A1_(i) = a1_(i)*a1_(i);
        A2_(i) = a2_(i)*a2_(i);
        for (int j=0; j<3; j++)
        {
          A1A2_(j,i) = a1_(j)*a2_(i);
        }
      }
      A1_(3) = a1_(0)*a1_(1); A1_(4) = a1_(1)*a1_(2); A1_(5) = a1_(0)*a1_(2);
      A2_(3) = a2_(0)*a2_(1); A2_(4) = a2_(1)*a2_(2); A2_(5) = a2_(0)*a2_(2);
    }

    else
    {
      dserror("Reading of element local cosy for anisotropic materials failed");
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsPrincipal(
    LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<6,1>& rcg)
{
  // 1st invariant, trace
  prinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  prinv(1) = 0.5*( prinv(0)*prinv(0)
                   - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
                   - .5*rcg(3)*rcg(3) - .5*rcg(4)*rcg(4) - .5*rcg(5)*rcg(5) );
  // 3rd invariant, determinant
  prinv(2) = rcg(0)*rcg(1)*rcg(2)
    + 0.25 * rcg(3)*rcg(4)*rcg(5)
    - 0.25 * rcg(1)*rcg(5)*rcg(5)
    - 0.25 * rcg(2)*rcg(3)*rcg(3)
    - 0.25 * rcg(0)*rcg(4)*rcg(4);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsModified(
    LINALG::Matrix<3,1>& modinv,  ///< modified invariants
    const LINALG::Matrix<3,1>& prinv  ///< principal invariants
    )
{
  // 1st invariant, trace
  modinv(0) = prinv(0)*std::pow(prinv(2),-1./3.);
  // 2nd invariant
  modinv(1) = prinv(1)*std::pow(prinv(2),-2./3.);
  // J
  modinv(2) = std::pow(prinv(2),1./2.);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsPrincipalAniso(
          LINALG::Matrix<6,1>& pranisoinv, ///< principal invariants
          const LINALG::Matrix<6,1>& rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
          )
{
  // 1st invariant, trace
  pranisoinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  pranisoinv(1) = 0.5*( pranisoinv(0)*pranisoinv(0)
                   - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
                   - .5*rcg(3)*rcg(3) - .5*rcg(4)*rcg(4) - .5*rcg(5)*rcg(5) );
  // 3rd invariant, determinant
  pranisoinv(2) = rcg(0)*rcg(1)*rcg(2)
    + 0.25 * rcg(3)*rcg(4)*rcg(5)
    - 0.25 * rcg(1)*rcg(5)*rcg(5)
    - 0.25 * rcg(2)*rcg(3)*rcg(3)
    - 0.25 * rcg(0)*rcg(4)*rcg(4);

  pranisoinv(3) =  A1_(0)*rcg(0) + A1_(1)*rcg(1) + A1_(2)*rcg(2)
            + A1_(3)*rcg(3) + A1_(4)*rcg(4) + A1_(5)*rcg(5);

  pranisoinv(4) =  A2_(0)*rcg(0) + A2_(1)*rcg(1) + A2_(2)*rcg(2)
            + A2_(3)*rcg(3) + A2_(4)*rcg(4) + A2_(5)*rcg(5);

  pranisoinv(5) =  A1A2_(0,0)*rcg(0) + A1A2_(1,1)*rcg(1) + A1A2_(2,2)*rcg(2)
            + 0.5*(A1A2_(0,1)*rcg(3) + A1A2_(1,2)*rcg(4) + A1A2_(0,2)*rcg(5))
            + 0.5*(A1A2_(1,0)*rcg(3) + A1A2_(2,1)*rcg(4) + A1A2_(2,0)*rcg(5));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StretchesPrincipal(
  LINALG::Matrix<3,1>& prstr,
  LINALG::Matrix<3,3>& prdir,
  const LINALG::Matrix<6,1>& rcg
  )
{
  // create right Cauchy-Green 2-tensor
  LINALG::Matrix<3,3> rcgt(false);
  rcgt(0,0) = rcg(0);
  rcgt(1,1) = rcg(1);
  rcgt(2,2) = rcg(2);
  rcgt(0,1) = rcgt(1,0) = 0.5*rcg(3);
  rcgt(1,2) = rcgt(2,1) = 0.5*rcg(4);
  rcgt(2,0) = rcgt(0,2) = 0.5*rcg(5);

  // eigenvalue decomposition
  LINALG::Matrix<3,3> prstr2;  // squared principal stretches
  LINALG::SYEV(rcgt,prstr2,prdir);

  // THE principal stretches
  for (int al=0; al<3; ++al) prstr(al) = std::sqrt(prstr2(al,al));

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StretchesModified(
  LINALG::Matrix<3,1>& modstr,
  const LINALG::Matrix<3,1>& prstr
  )
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0)*prstr(1)*prstr(2);

  // determine modified principal stretches
  modstr.Update(std::pow(detdefgrad,-1.0/3.0),prstr);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ElastHyper::HaveCoefficientsStretchesPrincipal()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      havecoeff = havecoeff or p->second->HaveCoefficientsStretchesPrincipal();
    }
  }

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ElastHyper::HaveCoefficientsStretchesModified()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      havecoeff = havecoeff or p->second->HaveCoefficientsStretchesModified();
    }
  }

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  Teuchos::ParameterList& params
  )
{

  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> rcg(true) ;
  LINALG::Matrix<6,1> scg(true) ;
  LINALG::Matrix<6,1> icg(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> modinv(true);
  LINALG::Matrix<6,1> pranisoinv(true);

  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);
  LINALG::Matrix<3,1> anisogamma(true);
  LINALG::Matrix<15,1> anisodelta(true);

  EvaluateKinQuant(glstrain,id2,scg,rcg,icg,id4,id4sharp,prinv,modinv,pranisoinv);
  EvaluateGammaDelta(rcg,prinv,modinv,pranisoinv,gamma,delta,modgamma,moddelta,anisogamma,anisodelta);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress.Clear();
  cmat.Clear();

  // build stress response and elasticity tensor
  // for potentials based on principal invariants
  if (isoprinc_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisoprinc(true) ;
    EvaluateIsotropicPrinc(stressisoprinc,cmatisoprinc,scg,id2,icg,id4sharp,gamma,delta);
    stress.Update(1.0, stressisoprinc, 1.0);
    cmat.Update(1.0,cmatisoprinc,1.0);
  }

  if (isomod_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodiso(true);
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvol(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvol(true) ;
    EvaluateIsotropicMod(stressisomodiso,stressisomodvol,cmatisomodiso,cmatisomodvol,rcg,id2,icg,id4,id4sharp,modinv,prinv,modgamma,moddelta);
    stress.Update(1.0, stressisomodiso, 1.0);
    stress.Update(1.0, stressisomodvol, 1.0);
    cmat.Update(1.0,cmatisomodiso,1.0);
    cmat.Update(1.0,cmatisomodvol,1.0);
  }

  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  const bool havecoeffstrpr = HaveCoefficientsStretchesPrincipal();
  const bool havecoeffstrmod = HaveCoefficientsStretchesModified();
  if (havecoeffstrpr or havecoeffstrmod) {
    ResponseStretches(cmat,stress,rcg,havecoeffstrpr,havecoeffstrmod);
  }

  /*----------------------------------------------------------------------*/
  //Do all the anisotropic stuff!
  if (anisoprinc_)
  {
      LINALG::Matrix<NUM_STRESS_3D,1> stressanisoprinc(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoprinc(true) ;
      EvaluateAnisotropicPrinc(stressanisoprinc,cmatanisoprinc,scg,id2,icg,anisogamma,anisodelta);
      stress.Update(1.0, stressanisoprinc, 1.0);
      cmat.Update(1.0, cmatanisoprinc, 1.0);
  }

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateKinQuant(
    const LINALG::Matrix<6,1>& glstrain,
    LINALG::Matrix<6,1>& id2,
    LINALG::Matrix<6,1>& scg,
    LINALG::Matrix<6,1>& rcg,
    LINALG::Matrix<6,1>& icg,
    LINALG::Matrix<6,6>& id4,
    LINALG::Matrix<6,6>& id4sharp,
    LINALG::Matrix<3,1>& prinv,
    LINALG::Matrix<3,1>& modinv,
    LINALG::Matrix<6,1>& pranisoinv)

{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  rcg.Update(2.0,glstrain,1.0);
  rcg.Update(1.0, id2, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor C^{AB}
  // REMARK: stress-like 6-Voigt vector of right CG
  scg.Update(1.0,rcg,1.0);
  for (int i=3; i<6; i++) scg(i) *= 0.5;

  // principal invariants of right Cauchy-Green strain
  InvariantsPrincipal(prinv,rcg);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    icg(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / prinv(2);
    icg(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / prinv(2);
    icg(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / prinv(2);
    icg(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / prinv(2);
    icg(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / prinv(2);
    icg(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / prinv(2);
  }

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  for (int i=0; i<6; i++) id4(i,i) = 1.0;

  // modified invariants of right Cauchy-Green strain
  InvariantsModified(modinv,prinv);

  // Anisotropic contribution
  // principal anisotropic invariants of right Cauchy-Green strain
  // we only use the set {I1,I2,I3,I4,I6,I8},
  // because the quadratic invariants I5 and I7 don't give any new info, I9 is constant
  InvariantsPrincipalAniso(pranisoinv,rcg);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateGammaDelta(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<6,1> pranisoinv,
    LINALG::Matrix<3,1>& gamma,
    LINALG::Matrix<8,1>& delta,
    LINALG::Matrix<3,1>& modgamma,
    LINALG::Matrix<5,1>& moddelta,
    LINALG::Matrix<3,1>& anisogamma,
    LINALG::Matrix<15,1>& anisodelta
    )

{
  // principal coefficients
  if (isoprinc_)
  {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddCoefficientsPrincipal(gamma,delta,prinv);
    }
  }

  // modified coefficients
  if (isomod_)
  {

    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddCoefficientsModified(modgamma,moddelta,modinv);
    }

  }

  // modified coefficients
  if (anisoprinc_)
  {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p)
    {
      p->second->AddCoefficientsPrincipalAniso(anisogamma,anisodelta,pranisoinv);
    }

    return ;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateIsotropicPrinc(
    LINALG::Matrix<6,1>& stressisoprinc,
    LINALG::Matrix<6,6>& cmatisoprinc,
    LINALG::Matrix<6,1> scg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4sharp,
    LINALG::Matrix<3,1> gamma,
    LINALG::Matrix<8,1> delta
    )
{

  // 2nd Piola Kirchhoff stresses
  stressisoprinc.Update(gamma(0), id2, 1.0);
  stressisoprinc.Update(gamma(1), scg, 1.0);
  stressisoprinc.Update(gamma(2), icg, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmatisoprinc.MultiplyNT(delta(0), id2, id2, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmatisoprinc.MultiplyNT(delta(1), id2, scg, 1.0);
  cmatisoprinc.MultiplyNT(delta(1), scg, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmatisoprinc.MultiplyNT(delta(2), id2, icg, 1.0);
  cmatisoprinc.MultiplyNT(delta(2), icg, id2, 1.0);
  // contribution: C \otimes C
  cmatisoprinc.MultiplyNT(delta(3), scg, scg, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmatisoprinc.MultiplyNT(delta(4), scg, icg, 1.0);
  cmatisoprinc.MultiplyNT(delta(4), icg, scg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmatisoprinc.MultiplyNT(delta(5), icg, icg, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmatisoprinc, icg, delta(6));
  // contribution: Id4^#
  cmatisoprinc.Update(delta(7), id4sharp, 1.0);

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateIsotropicMod(
    LINALG::Matrix<6,1>& stressisomodiso,
    LINALG::Matrix<6,1>& stressisomodvol,
    LINALG::Matrix<6,6>& cmatisomodiso,
    LINALG::Matrix<6,6>& cmatisomodvol,
    LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4,
    LINALG::Matrix<6,6> id4sharp,
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> modgamma,
    LINALG::Matrix<5,1> moddelta
    )
{
  // define necessary variables
  const double modscale = std::pow(prinv(2),-1./3.);
  // modified right Cauchy-Green
  LINALG::Matrix<6,1> modrcg(true);
  modrcg.Update(modscale,rcg);

  // 2nd Piola Kirchhoff stresses

  // isochoric contribution
  LINALG::Matrix<6,1> modstress(true);
  modstress.Update(modgamma(0), id2);
  modstress.Update(modgamma(1), modrcg, 1.0);
  // build 4-tensor for projection as 6x6 tensor
  LINALG::Matrix<6,6> Projection;
  Projection.MultiplyNT(1./3., icg, rcg);
  Projection.Update(1.0, id4, -1.0);
  // isochoric stress
  LINALG::Matrix<6,1> isostress(true);
  stressisomodiso.MultiplyNN(modscale,Projection,modstress,1.0);

  // volumetric contribution
  stressisomodvol.Update(modgamma(2)*modinv(2), icg, 1.0);

  // constitutive tensor

  //isochoric contribution
  // modified constitutive tensor
  LINALG::Matrix<6,6> modcmat;
  LINALG::Matrix<6,6> modcmat2(true);
  // contribution: Id \otimes Id
  modcmat.MultiplyNT(moddelta(0), id2, id2);
  // contribution: Id \otimes C + C \otimes Id
  modcmat.MultiplyNT(moddelta(1), id2, modrcg, 1.0);
  modcmat.MultiplyNT(moddelta(1), rcg, id2, 1.0);
  // contribution: C \otimes C
  modcmat.MultiplyNT(moddelta(2), rcg, modrcg, 1.0);
  // contribution: Id4^#
  modcmat.Update(moddelta(3), id4sharp, 1.0);
  //scaling
  modcmat.Scale(std::pow(modinv(2),-4./3.));
  //contribution: P:modC:P
  modcmat2.MultiplyNN(Projection,modcmat);
  cmatisomodiso.MultiplyNT(1.0,modcmat2,Projection,1.0);
  // contribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
  modcmat.Clear();
  modcmat.MultiplyNT(-1.0/3.0,icg,icg);
  AddtoCmatHolzapfelProduct(modcmat, icg, 1.0);
  LINALG::Matrix<1,1> tracemat;
  tracemat.MultiplyTN(2./3.*std::pow(modinv(2),-2./3.),modstress,rcg);
  cmatisomodiso.Update(tracemat(0,0),modcmat,1.0);
  //contribution: -2/3 (Cinv \otimes S_iso + S_iso \otimes Cinv)
  cmatisomodiso.MultiplyNT(-2./3.,icg,stressisomodiso,1.0);
  cmatisomodiso.MultiplyNT(-2./3.,stressisomodiso,icg,1.0);

  //volumetric contribution
  //contribution: 2 \tilde p Cinv \otimes Cinv
  cmatisomodvol.MultiplyNT(modinv(2)* moddelta(4),icg,icg,1.0);
  //contribution: -2 J*p Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmatisomodvol, icg, -2*modinv(2)*modgamma(2));

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateAnisotropicPrinc(
    LINALG::Matrix<6,1>& stressanisoprinc,
    LINALG::Matrix<6,6>& cmatanisoprinc,
    LINALG::Matrix<6,1> scg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<3,1> anisogamma,
    LINALG::Matrix<15,1> anisodelta
    )
{
  // build Voigt (stress-like) version of a1 \otimes a2 + a2 \otimes a1
  LINALG::Matrix<6,1> A1A2sym;
  A1A2sym(0)=2*A1A2_(0,0);
  A1A2sym(1)=2*A1A2_(1,1);
  A1A2sym(2)=2*A1A2_(2,2);
  A1A2sym(3)=A1A2_(0,1)+A1A2_(1,0);
  A1A2sym(4)=A1A2_(1,2)+A1A2_(2,1);
  A1A2sym(5)=A1A2_(0,2)+A1A2_(2,0);

  // 2nd Piola Kirchhoff stresses
  stressanisoprinc.Update(anisogamma(0), A1_, 1.0);
  stressanisoprinc.Update(anisogamma(1), A2_, 1.0);
  stressanisoprinc.Update(anisogamma(2), A1A2sym, 1.0);

  // constitutive tensor
  // contribution: A1_ \otimes A1_
  cmatanisoprinc.MultiplyNT(anisodelta(0), A1_, A1_, 1.0);
  // contribution: A2_ \otimes A2_
  cmatanisoprinc.MultiplyNT(anisodelta(1), A2_, A2_, 1.0);
  // contribution: A1_ \otimes Id + Id \otimes A1_
  cmatanisoprinc.MultiplyNT(anisodelta(2), A1_, id2, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(2), id2, A1_, 1.0);
  // contribution: A2_ \otimes Id + Id \otimes A2_
  cmatanisoprinc.MultiplyNT(anisodelta(3), A2_, id2, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(3), id2, A2_, 1.0);
  // contribution: A1_ \otimes C + C \otimes A1_
  cmatanisoprinc.MultiplyNT(anisodelta(4), A1_, scg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(4), scg, A1_, 1.0);
  // contribution: A2_ \otimes C + C \otimes A2_
  cmatanisoprinc.MultiplyNT(anisodelta(5), A2_, scg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(5), scg, A2_, 1.0);
  // contribution: A1_ \otimes Cinv + Cinv \otimes A1_
  cmatanisoprinc.MultiplyNT(anisodelta(6), A1_, icg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(6), icg, A1_, 1.0);
  // contribution: A2_ \otimes Cinv + Cinv \otimes A2_
  cmatanisoprinc.MultiplyNT(anisodelta(7), A2_, icg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(7), icg, A2_, 1.0);
  // contribution: A1_ \otimes A2_ + A2_ \otimes A1_
  cmatanisoprinc.MultiplyNT(anisodelta(8), A1_, A2_, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(8), A2_, A1_, 1.0);
  // contribution: A1A2sym \otimes Id + Id \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(9), A1A2sym, id2, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(9), id2, A1A2sym, 1.0);
  // contribution: A1A2sym \otimes C + C \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(10), A1A2sym, scg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(10), scg, A1A2sym, 1.0);
  // contribution: A1A2sym \otimes Cinv + Cinv \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(11), A1A2sym, icg, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(11), icg, A1A2sym, 1.0);
  // contribution: A1A2sym \otimes A1_ + A1_ \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(12), A1A2sym, A1_, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(12), A1_, A1A2sym, 1.0);
  // contribution: A1A2sym \otimes A2_ + A2_ \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(13), A1A2sym, A2_, 1.0);
  cmatanisoprinc.MultiplyNT(anisodelta(13), A2_, A1A2sym, 1.0);
  // contribution: A1A2sym \otimes A1A2sym
  cmatanisoprinc.MultiplyNT(anisodelta(14), A1A2sym, A1A2sym, 1.0);

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::ResponseStretches(
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  const LINALG::Matrix<6,1> rcg,
  const bool& havecoeffstrpr,
  const bool& havecoeffstrmod
  )
{
  // get principal stretches and directions
  LINALG::Matrix<3,1> prstr;
  LINALG::Matrix<3,3> prdir;
  StretchesPrincipal(prstr,prdir,rcg);
  // modified stretches
  LINALG::Matrix<3,1> modstr;
  StretchesModified(modstr,prstr);
  // determinant of deformation gradient
  const double detdefgrad = prstr(0)*prstr(1)*prstr(2);

  // get coefficients
  LINALG::Matrix<3,1> gamma_(true);
  LINALG::Matrix<6,1> delta_(true);
  if (havecoeffstrpr) {
    // loop map of associated potential summands
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
    std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
    for (p=pot.begin(); p!=pot.end(); ++p) {
      p->second->AddCoefficientsStretchesPrincipal(gamma_,delta_,prstr);
    }
  }
  if (havecoeffstrmod) {
    // reciprocal of cubic root of determinant of deformation gradient (convenience)
    const double detdefgrad13 = std::pow(detdefgrad,-1.0/3.0);
    // retrieve coefficients with respect to modified principal stretches
    LINALG::Matrix<3,1> modgamma(true);
    LINALG::Matrix<6,1> moddelta(true);
    {
      // loop map of associated potential summands
      std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
      std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
      for (p=pot.begin(); p!=pot.end(); ++p) {
        p->second->AddCoefficientsStretchesModified(modgamma,moddelta,modstr);
      }
    }
    // convert modified coefficients to oridinary counterparts
    //
    // derivatives of modified pr. stretches WRT pr. stretches
    LINALG::Matrix<3,3> modbypr(false);
    for (int al=0; al<3; ++al) {
      for (int be=0; be<3; ++be) {
        modbypr(al,be) = -modstr(al)/modstr(be);
      }
      modbypr(al,al) += 3.0;
    }
    modbypr.Scale(detdefgrad13/3.0);
    // determine unmodified coefficients gamma and add them
    gamma_.MultiplyTN(1.0,modbypr,modgamma,1.0);
    // determine unmodified coefficients delta and add them
    //
    // rewrite mod.coeff. as 2-tensor
    LINALG::Matrix<3,3> moddeltat(false);
    moddeltat(0,0) = moddelta(0);
    moddeltat(1,1) = moddelta(1);
    moddeltat(2,2) = moddelta(2);
    moddeltat(0,1) = moddeltat(1,0) = moddelta(3);
    moddeltat(1,2) = moddeltat(2,1) = moddelta(4);
    moddeltat(2,0) = moddeltat(0,2) = moddelta(5);
    // Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    LINALG::Matrix<3,3> aux(false);
    aux.MultiplyTN(modbypr,moddeltat);
    LINALG::Matrix<3,3> deltat(false);
    deltat.MultiplyNN(aux,modbypr);
    // Psi_{,barlam} barlam_{,lam lam}
    for (int be=0; be<3; ++be) {
      for (int ga=0; ga<3; ++ga) {
        double deltat_bega = 0.0;
        for (int al=0; al<3; ++al) {
          deltat_bega += -modgamma(al)*modbypr(al,be)/(3.0*prstr(ga));
          if (ga==al)
            deltat_bega += -modgamma(al)*detdefgrad13/(3.0*prstr(be));
          if (be==ga)
            deltat_bega += modgamma(al)*detdefgrad13*prstr(al)/(3.0*prstr(be)*prstr(be));
        }
        deltat(be,ga) += deltat_bega;
      }
    }
    // add to delta
    // Psi_{lam lam} = Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    //               + Psi_{,barlam} barlam_{,lam lam}
    delta_(0) += deltat(0,0);
    delta_(1) += deltat(1,1);
    delta_(2) += deltat(2,2);
    delta_(3) += deltat(0,1);
    delta_(4) += deltat(1,2);
    delta_(5) += deltat(2,0);
  }

  // principal 2nd Piola--Kirchhoff stress tensor, cf [1] Eq (6.47)
  LINALG::Matrix<3,1> prsts(true);
  for (int al=0; al<3; ++al) {
    // PK2 principal stresses
    prsts(al) = gamma_(al)/prstr(al);
    // PK2 tensor in Voigt notation
    stress(0) += prsts(al)*prdir(0,al)*prdir(0,al);  // S^11
    stress(1) += prsts(al)*prdir(1,al)*prdir(1,al);  // S^22
    stress(2) += prsts(al)*prdir(2,al)*prdir(2,al);  // S^33
    stress(3) += prsts(al)*prdir(0,al)*prdir(1,al);  // S^12
    stress(4) += prsts(al)*prdir(1,al)*prdir(2,al);  // S^23
    stress(5) += prsts(al)*prdir(2,al)*prdir(0,al);  // S^31
  }

  // integration factor prfact_{al be}
  LINALG::Matrix<6,1> prfact1(true);
  LINALG::Matrix<6,1> prfact2(true);
  for (int albe=0; albe<6; ++albe) {
    const int al = VOIGT6ROW_[albe];
    const int be = VOIGT6COL_[albe];
    double prfact1_albe = delta_(albe)/(prstr(al)*prstr(be));
    if (albe<3) prfact1_albe -= gamma_(al)/(prstr(be)*prstr(al)*prstr(al));
    prfact1(albe) = prfact1_albe;
    if (al != be) {
      if (fabs(prstr(al)-prstr(be)) < EPS6)
        prfact2(albe) = (prfact1(be) - prfact1(albe))/2.0;
      else
        prfact2(albe) = (prsts(be)-prsts(al))/(prstr(be)*prstr(be)-prstr(al)*prstr(al));
    }
  }

  // add elasticity 4-tensor, cf Holzapfel [1] Eq (6.180),(6.196)
  for (int kl=0; kl<6; ++kl) {
    const int k = VOIGT6ROW_[kl];
    const int l = VOIGT6COL_[kl];
    for (int ij=0; ij<6; ++ij) {
      const int i = VOIGT6ROW_[ij];
      const int j = VOIGT6COL_[ij];
      double c_ijkl = 0.0;
      for (int albe=0; albe<6; ++albe) {
        const int al = VOIGT6ROW_[albe];
        const int be = VOIGT6COL_[albe];
        const double fact1 = prfact1(albe);
        c_ijkl += fact1*prdir(i,al)*prdir(j,al)*prdir(k,be)*prdir(l,be);
        if (albe>=3) { // al!=be
          c_ijkl += fact1*prdir(i,be)*prdir(j,be)*prdir(k,al)*prdir(l,al);
          const double fact2 = prfact2(albe);
          c_ijkl += fact2*prdir(i,al)*prdir(j,be)*prdir(k,al)*prdir(l,be)
                  + fact2*prdir(i,al)*prdir(j,be)*prdir(k,be)*prdir(l,al)
                  + fact2*prdir(i,be)*prdir(j,al)*prdir(k,be)*prdir(l,al)
                  + fact2*prdir(i,be)*prdir(j,al)*prdir(k,al)*prdir(l,be);
        }
      }
      cmat(ij,kl) += c_ijkl;
    }
  }

//  cout << "rcg=" << rcg << endl;
//  cout << "stress=" << stress << endl;
//  cout << "cmat=" << cmat << endl;
//  exit(0);


  // ready
  return;
}

#endif
