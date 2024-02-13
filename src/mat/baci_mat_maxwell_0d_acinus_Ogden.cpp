/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional
acinus elements with non-linear Ogden-like spring, inherits from Maxwell_0d_acinus

The originally linear spring (Stiffness1) of the 4-element Maxwell model is substituted by a
non-linear pressure-volume relation derived from the Ogden strain energy function considering pure
volumetric expansion (derivation: see Christian Roth's dissertation, Appendix B)

Necessary input lines:
(material section)
MAT 3 MAT_0D_MAXWELL_ACINUS_OGDEN Stiffness1 1.0 Stiffness2 5249.1 Viscosity1 3221.86 Viscosity2
1000.0 // acinus properties;
(element section)
1 RED_ACINUS LINE2 2 3 MAT 3 TYPE VolumetricOgden AcinusVolume 300 AlveolarDuctVolume 0.03711 KAPPA
2000.0 BETA -3.1


\level 3

*/
/*----------------------------------------------------------------------*/


#include "baci_mat_maxwell_0d_acinus_Ogden.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_red_airways_elem_params.hpp"
#include "baci_red_airways_elementbase.hpp"

#include <vector>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell_0d_acinus_Ogden::Maxwell_0d_acinus_Ogden(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Maxwell_0d_acinus(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell_0d_acinus_Ogden::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell_0d_acinus_Ogden(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

MAT::Maxwell_0d_acinusOgdenType MAT::Maxwell_0d_acinusOgdenType::instance_;


CORE::COMM::ParObject* MAT::Maxwell_0d_acinusOgdenType::Create(const std::vector<char>& data)
{
  MAT::Maxwell_0d_acinus_Ogden* mxwll_0d_acin = new MAT::Maxwell_0d_acinus_Ogden();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_Ogden::Maxwell_0d_acinus_Ogden() : Maxwell_0d_acinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus_Ogden::Maxwell_0d_acinus_Ogden(MAT::PAR::Maxwell_0d_acinus* params)
    : Maxwell_0d_acinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // Pack type of this instance of ParObject
  int type = UniqueParObjectId();

  AddtoPack(data, type);
  AddtoPack(data, kappa_);
  AddtoPack(data, beta_);

  // Pack matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // Extract kappa and beta
  ExtractfromPack(position, data, kappa_);
  ExtractfromPack(position, data, beta_);

  // Extract matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Maxwell_0d_acinus_Ogden*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | Setup routine to add Ogden material specific parameters kappa and    |
 | beta to material                                         roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::Setup(INPUT::LineDefinition* linedef)
{
  linedef->ExtractDouble("KAPPA", kappa_);
  linedef->ExtractDouble("BETA", beta_);
  // TODO bool -variable init, in Evaluate abfragen ob init=true
}


/*----------------------------------------------------------------------*
 | Evaluate Ogden material and build system matrix and rhs. Acinus Type |
 | "VolumetricOgden": continuum mechanics derivation of cauchy stress   |
 | (=hydrostatic pressure) for Ogden material for purely volumetric     |
 | deformation                                                          |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::Evaluate(CORE::LINALG::SerialDenseVector& epnp,
    CORE::LINALG::SerialDenseVector& epn, CORE::LINALG::SerialDenseVector& epnm,
    CORE::LINALG::SerialDenseMatrix& sysmat, CORE::LINALG::SerialDenseVector& rhs,
    const DRT::REDAIRWAYS::ElemParams& params, const double NumOfAcini, const double Vo,
    double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.putScalar(0.0);
  rhs.putScalar(0.0);

  // Get acinar volume in current timestep
  double acin_vn = params.acin_vn;

  // Get flow in current and next timestep
  double qnp = params.qin_np;
  double qn = params.qin_n;

  // Get acini pressure and beginning and end of acinus element
  double p1n = epn(0);
  double p2n = epn(1);

  // Safety check for NumOfAcini
  if (NumOfAcini < 1.0)
  {
    dserror("Acinus condition at node (%d) has zero acini");
  }
  // Calculate volume and flow per acinuar duct
  double vi_n = (acin_vn / NumOfAcini);
  double qi_n = (qn / NumOfAcini);
  double qi_np = (qnp / NumOfAcini);

  // Linear branches of the Maxwell model (Stiffness2(), B=R_t, B_a=R_a), notation according to
  // interacinar dependency paper
  double Kp_np = Viscosity1() / (Stiffness2() * dt) + 1.0;
  double Kp_n = -Viscosity1() / (Stiffness2() * dt);
  double Kq_np = Viscosity1() * Viscosity2() / (Stiffness2() * dt) + Viscosity1() + Viscosity2();
  double Kq_n = -Viscosity1() * Viscosity2() / (Stiffness2() * dt);
  double rhsLin = -Kp_n * (p1n - p2n) + Kq_n * qi_n;

  // Branch E_1 of the Maxwell model: Hydrostatic pressure (=Cauchy stress) for Ogden material
  // P_1  = P_c + P_d
  // where P_c = (kappa/beta_)*(lambda^(-3))
  //       P_d =-(kappa/beta_)*(lambda^(-3-3*beta_))
  //       \lambda is the volumetric strain ratio, \lambda = (V/Vo)^(1/3)
  double vi_np = qi_np * dt + vi_n;
  double Kq_npNL = (Viscosity1() / Stiffness2()) *
                   (-kappa_ * Vo / (pow(vi_np, 2.0) * beta_) +
                       (beta_ + 1.0) * kappa_ * (pow(Vo / (vi_np), beta_ + 1.0)) / ((vi_np)*beta_));
  double rhsNL = (Vo / vi_n) * (kappa_ / beta_) * (1 - pow((Vo / vi_n), beta_));

  // add linearisation part to system matrix
  Kq_np += Kq_npNL;

  // Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
  sysmat(0, 0) = -1.0 * (Kp_np / Kq_np) * NumOfAcini;
  sysmat(0, 1) = 1.0 * (Kp_np / Kq_np) * NumOfAcini;
  sysmat(1, 0) = 1.0 * (Kp_np / Kq_np) * NumOfAcini;
  sysmat(1, 1) = -1.0 * (Kp_np / Kq_np) * NumOfAcini;

  // Build the corresponding right hand side
  rhs(0) = -1.0 * ((rhsLin + rhsNL) * NumOfAcini / Kq_np);
  rhs(1) = 1.0 * ((rhsLin + rhsNL) * NumOfAcini / Kq_np);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Maxwell_0d_acinus_Ogden::GetParams(std::string parametername)
{
  if (parametername == "kappa")
    return kappa_;
  else if (parametername == "beta")
    return beta_;
  else
  {
    dserror("Chosen Parameter can not be returned with this function!");
    return 0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::SetParams(std::string parametername, double new_value)
{
  if (parametername == "kappa")
    kappa_ = new_value;
  else
    dserror("Chosen Parameter can not be set with this function yet!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus_Ogden::VisNames(std::map<std::string, int>& names)
{
  std::string fiber = "kappa";
  names[fiber] = 1;  // scalar
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::Maxwell_0d_acinus_Ogden::VisData(
    const std::string& name, std::vector<double>& data, int eleGID)
{
  if (name == "kappa")
  {
    if ((int)data.size() != 1) dserror("size mismatch");

    data[0] = GetParams("kappa");
  }
  else
  {
    return false;
  }
  return true;
}

BACI_NAMESPACE_CLOSE
