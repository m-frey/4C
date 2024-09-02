/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional acinus elements with non-linear
spring with double-exponential behaviour, inherits from Maxwell_0d_acinus

The originally linear spring (Stiffness1) of the 4-element Maxwell model is substituted by a
double-exponential pressure-volume relation (derivation: see Ismail Mahmoud's dissertation,
chapter 3.4)

Necessary input lines:
(material section)
MAT 3 MAT_0D_MAXWELL_ACINUS_OGDEN Stiffness1 1.0 Stiffness2 5249.1 Viscosity1 3221.86 Viscosity2
1000.0 // acinus properties;
(element section)
1 RED_ACINUS  MAT 3 TYPE DoubleExponential AcinusVolume 300 AlveolarDuctVolume 0.03711 E1_01
0 E1_02 0 E1_EXP1 0 E1_EXP2 0 E1_LIN1 0 E1_LIN2 0 TAU1 0 TAU2 0


\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mat_maxwell_0d_acinus_DoubleExponential.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_red_airways_elem_params.hpp"
#include "4C_red_airways_elementbase.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::Maxwell0dAcinusDoubleExponential::Maxwell0dAcinusDoubleExponential(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Maxwell0dAcinus(matdata)
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::Maxwell0dAcinusDoubleExponential::create_material()
{
  return Teuchos::rcp(new Mat::Maxwell0dAcinusDoubleExponential(this));
}


Mat::Maxwell0dAcinusDoubleExponentialType Mat::Maxwell0dAcinusDoubleExponentialType::instance_;


Core::Communication::ParObject* Mat::Maxwell0dAcinusDoubleExponentialType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::Maxwell0dAcinusDoubleExponential* mxwll_0d_acin =
      new Mat::Maxwell0dAcinusDoubleExponential();
  mxwll_0d_acin->unpack(buffer);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Maxwell0dAcinusDoubleExponential::Maxwell0dAcinusDoubleExponential() : Maxwell0dAcinus() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Maxwell0dAcinusDoubleExponential::Maxwell0dAcinusDoubleExponential(
    Mat::PAR::Maxwell0dAcinus* params)
    : Maxwell0dAcinus(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusDoubleExponential::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // Pack type of this instance of ParObject
  int type = unique_par_object_id();

  add_to_pack(data, type);
  add_to_pack(data, e1_01_);
  add_to_pack(data, e1_lin1_);
  add_to_pack(data, e1_exp1_);
  add_to_pack(data, tau1_);

  add_to_pack(data, e1_02_);
  add_to_pack(data, e1_lin2_);
  add_to_pack(data, e1_exp2_);
  add_to_pack(data, tau2_);

  // Pack matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusDoubleExponential::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // Extract e1_01_, e1_lin1_, e1_exp1_, tau1_
  extract_from_pack(buffer, e1_01_);
  extract_from_pack(buffer, e1_lin1_);
  extract_from_pack(buffer, e1_exp1_);
  extract_from_pack(buffer, tau1_);

  // Extract e1_02_, e1_lin2_, e1_exp2_, tau2_
  extract_from_pack(buffer, e1_02_);
  extract_from_pack(buffer, e1_lin2_);
  extract_from_pack(buffer, e1_exp2_);
  extract_from_pack(buffer, tau2_);

  // Extract matid
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::Maxwell0dAcinusDoubleExponential*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}


/*----------------------------------------------------------------------*
 | Setup routine to add DoubleExponential material specific parameters  |
 | E1_01, E1_LIN1, E1_EXP1, TAU1 and                                    |
 | E1_02, E1_LIN2, E1_EXP2, TAU2 to material                roth 10/2014|
 *----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusDoubleExponential::setup(
    const Core::IO::InputParameterContainer& container)
{
  e1_01_ = container.get<double>("E1_01");
  e1_lin1_ = container.get<double>("E1_LIN1");
  e1_exp1_ = container.get<double>("E1_EXP1");
  tau1_ = container.get<double>("TAU1");

  e1_02_ = container.get<double>("E1_02");
  e1_lin2_ = container.get<double>("E1_LIN2");
  e1_exp2_ = container.get<double>("E1_EXP2");
  tau2_ = container.get<double>("TAU2");
  // TODO bool -variable init, in Evaluate abfragen ob init=true
}


/*----------------------------------------------------------------------*
 | Evaluate DoubleExponential material and build system matrix and rhs. |
 |                                                          roth 10/2014|
 *----------------------------------------------------------------------*/
void Mat::Maxwell0dAcinusDoubleExponential::evaluate(Core::LinAlg::SerialDenseVector& epnp,
    Core::LinAlg::SerialDenseVector& epn, Core::LinAlg::SerialDenseVector& epnm,
    Core::LinAlg::SerialDenseMatrix& sysmat, Core::LinAlg::SerialDenseVector& rhs,
    const Discret::ReducedLung::ElemParams& params, const double NumOfAcini, const double Vo,
    double time, double dt)
{
  // Set sysmat and rhs to zero
  sysmat.putScalar(0.0);
  rhs.putScalar(0.0);

  // Get acinar volume in next and current timestep
  double acin_vnp = params.acin_vnp;
  double acin_vn = params.acin_vn;

  // Get flow in next and current timestep
  double qnp = params.qin_np;
  double qn = params.qin_n;

  // Get acini pressure and beginning and end of acinus element
  double p1n = epn(0);
  double p2n = epn(1);

  // Safety check for NumOfAcini
  if (NumOfAcini < 1.0)
  {
    FOUR_C_THROW("Acinus condition at node (%d) has zero acini");
  }

  // Calculate volume and flow difference per acinuar duct
  double dvnp = (acin_vnp / NumOfAcini) - Vo;
  double dvn = (acin_vn / NumOfAcini) - Vo;

  //------------------------------------------------------------
  // V  = A + B*exp(-K*P)
  //
  // The P-V curve is fitted to create the following
  // P1 = E1.(V-Vo)
  //
  // E1 = e1_01_ + b.(V-Vo) + c.exp(d.(V-Vo))
  //------------------------------------------------------------
  double kp_np = viscosity1() / (stiffness2() * dt) + 1.0;
  double kp_n = -viscosity1() / (stiffness2() * dt);
  double kq_np = viscosity1() * viscosity2() / (stiffness2() * dt) + (viscosity2() + viscosity1());
  double kq_n = -viscosity1() * viscosity2() / (stiffness2() * dt);

  // Get the terms assosciated with the nonlinear behavior of
  double term_nonlin = 0.0;
  double pnpi = 0.0;
  double pnpi2 = 0.0;
  double dpnpi_dt = 0.0;
  double dpnpi2_dt = 0.0;

  // Components of linearized E1
  pnpi = (e1_01_ + e1_lin1_ * dvnp + e1_exp1_ * exp(tau1_ * dvnp)) * dvnp;
  pnpi += (e1_02_ + e1_lin2_ * dvnp + e1_exp2_ * exp(tau2_ * dvnp)) * dvnp;

  pnpi2 = (e1_01_ + 2.0 * e1_lin1_ * dvnp + e1_exp1_ * exp(tau1_ * dvnp) * (tau1_ * dvnp + 1.0));
  pnpi2 += (e1_02_ + 2.0 * e1_lin2_ * dvnp + e1_exp2_ * exp(tau2_ * dvnp) * (tau2_ * dvnp + 1.0));

  // Components of linearized d(E1)/dt
  dpnpi_dt =
      (e1_01_ + 2.0 * e1_lin1_ * dvnp + e1_exp1_ * exp(tau1_ * dvnp) * (1.0 + tau1_ * dvnp)) *
      (dvnp - dvn) / dt;
  dpnpi_dt +=
      (e1_02_ + 2.0 * e1_lin2_ * dvnp + e1_exp2_ * exp(tau2_ * dvnp) * (1.0 + tau2_ * dvnp)) *
      (dvnp - dvn) / dt;

  dpnpi2_dt =
      (2.0 * e1_lin1_ + tau1_ * e1_exp1_ * exp(tau1_ * dvnp) * (1.0 + tau1_ * dvnp) +
          e1_exp1_ * tau1_ * exp(tau1_ * dvnp)) *
          (dvnp - dvn) / dt +
      (e1_01_ + 2.0 * e1_lin1_ * dvnp + e1_exp1_ * exp(tau1_ * dvnp) * (1.0 + tau1_ * dvnp)) / dt;
  dpnpi2_dt +=
      (2.0 * e1_lin2_ + tau2_ * e1_exp2_ * exp(tau2_ * dvnp) * (1.0 + tau2_ * dvnp) +
          e1_exp2_ * tau2_ * exp(tau2_ * dvnp)) *
          (dvnp - dvn) / dt +
      (e1_02_ + 2.0 * e1_lin2_ * dvnp + e1_exp2_ * exp(tau2_ * dvnp) * (1.0 + tau2_ * dvnp)) / dt;

  // Add up the nonlinear terms
  term_nonlin = pnpi + pnpi2 * (-(dvnp) + (qn / NumOfAcini) * dt / 2.0 + dvn);
  kq_np = kq_np + pnpi2 / 2.0 * dt;
  term_nonlin =
      term_nonlin + dpnpi_dt * viscosity1() / stiffness2() +
      dpnpi2_dt * viscosity1() / stiffness2() * (-(dvnp) + (qnp / NumOfAcini) * dt / 2.0 + dvn);
  kq_np = kq_np + dpnpi2_dt * viscosity1() / stiffness2() / 2.0 * dt;

  // Build the system matrix for \boldsymbol{K} * \boldsymbol{P} = \boldsymbol{Q}
  sysmat(0, 0) = -1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(0, 1) = 1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(1, 0) = 1.0 * (kp_np / kq_np) * NumOfAcini;
  sysmat(1, 1) = -1.0 * (kp_np / kq_np) * NumOfAcini;

  // Build the corresponding right hand side
  rhs(0) = -1.0 * ((-kp_n * (p1n - p2n) + term_nonlin) * NumOfAcini / kq_np + (kq_n * qn) / kq_np);
  rhs(1) = 1.0 * ((-kp_n * (p1n - p2n) + term_nonlin) * NumOfAcini / kq_np + (kq_n * qn) / kq_np);
}

FOUR_C_NAMESPACE_CLOSE
