/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a logarithmic neo-Hooke material according to Bonet and Wood, "Nonlinear
continuum mechanics for finite element analysis", Cambridge, 1997.

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_couplogneohooke.H"

#include "baci_mat_par_material.H"


MAT::ELASTIC::PAR::CoupLogNeoHooke::CoupLogNeoHooke(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata)
{
  std::string parmode = *(matdata->Get<std::string>("MODE"));
  double c1 = matdata->GetDouble("C1");
  double c2 = matdata->GetDouble("C2");

  if (parmode == "YN")
  {
    if (c2 <= 0.5 and c2 > -1.0)
    {
      lambda_ = (c2 == 0.5) ? 0.0 : c1 * c2 / ((1.0 + c2) * (1.0 - 2.0 * c2));
      mue_ = c1 / (2.0 * (1.0 + c2));  // shear modulus
    }
    else
      dserror("Poisson's ratio must be between -1.0 and 0.5!");
  }
  else if (parmode == "Lame")
  {
    mue_ = c1;
    lambda_ = c2;
  }
  else
    dserror(
        "unknown parameter set for NeoHooke material!\n Must be either YN (Young's modulus and "
        "Poisson's ratio) or Lame");
}

MAT::ELASTIC::CoupLogNeoHooke::CoupLogNeoHooke(MAT::ELASTIC::PAR::CoupLogNeoHooke* params)
    : params_(params)
{
}

void MAT::ELASTIC::CoupLogNeoHooke::AddShearMod(
    bool& haveshearmod,  ///< non-zero shear modulus was added
    double& shearmod     ///< variable to add upon
) const
{
  haveshearmod = true;

  shearmod += params_->mue_;
}

void MAT::ELASTIC::CoupLogNeoHooke::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double mue = params_->mue_;
  const double lambda = params_->lambda_;

  // strain energy: Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
  //                     - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
  //                     + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
  // add to overall strain energy

  psi += mue * 0.5 * (prinv(0) - 3.) - mue * log(sqrt(prinv(2))) +
         lambda * 0.5 * pow(log(sqrt(prinv(2))), 2.);
}

void MAT::ELASTIC::CoupLogNeoHooke::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  // ln of determinant of deformation gradient
  const double logdetf = std::log(std::sqrt(prinv(2)));
  const double mue = params_->mue_;
  const double lambda = params_->lambda_;

  dPI(0) += mue * 0.5;
  dPI(2) += (lambda * logdetf) / (2. * prinv(2)) - mue / (2. * prinv(2));

  ddPII(2) += lambda / (4. * prinv(2) * prinv(2)) + mue / (2. * prinv(2) * prinv(2)) -
              (lambda * logdetf) / (2. * prinv(2) * prinv(2));
}