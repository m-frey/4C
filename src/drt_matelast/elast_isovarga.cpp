/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for Varga's material
The input line should read
  MAT 1 ELAST_IsoVarga MUE 200.0

\level 2

*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isovarga.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVarga::IsoVarga(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), mue_(matdata->GetDouble("MUE")), beta_(matdata->GetDouble("BETA"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVarga::IsoVarga(MAT::ELASTIC::PAR::IsoVarga* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVarga::AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
    double& shearmod                                          ///< variable to add upon
    ) const
{
  // indeed, a shear modulus is provided
  haveshearmod = haveshearmod or true;

  // material parameters for isochoric part
  shearmod += params_->mue_;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVarga::AddCoefficientsStretchesModified(LINALG::Matrix<3, 1>& modgamma,
    LINALG::Matrix<6, 1>& moddelta, const LINALG::Matrix<3, 1>& modstr)
{
  // parameters
  const double alpha = 2.0 * params_->mue_ - params_->beta_;
  const double beta = params_->beta_;

  // first derivatives
  // \frac{\partial Psi}{\partial \bar{\lambda}_1}
  modgamma(0) += alpha - beta / (modstr(0) * modstr(0));
  // \frac{\partial Psi}{\partial \bar{\lambda}_2}
  modgamma(1) += alpha - beta / (modstr(1) * modstr(1));
  // \frac{\partial Psi}{\partial \bar{\lambda}_3}
  modgamma(2) += alpha - beta / (modstr(2) * modstr(2));

  // second derivatives
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_1^2}
  moddelta(0) += 2.0 * beta / (modstr(0) * modstr(0) * modstr(0));
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_2^2}
  moddelta(1) += 2.0 * beta / (modstr(1) * modstr(1) * modstr(1));
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_3^2}
  moddelta(2) += 2.0 * beta / (modstr(2) * modstr(2) * modstr(2));
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_1 \partial\bar{\lambda}_2}
  moddelta(3) += 0.0;
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_2 \partial\bar{\lambda}_3}
  moddelta(4) += 0.0;
  // \frac{\partial^2 Psi}{\partial\bar{\lambda}_3 \partial\bar{\lambda}_1}
  moddelta(5) += 0.0;

  // done
  return;
}

/*----------------------------------------------------------------------*/
