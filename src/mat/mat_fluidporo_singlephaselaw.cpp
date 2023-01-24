/*----------------------------------------------------------------------*/
/*! \file
 \brief a material defining pressure-saturation relationship
        for a single phase within a multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/


#include "mat_fluidporo_singlephaselaw.H"

#include "mat_par_bundle.H"

#include "lib_globalproblem.H"
#include "lib_get_functionofanything.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLaw::FluidPoroPhaseLaw(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  return;
}


/*----------------------------------------------------------------------*
 *  factory method for phase law                       vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLaw* MAT::PAR::FluidPoroPhaseLaw::CreatePhaseLaw(int phaselawId)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(phaselawId);

  // phase law
  MAT::PAR::FluidPoroPhaseLaw* phaselaw = NULL;

  // build the pressure-saturation law
  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluidporo_phaselaw_linear:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawLinear(curmat));
      phaselaw = static_cast<MAT::PAR::FluidPoroPhaseLawLinear*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_phaselaw_tangent:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawTangent(curmat));
      phaselaw = static_cast<MAT::PAR::FluidPoroPhaseLawTangent*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_phaselaw_constraint:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawConstraint(curmat));
      phaselaw = static_cast<MAT::PAR::FluidPoroPhaseLawConstraint*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_phaselaw_byfunction:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawByFunction(curmat));
      phaselaw = static_cast<MAT::PAR::FluidPoroPhaseLawByFunction*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid pressure-saturation law for material %d", curmat->Type());
      break;
  }

  return phaselaw;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawLinear::FluidPoroPhaseLawLinear(Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroPhaseLaw(matdata),
      numdof_(matdata->GetInt("NUMDOF")),
      presids_(matdata->Get<std::vector<int>>("PRESCOEFF")),
      reltensions_(matdata->GetDouble("RELTENSION")),
      sat0_(matdata->GetDouble("SATURATION_0"))
{
  // check if sizes fit
  if (numdof_ != (int)presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", numdof_, presids_->size());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateSaturation(const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);

  double saturation = sat0_ + reltensions_ * presval;

  return saturation;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  if ((*presids_)[doftoderive] == 0) return 0.0;

  double deriv = reltensions_;

  return deriv * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateSecondDerivOfSaturationWrtPressure(
    int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  // second derivative is zero
  return 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive, double saturation)
{
  if ((*presids_)[doftoderive] == 0) return 0.0;

  double deriv = 1.0 / reltensions_;

  return deriv * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateGenPressure(double saturation)
{
  double presval = 1.0 / reltensions_ * (saturation - sat0_);

  return presval;
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawTangent::FluidPoroPhaseLawTangent(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroPhaseLaw(matdata),
      numdof_(matdata->GetInt("NUMDOF")),
      presids_(matdata->Get<std::vector<int>>("PRESCOEFF")),
      reltensions_(matdata->GetDouble("RELTENSION")),
      exp_(matdata->GetDouble("EXP")),
      sat0_(matdata->GetDouble("SATURATION_0"))
{
  // check if sizes fit
  if (numdof_ != (int)presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", numdof_, presids_->size());
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateSaturation(const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);

  double saturation = sat0_ - std::pow(2 / M_PI * std::atan(reltensions_ * presval), exp_);

  return saturation;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  if ((*presids_)[doftoderive] == 0) return 0.0;

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);

  double deriv = -exp_ * std::pow(2 / M_PI * std::atan(reltensions_ * presval), exp_ - 1.0) * 2.0 *
                 reltensions_ /
                 (M_PI * (1.0 + (reltensions_ * presval) * (reltensions_ * presval)));

  return deriv * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateSecondDerivOfSaturationWrtPressure(
    int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  if ((*presids_)[firstdoftoderive] == 0 || (*presids_)[seconddoftoderive] == 0) return 0.0;

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);

  double secondderiv = 0.0;
  // avoid division by zero in case of small presval --> second deriv goes to zero for presval --> 0
  if (fabs(presval) > 1.0e-12)
  {
    secondderiv = -exp_ * reltensions_ * reltensions_ *
                  (exp_ - 2.0 * reltensions_ * presval * std::atan(reltensions_ * presval) - 1.0) *
                  std::pow(2.0 / M_PI * std::atan(reltensions_ * presval), exp_) /
                  ((1.0 + (reltensions_ * presval) * (reltensions_ * presval)) *
                      (1.0 + (reltensions_ * presval) * (reltensions_ * presval))) /
                  (std::atan(reltensions_ * presval) * std::atan(reltensions_ * presval));
  }

  // second derivative
  return secondderiv * (*presids_)[firstdoftoderive] * (*presids_)[seconddoftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive, double saturation)
{
  if ((*presids_)[doftoderive] == 0) return 0.0;

  double deriv =
      -0.5 * M_PI / (reltensions_ * exp_) * std::pow(sat0_ - saturation, 1.0 / exp_ - 1.0) *
      (1.0 + std::pow(std::tan(0.5 * M_PI * std::pow(sat0_ - saturation, 1.0 / exp_)), 2));

  return deriv * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateGenPressure(double saturation)
{
  double presval =
      1.0 / reltensions_ * std::tan(0.5 * M_PI * std::pow(sat0_ - saturation, 1.0 / exp_));

  return presval;
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawByFunction::FluidPoroPhaseLawByFunction(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroPhaseLaw(matdata),
      numdof_(matdata->GetInt("NUMDOF")),
      presids_(matdata->Get<std::vector<int>>("PRESCOEFF")),
      functionID_saturation_(matdata->GetInt("FUNCTSAT")),
      functionID_pressure_(matdata->GetInt("FUNCTPRES"))
{
  // check if sizes fit
  if (numdof_ != (int)presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", numdof_, presids_->size());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseLawByFunction::Initialize()
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return InitializeInternal<1>();
    case 2:
      return InitializeInternal<2>();
    case 3:
      return InitializeInternal<3>();
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void MAT::PAR::FluidPoroPhaseLawByFunction::InitializeInternal()
{
  if (DRT::UTILS::GetFunctionOfAnything(functionID_saturation_ - 1).NumberComponents() != 1)
    dserror("expected only one component for the saturation evaluation");
  if (DRT::UTILS::GetFunctionOfAnything(functionID_pressure_ - 1).NumberComponents() != 1)
    dserror("expected only one component for the pressure evaluation");


  // initialize pressure vector for function evaluation
  dp_.clear();
  dp_.emplace_back("dp", 0.0);

  // initialize saturation vector for function evaluation
  S_.clear();
  S_.emplace_back("S", 0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateSaturation(
    const std::vector<double>& pressure)
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return EvaluateSaturationInternal<1>(pressure);
    case 2:
      return EvaluateSaturationInternal<2>(pressure);
    case 3:
      return EvaluateSaturationInternal<3>(pressure);
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
      return 0.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateSaturationInternal(
    const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);

  // directly write into entry without checking the name for performance reasons
  dp_[0].second = presval;

  return DRT::UTILS::GetFunctionOfAnything(functionID_saturation_ - 1).Evaluate(dp_, {}, 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive, const std::vector<double>& pressure)
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return EvaluateDerivOfSaturationWrtPressureInternal<1>(doftoderive, pressure);
    case 2:
      return EvaluateDerivOfSaturationWrtPressureInternal<2>(doftoderive, pressure);
    case 3:
      return EvaluateDerivOfSaturationWrtPressureInternal<3>(doftoderive, pressure);
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
      return 0.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfSaturationWrtPressureInternal(
    int doftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  if ((*presids_)[doftoderive] == 0) return 0.0;

  double presval = std::inner_product(presids_->begin(), presids_->end(), pressure.begin(), 0.0);
  // directly write into entry without checking the name for performance reasons
  dp_[0].second = presval;

  std::vector<double> deriv =
      DRT::UTILS::GetFunctionOfAnything(functionID_saturation_ - 1).EvaluateDerivative(dp_, {}, 0);

  return deriv[0] * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateSecondDerivOfSaturationWrtPressure(
    int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure)
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(),
        presids_->size());

  // TODO: implementation for phaselaw by function --> really necessary???
  return 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive, double saturation)
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return EvaluateDerivOfPressureWrtSaturationInternal<1>(doftoderive, saturation);
    case 2:
      return EvaluateDerivOfPressureWrtSaturationInternal<2>(doftoderive, saturation);
    case 3:
      return EvaluateDerivOfPressureWrtSaturationInternal<3>(doftoderive, saturation);
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
      return 0.0;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfPressureWrtSaturationInternal(
    int doftoderive, double saturation)
{
  if ((*presids_)[doftoderive] == 0) return 0.0;

  // directly write into entry without checking the name for performance reasons
  S_[0].second = saturation;

  std::vector<double> deriv =
      DRT::UTILS::GetFunctionOfAnything(functionID_pressure_ - 1).EvaluateDerivative(S_, {}, 0);

  return deriv[0] * (*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateGenPressure(double saturation)
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return EvaluateGenPressureInternal<1>(saturation);
    case 2:
      return EvaluateGenPressureInternal<2>(saturation);
    case 3:
      return EvaluateGenPressureInternal<3>(saturation);
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
      return 0.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateGenPressureInternal(double saturation)
{
  // directly write into entry without checking the name for performance reasons
  S_[0].second = saturation;

  return DRT::UTILS::GetFunctionOfAnything(functionID_pressure_ - 1).Evaluate(S_, {}, 0);
}
