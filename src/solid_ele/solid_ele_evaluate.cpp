/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the solid element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "structure_new_elements_paramsinterface.H"
#include "solid_ele_neumann_evaluator.H"
#include "structure_new_elements_paramsinterface.H"
#include "lib_dserror.H"
#include "solid_ele.H"
#include "solid_ele_factory.H"
#include "solid_ele_calc_interface.H"

namespace
{
  void LumpMatrix(Epetra_SerialDenseMatrix& matrix)
  {
    dsassert(matrix.N() == matrix.M(), "The provided mass matrix is not a square matrix!");

    // we assume mass is a square matrix
    for (int c = 0; c < matrix.N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < matrix.M(); ++r)  // parse rows
      {
        d += matrix(r, c);  // accumulate row entries
        matrix(r, c) = 0.0;
      }
      matrix(c, c) = d;  // apply sum of row entries on diagonal
    }
  }

  inline std::vector<char>& GetMutableStressData(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().MutableStressDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  inline std::vector<char>& GetMutableStrainData(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().MutableStrainDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  inline INPAR::STR::StressType GetIOStressType(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStressOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress");
    }
  }

  inline INPAR::STR::StrainType GetIOStrainType(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStrainOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain");
    }
  }
}  // namespace

int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->MaterialPostSetup(*this, *SolidMaterial());
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  SetParamsInterfacePtr(params);

  const ELEMENTS::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return ParamsInterface().GetActionType();
        else
          return String2ActionType(params.get<std::string>("action", "none"));
      });

  switch (action)
  {
    case DRT::ELEMENTS::struct_calc_nlnstiff:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
      return 0;
    }
    case struct_calc_internalforce:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, nullptr, nullptr);
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      LumpMatrix(elemat2);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Update(
          *this, *SolidMaterial(), discretization, lm, params);
      return 0;
    case DRT::ELEMENTS::struct_calc_recover:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Recover(*this, discretization, lm, params);
      return 0;
    }
    case struct_calc_stress:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->CalculateStress(*this, *SolidMaterial(),
          StressIO{GetIOStressType(*this, params), GetMutableStressData(*this, params)},
          StrainIO{GetIOStrainType(*this, params), GetMutableStrainData(*this, params)},
          discretization, lm, params);
      return 0;
    }
    case struct_init_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->InitializeGaussPointDataOutput(
          *this, *SolidMaterial(), *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case struct_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateGaussPointDataOutput(
          *this, *SolidMaterial(), *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case ELEMENTS::struct_calc_reset_all:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->ResetAll(*this, *SolidMaterial());
      return 0;
    case ELEMENTS::struct_calc_reset_istep:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->ResetToLastConverged(*this, *SolidMaterial());
      return 0;
    case DRT::ELEMENTS::struct_calc_predict:
      // do nothing for now
      return 0;
    default:
      dserror("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
  }

  return 0;
}
int DRT::ELEMENTS::Solid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return ParamsInterface().GetTotalTime();
        else
          return params.get("total time", -1.0);
      });

  DRT::ELEMENTS::EvaluateNeumannByElement(*this, discretization, condition, lm, elevec1, time);
  return 0;
}
