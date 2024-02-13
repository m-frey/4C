/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_vele3.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  return 0;
}  // DRT::ELEMENTS::Vele3Line::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}



/*----------------------------------------------------------------------*
 |  Optimal Gauss rule (public)                              u.may 04/09|
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule1D DRT::ELEMENTS::Vele3Line::getOptimalGaussrule(
    const CORE::FE::CellType& distype)
{
  CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::line2:
      rule = CORE::FE::GaussRule1D::line_2point;
      break;
    case CORE::FE::CellType::line3:
      rule = CORE::FE::GaussRule1D::line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

BACI_NAMESPACE_CLOSE
