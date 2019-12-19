/*----------------------------------------------------------------------*/
/*! \file

\brief volume element

\maintainer Jonas Eichinger

\level 2
*/
/*----------------------------------------------------------------------*/

#include "vele3.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      u.may 05/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}
