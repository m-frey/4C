/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_bele3.hpp"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_utils.hpp"
#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_function.hpp"

BACI_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Bele3Line::ActionType act = Bele3Line::none;
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "integrate_Shapefunction")
    act = Bele3Line::integrate_Shapefunction;
  else
    dserror("Unknown type of action for Bele3Line");

  switch (act)
  {
    case integrate_Shapefunction:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      //      if (parent_->IsMoving())
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp != Teuchos::null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
        }
        else
        {
          std::cout << "could not get displacement vector to compute current positions!"
                    << std::endl;
        }
      }

      IntegrateShapeFunction(params, discretization, lm, elevec1, mydispnp);
      break;
    }
    default:
      dserror("Unknown type of action for Bele3Line");
      break;
  }  // end of switch(act)

  return 0;

}  // DRT::ELEMENTS::Bele3Line::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  const double thsl = params.get("thsl", 0.0);

  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  const std::vector<int>* functions = condition.Get<std::vector<int>>("funct");

  // set number of nodes
  const size_t iel = this->NumNode();

  const CORE::FE::CellType distype = this->Shape();

  // gaussian points
  const CORE::FE::GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const CORE::FE::IntegrationPoints1D intpoints(gaussrule);


  // allocate vector for shape functions and for derivatives
  CORE::LINALG::SerialDenseVector funct(iel);
  CORE::LINALG::SerialDenseMatrix deriv(1, iel);

  // node coordinates
  CORE::LINALG::SerialDenseMatrix xye(2, iel);

  // get node coordinates
  for (size_t i = 0; i < iel; ++i)
  {
    xye(0, i) = this->Nodes()[i]->X()[0];
    xye(1, i) = this->Nodes()[i]->X()[1];
  }


  // loop over integration points
  for (int gpid = 0; gpid < intpoints.nquad; gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];
    // get shape functions and derivatives in the line
    CORE::FE::shape_function_1D(funct, e1, distype);
    CORE::FE::shape_function_1D_deriv1(deriv, e1, distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye, deriv, iel);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)

    const double fac = intpoints.qwgt[gpid] * dr * thsl;

    // factor given by spatial function
    double functionfac = 1.0;
    // determine coordinates of current Gauss point
    double coordgp[3];
    coordgp[0] = 0.0;
    coordgp[1] = 0.0;
    coordgp[2] = 0.0;
    for (size_t i = 0; i < iel; i++)
    {
      coordgp[0] += xye(0, i) * funct[i];
      coordgp[1] += xye(1, i) * funct[i];
    }

    int functnum = -1;
    const double* coordgpref = coordgp;  // needed for function evaluation

    for (size_t node = 0; node < iel; ++node)
    {
      for (size_t dim = 0; dim < 3; dim++)
      {
        // factor given by spatial function
        if (functions) functnum = (*functions)[dim];
        {
          if (functnum > 0)
            // evaluate function at current gauss point (3D position vector required!)
            functionfac = GLOBAL::Problem::Instance()
                              ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                              .Evaluate(coordgpref, time, dim);
          else
            functionfac = 1.0;
        }

        elevec1[node * numdf + dim] +=
            funct[node] * (*onoff)[dim] * (*val)[dim] * fac * functionfac;
      }
    }
  }  // end of loop over integrationen points


  // dserror("Line Neumann condition not yet implemented for Bele3");
  return 0;
}

CORE::FE::GaussRule1D DRT::ELEMENTS::Bele3Line::getOptimalGaussrule(
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
      break;
  }
  return rule;
}


double DRT::ELEMENTS::Bele3Line::f2_substitution(const CORE::LINALG::SerialDenseMatrix xye,
    const CORE::LINALG::SerialDenseMatrix deriv, const int iel)
{
  // compute derivative of parametrization
  double dr = 0.0;
  CORE::LINALG::SerialDenseVector der_par(iel);
  CORE::LINALG::multiplyNT(der_par, xye, deriv);
  dr = CORE::LINALG::Norm2(der_par);
  return dr;
}

/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over line (public)              g.bau 07/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Line::IntegrateShapeFunction(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, const std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, const std::vector<double>& edispnp)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  //  const double thsl = params.get("thsl",1.0);

  /*
    // find out whether we will use a time curve
    bool usetime = true;
    const double time = params.get("total time",-1.0);
    if (time<0.0) usetime = false;
  */

  // set number of nodes
  const size_t iel = this->NumNode();

  // gaussian points
  const CORE::FE::CellType distype = this->Shape();
  const CORE::FE::GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const CORE::FE::IntegrationPoints1D intpoints(gaussrule);

  // allocate vector for shape functions and for derivatives
  CORE::LINALG::SerialDenseVector funct(iel);
  CORE::LINALG::SerialDenseMatrix deriv(1, iel);

  // node coordinates
  CORE::LINALG::SerialDenseMatrix xye(2, iel);

  // get node coordinates
  for (size_t i = 0; i < iel; ++i)
  {
    xye(0, i) = this->Nodes()[i]->X()[0];
    xye(1, i) = this->Nodes()[i]->X()[1];
  }

  //  if (parent_->IsMoving())
  {
    dsassert(edispnp.size() != 0, "paranoid");

    for (size_t i = 0; i < iel; i++)
    {
      xye(0, i) += edispnp[3 * i];
      xye(1, i) += edispnp[3 * i + 1];
    }
  }

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.nquad; gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];
    // get shape functions and derivatives in the line
    CORE::FE::shape_function_1D(funct, e1, distype);
    CORE::FE::shape_function_1D_deriv1(deriv, e1, distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye, deriv, iel);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)

    // double fac = intpoints.qwgt[gpid] *dr * thsl;
    const double fac = intpoints.qwgt[gpid] * dr;

    for (size_t node = 0; node < iel; ++node)
    {
      for (size_t dim = 0; dim < 3; dim++)
      {
        elevec1[node * numdf + dim] += funct[node] * fac;
      }
    }
  }  // end of loop over integrationen points

  return;
}  // DRT::ELEMENTS::Bele3Line::IntegrateShapeFunction

BACI_NAMESPACE_CLOSE
