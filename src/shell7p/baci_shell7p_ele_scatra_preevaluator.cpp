/*! \file

\brief PreEvaluator of Shell7p-ScaTra elements

\level 1
*/


#include "baci_shell7p_ele_scatra_preevaluator.H"

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_shell7p_ele_calc_lib.H"
#include "baci_shell7p_ele_scatra.H"


void DRT::ELEMENTS::SHELL::PreEvaluateScatraByElement(DRT::Element& ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& dof_index_array)
{
  switch (ele.Shape())
  {
    case DRT::Element::DiscretizationType::quad4:
      return PreEvaluateScatra<DRT::Element::DiscretizationType::quad4>(
          ele, params, discretization, dof_index_array);
    case DRT::Element::DiscretizationType::quad8:
      return PreEvaluateScatra<DRT::Element::DiscretizationType::quad8>(
          ele, params, discretization, dof_index_array);
    case DRT::Element::DiscretizationType::quad9:
      return PreEvaluateScatra<DRT::Element::DiscretizationType::quad9>(
          ele, params, discretization, dof_index_array);
    case DRT::Element::DiscretizationType::tri3:
      return PreEvaluateScatra<DRT::Element::DiscretizationType::tri3>(
          ele, params, discretization, dof_index_array);
    case DRT::Element::DiscretizationType::tri6:
      return PreEvaluateScatra<DRT::Element::DiscretizationType::tri6>(
          ele, params, discretization, dof_index_array);
    default:
      dserror(
          "The discretization type you are trying to pre-evaluate for shell7p scatra is not yet "
          "implemented.");
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SHELL::PreEvaluateScatra(DRT::Element& ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& dof_index_array)
{
  CORE::DRT::UTILS::IntegrationPoints2D intpoints_midsurface_ =
      CreateGaussIntegrationPoints<distype>(GetGaussRule<distype>());

  if (dof_index_array.Size() > 1)
  {
    // ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1, ele.Nodes()[0]);

    if (dof_index_array[1].Size() != SHELL::DETAIL::num_node<distype> * numscal)
      dserror("location vector length does not match!");

    // name of scalarfield
    std::string scalarfield = "scalarfield";

    if (discretization.HasState(1, scalarfield))
    {
      // get the scalar state
      Teuchos::RCP<const Epetra_Vector> scalarnp = discretization.GetState(1, scalarfield);

      if (scalarnp == Teuchos::null) dserror("can not get state vector %s", scalarfield.c_str());

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myscalar =
          Teuchos::rcp(new std::vector<double>(dof_index_array[1].lm_.size(), 0.0));

      DRT::UTILS::ExtractMyValues(*scalarnp, *myscalar, dof_index_array[1].lm_);

      // element vector for k-th scalar
      std::vector<CORE::LINALG::Matrix<SHELL::DETAIL::num_node<distype>, 1>> elescalar(numscal);
      for (int k = 0; k < numscal; ++k)
      {
        for (int i = 0; i < SHELL::DETAIL::num_node<distype>; ++i)
        {
          (elescalar.at(k))(i, 0) = myscalar->at(numscal * i + k);
        }
      }

      // create vector of gauss point values to be set in params list
      Teuchos::RCP<std::vector<std::vector<double>>> gpscalar =
          Teuchos::rcp(new std::vector<std::vector<double>>(
              intpoints_midsurface_.NumPoints(), std::vector<double>(numscal, 0.0)));

      // allocate vector for shape functions and matrix for derivatives at gp
      CORE::LINALG::Matrix<SHELL::DETAIL::num_node<distype>, 1> shapefunctions(true);

      // loop over gauss points
      for (int gp = 0; gp < intpoints_midsurface_.NumPoints(); ++gp)
      {
        // get gauss points from integration rule
        double xi_gp = intpoints_midsurface_.qxg[gp][0];
        double eta_gp = intpoints_midsurface_.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        CORE::DRT::UTILS::shape_function_2D(shapefunctions, xi_gp, eta_gp, distype);

        // scalar at current gp
        std::vector<double> scalar_curr_gp(numscal, 0.0);

        for (int k = 0; k < numscal; ++k)
        {
          // identical shapefunctions for displacements and scalar fields
          scalar_curr_gp.at(k) = shapefunctions.Dot(elescalar.at(k));
        }

        gpscalar->at(gp) = scalar_curr_gp;
      }

      // set scalar states at gp to params list
      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpscalar);
    }
  }
  std::vector<double> center = {0.0, 0.0};
  auto xrefe = Teuchos::rcp(new std::vector<double>(center));
  params.set<Teuchos::RCP<std::vector<double>>>("position", xrefe);
}
