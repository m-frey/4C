/*! \file

\brief Evaluation of neumann loads for shell7p elements

\level 3
*/

#include "4C_shell7p_ele_neumann_evaluator.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_shell7p_ele_calc_lib.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

void Discret::ELEMENTS::Shell::EvaluateNeumannByElement(Core::Elements::Element& ele,
    const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, Core::LinAlg::SerialDenseVector& element_force_vector,
    Core::LinAlg::SerialDenseMatrix* element_stiffness_matrix, double total_time)
{
  switch (ele.Shape())
  {
    case Core::FE::CellType::quad4:
      return evaluate_neumann<Core::FE::CellType::quad4>(ele, discretization, condition,
          dof_index_array, element_force_vector, element_stiffness_matrix, total_time);
    case Core::FE::CellType::quad8:
      return evaluate_neumann<Core::FE::CellType::quad8>(ele, discretization, condition,
          dof_index_array, element_force_vector, element_stiffness_matrix, total_time);
    case Core::FE::CellType::quad9:
      return evaluate_neumann<Core::FE::CellType::quad9>(ele, discretization, condition,
          dof_index_array, element_force_vector, element_stiffness_matrix, total_time);
    case Core::FE::CellType::tri3:
      return evaluate_neumann<Core::FE::CellType::tri3>(ele, discretization, condition,
          dof_index_array, element_force_vector, element_stiffness_matrix, total_time);
    case Core::FE::CellType::tri6:
      return evaluate_neumann<Core::FE::CellType::tri6>(ele, discretization, condition,
          dof_index_array, element_force_vector, element_stiffness_matrix, total_time);
    default:
      FOUR_C_THROW(
          "The discretization type you are trying to evaluate the Neumann condition for is not yet "
          "implemented in EvaluateNeumannByElement.");
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Shell::evaluate_neumann(Core::Elements::Element& ele,
    const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, Core::LinAlg::SerialDenseVector& element_force_vector,
    Core::LinAlg::SerialDenseMatrix* element_stiffness_matrix, double total_time)
{
  constexpr auto num_dim = Shell::DETAIL::num_dim;
  constexpr auto numnod = Shell::DETAIL::num_node<distype>;
  constexpr auto noddof = Shell::DETAIL::node_dof;

  Core::FE::IntegrationPoints2D intpoints =
      create_gauss_integration_points<distype>(get_gauss_rule<distype>());

  // IMPORTANT: The 'neum_orthopressure' case represents a truly nonlinear follower-load
  // acting on the spatial configuration. Therefore, it needs to be linearized. On the
  // contrary, the simplified 'neum_pseudo_orthopressure' option allows for an approximative
  // modeling of an orthopressure load without the need to do any linearization. However,
  // this can only be achieved by referring the 'neum_pseudo_orthopressure' load to the last
  // converged configuration, which introduces an error as compared with 'neum_orthopressure'.
  bool loadlin = (element_stiffness_matrix != nullptr);

  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live,
    neum_live_FSI,
    neum_orthopressure,
    neum_pseudo_orthopressure,
    neum_consthydro_z,
    neum_increhydro_z,
    neum_opres_FSI
  };
  LoadType ltype = neum_none;

  // type of configurations
  enum Configuration
  {
    config_none,
    config_material,       // material configuration
    config_lastconverged,  // last converged configuration
    config_spatial         // spatial configuration
  };

  Configuration config = config_none;

  const std::string& type = condition.parameters().get<std::string>("type");
  if (type == "neum_live")
  {
    ltype = neum_live;
    config = config_material;
  }
  else if (type == "neum_live_FSI")
  {
    ltype = neum_live_FSI;
    config = config_material;
  }
  else if (type == "neum_pseudo_orthopressure")
  {
    ltype = neum_pseudo_orthopressure;
    config = config_lastconverged;
  }
  else if (type == "neum_orthopressure")
  {
    ltype = neum_orthopressure;
    config = config_spatial;
  }
  else if (type == "neum_consthydro_z")
  {
    ltype = neum_consthydro_z;
    config = config_lastconverged;
  }
  else if (type == "neum_increhydro_z")
  {
    ltype = neum_increhydro_z;
    config = config_lastconverged;
  }
  else
  {
    FOUR_C_THROW("Unknown type of SurfaceNeumann condition");
  }
  // get values and switches from the condition
  const auto& onoff = condition.parameters().get<std::vector<int>>("onoff");
  const auto& value = condition.parameters().get<std::vector<double>>("val");

  // ensure that at least as many curves/functs as dofs are available
  if (onoff.size() < DETAIL::node_dof)
    FOUR_C_THROW(
        "Fewer functions or curves defined than the element's nodal degree of freedoms (6).");

  for (std::size_t checkdof = num_dim; checkdof < onoff.size(); ++checkdof)
  {
    if (onoff[checkdof] != 0)
    {
      FOUR_C_THROW(
          "You have activated more than %d dofs in your Neumann boundary condition. This is higher "
          "than the dimension of the element.",
          num_dim);
    }
  }
  // get ids of functions of space and time
  const auto* function_ids = &condition.parameters().get<std::vector<int>>("funct");

  // integration loops
  std::array<double, 2> xi_gp;
  for (int gp = 0; gp < intpoints.NumPoints(); ++gp)
  {
    // get gauss points from integration rule
    xi_gp[0] = intpoints.qxg[gp][0];
    xi_gp[1] = intpoints.qxg[gp][1];
    // get gauss weight at current gp
    double gpweights = intpoints.qwgt[gp];
    // get shape functions and derivatives at gaussian points
    ShapefunctionsAndDerivatives<distype> shapefunctions =
        EvaluateShapefunctionsAndDerivs<distype>(xi_gp);
    Core::LinAlg::Matrix<num_dim, num_dim> g_metrics_kovariant(true);

    Core::LinAlg::Matrix<numnod, num_dim> x_refe(true);
    Core::LinAlg::Matrix<numnod, num_dim> x_curr(true);

    for (int i = 0; i < numnod; ++i)
    {
      x_refe(i, 0) = ele.Nodes()[i]->X()[0];
      x_refe(i, 1) = ele.Nodes()[i]->X()[1];
      x_refe(i, 2) = ele.Nodes()[i]->X()[2];
    }

    switch (config)
    {
      case config_material:
      {
        // no linearization needed for load in material configuration
        loadlin = false;
        g_metrics_kovariant.MultiplyNN(1.0, shapefunctions.derivatives_, x_refe, 0.0);
      }
      break;
      case config_lastconverged:
      {
        // no linearization needed for load in last converged configuration
        loadlin = false;
        const Epetra_Vector& disp = *discretization.GetState("displacement");
        std::vector<double> displacements(dof_index_array.size());
        Core::FE::ExtractMyValues(disp, displacements, dof_index_array);

        spatial_configuration<distype>(x_curr, x_refe, displacements, 0);

        g_metrics_kovariant.MultiplyNN(1.0, shapefunctions.derivatives_, x_curr, 0.0);
      }
      break;
      case config_spatial:
      {
        const Epetra_Vector& disp = *discretization.GetState("displacement new");
        std::vector<double> displacements(dof_index_array.size());
        Core::FE::ExtractMyValues(disp, displacements, dof_index_array);

        spatial_configuration<distype>(x_curr, x_refe, displacements, 0);

        g_metrics_kovariant.MultiplyNN(1.0, shapefunctions.derivatives_, x_curr, 0.0);
      }
      break;
      default:
        FOUR_C_THROW("Unknown type of configuration");
    }

    // get thickness direction derivative perpendicular to g1 and g2 (with area as length)
    // -> g3 = (g1 x g2) / (|g1 x g2 |)
    Core::LinAlg::Matrix<num_dim, 1> g3(true);
    g3(0) = g_metrics_kovariant(0, 1) * g_metrics_kovariant(1, 2) -
            g_metrics_kovariant(0, 2) * g_metrics_kovariant(1, 1);
    g3(1) = g_metrics_kovariant(0, 2) * g_metrics_kovariant(1, 0) -
            g_metrics_kovariant(0, 0) * g_metrics_kovariant(1, 2);
    g3(2) = g_metrics_kovariant(0, 0) * g_metrics_kovariant(1, 1) -
            g_metrics_kovariant(0, 1) * g_metrics_kovariant(1, 0);
    // compute line increment ds
    double ds = g3.Norm2();
    if (ds <= 1.0e-14) FOUR_C_THROW("Element Area equal 0.0 or negativ detected");

    // material/reference coordinates of Gaussian point
    Core::LinAlg::Matrix<num_dim, 1> gauss_point_reference_coordinates;
    gauss_point_reference_coordinates.MultiplyTN(x_refe, shapefunctions.shapefunctions_);
    // current coordinates of Gaussian point
    Core::LinAlg::Matrix<num_dim, 1> gauss_point_current_coordinates;
    if (ltype != neum_live)
      gauss_point_current_coordinates.MultiplyTN(x_curr, shapefunctions.shapefunctions_);

    // evaluate the spatial function in the reference state
    std::vector<double> function_scale_factors(3, 1.0);
    for (auto dim = 0; dim < num_dim; dim++)
    {
      if (onoff[dim])
      {
        // function evaluation
        const int function_number = (function_ids != nullptr) ? (*function_ids)[dim] : -1;
        function_scale_factors[dim] =
            (function_number > 0)
                ? Global::Problem::Instance()
                      ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(function_number - 1)
                      .evaluate(gauss_point_reference_coordinates.A(), total_time, dim)
                : 1.0;
      }
    }

    double value_times_integration_factor = 0.0;
    for (int dim = 0; dim < Shell::DETAIL::num_dim; ++dim)
    {
      switch (ltype)
      {
        case neum_live:  // uniform load on reference configuration
        case neum_live_FSI:
        {
          value_times_integration_factor =
              value[dim] * function_scale_factors[dim] * gpweights * ds;
        }
        break;
        // hydrostatic pressure dependent on z-coordinate of gaussian point
        case neum_consthydro_z:
        {
          if (onoff[2] != 1) FOUR_C_THROW("hydropressure must be on third dof");
          value_times_integration_factor =
              -value[2] * function_scale_factors[2] * gpweights * g3(dim, 0);
        }
        break;
        // hydrostatic pressure dependent on z-coord of gaussian point increasing with time in
        // height
        case neum_increhydro_z:
        {
          if (onoff[2] != 1) FOUR_C_THROW("hydropressure must be on third dof");
          double height = total_time * 10.0;
          double pressure = 0.0;
          if (gauss_point_current_coordinates(2, 0) <= height)
            pressure = -value[2] * (height - gauss_point_current_coordinates(2, 0));
          value_times_integration_factor = pressure * gpweights * g3(dim, 0);
        }
        break;
        // pressure orthogonal to surface
        case neum_pseudo_orthopressure:
        case neum_orthopressure:
        case neum_opres_FSI:
        {
          if (onoff[2] != 1) FOUR_C_THROW("orthopressure must be on third dof");
          value_times_integration_factor =
              value[2] * function_scale_factors[2] * gpweights * g3(dim, 0);
        }
        break;
        default:
          FOUR_C_THROW("Unknown type of SurfaceNeumann load");
      }  // end switch(eletype)

      for (auto nodeid = 0; nodeid < Shell::DETAIL::num_node<distype>; ++nodeid)
      {
        // Evaluates the Neumann boundary condition: f_{x,y,z}^i=\sum_j N^i(xi^j) * value(t) *
        // integration_factor_j
        // assembles the element force vector [f_x^1, f_y^1, f_z^1, ..., f_x^n, f_y^n, f_z^n]
        element_force_vector[nodeid * noddof + dim] +=
            shapefunctions.shapefunctions_(nodeid) * value_times_integration_factor;
      }
    }
    // load linerization (if neccessary)
    if (loadlin)
    {
      constexpr auto numdof = DETAIL::numdofperelement<distype>;
      Core::LinAlg::Matrix<num_dim, numdof> dnormal;
      Core::LinAlg::Matrix<num_dim, numdof> dg1;
      Core::LinAlg::Matrix<num_dim, numdof> dg2;

      // linearization of basis vectors
      for (int nodeid = 0; nodeid < numnod; ++nodeid)
      {
        for (int dim = 0; dim < num_dim; ++dim)
        {
          dg1(dim, nodeid + num_dim + dim) = shapefunctions.derivatives_(0, nodeid);
          dg2(dim, nodeid + num_dim + dim) = shapefunctions.derivatives_(1, nodeid);
        }
      }
      // linearization of local surface normal vector
      for (int dim = 0; dim < num_dim; ++dim)
      {
        dnormal(0, dim) =
            dg1(1, dim) * g_metrics_kovariant(1, 2) + g_metrics_kovariant(0, 1) * dg2(2, dim) -
            dg1(2, dim) * g_metrics_kovariant(1, 1) - g_metrics_kovariant(0, 2) * dg2(1, dim);
        dnormal(1, dim) =
            dg1(2, dim) * g_metrics_kovariant(1, 0) + g_metrics_kovariant(0, 2) * dg2(0, dim) -
            dg1(0, dim) * g_metrics_kovariant(1, 2) - g_metrics_kovariant(0, 0) * dg2(2, dim);
        dnormal(2, dim) =
            dg1(0, dim) * g_metrics_kovariant(1, 1) + g_metrics_kovariant(0, 0) * dg2(1, dim) -
            dg1(1, dim) * g_metrics_kovariant(1, 0) - g_metrics_kovariant(0, 1) * dg2(0, dim);
      }
      // build surface element load linearization matrix
      // (CAREFUL: Minus sign due to the fact that external forces enter the global
      // residual vector with a minus sign, too! However, the load linaerization is
      // simply added to the global tangent stiffness matrix, thus we explicitly
      // need to set the minus sign here.)
      for (int nodeid = 0; nodeid < numnod; ++nodeid)
        for (int dim = 0; dim < num_dim; dim++)
          for (int dof = 0; dof < element_force_vector.numRows(); dof++)
            (*element_stiffness_matrix)(nodeid * noddof + dim, dof) -=
                shapefunctions.shapefunctions_(nodeid) * value[2] * function_scale_factors[2] *
                gpweights * dnormal(dim, dof);
    }
  }
}
FOUR_C_NAMESPACE_CLOSE
