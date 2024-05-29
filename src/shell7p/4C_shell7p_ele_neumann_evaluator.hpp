/*! \file

\brief Evaluation of neumann loads for shell7p element

\level 3
 */

#ifndef FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP
#define FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

}  // namespace DRT


namespace DRT::ELEMENTS::SHELL
{
  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element.
   *
   * @note This function determines the shape of the element at runtime and calls the respective
   * templated version of @p evaluate_neumann. If you already know the CORE::FE::CellType
   * of the element at compile-time, you could directly call @evaluate_neumann.
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  void EvaluateNeumannByElement(CORE::Elements::Element& element,
      const DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
      const std::vector<int>& dof_index_array,
      CORE::LINALG::SerialDenseVector& element_force_vector,
      CORE::LINALG::SerialDenseMatrix* element_stiffness_matrix, double total_time);

  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element with the
   * discretization type known at compile time.
   *
   *
   * @tparam distype discretization type known at compile time
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  template <CORE::FE::CellType distype>
  void evaluate_neumann(CORE::Elements::Element& element, const DRT::Discretization& discretization,
      CORE::Conditions::Condition& condition, const std::vector<int>& dof_index_array,
      CORE::LINALG::SerialDenseVector& element_force_vector,
      CORE::LINALG::SerialDenseMatrix* element_stiffness_matrix, double total_time);

}  // namespace DRT::ELEMENTS::SHELL

FOUR_C_NAMESPACE_CLOSE

#endif
