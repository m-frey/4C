/*! \file

\brief Definition of routines for calculation of solid poro element with pressure based
implementation

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_PRESSURE_BASED_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_PRESSURE_BASED_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class FluidPoroMultiPhase;
  class StructPoro;
}  // namespace MAT
namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{

  namespace ELEMENTS
  {
    template <CORE::FE::CellType celltype>
    class SolidPoroPressureBasedEleCalc
    {
     public:
      SolidPoroPressureBasedEleCalc();

      void evaluate_nonlinear_force_stiffness(const CORE::Elements::Element& ele,
          MAT::StructPoro& porostructmat, MAT::FluidPoroMultiPhase& porofluidmat,
          const INPAR::STR::KinemType& kinematictype, const DRT::Discretization& discretization,
          CORE::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix);

      void coupling_poroelast(const CORE::Elements::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& stiffness_matrix);

      void CouplingStress(const CORE::Elements::Element& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      void PoroSetup(MAT::StructPoro& porostructmat, INPUT::LineDefinition* linedef);

     private:
      /// static values for matrix sizes
      static constexpr int num_nodes_ = CORE::FE::num_nodes<celltype>;
      static constexpr int num_dim_ = CORE::FE::dim<celltype>;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      CORE::FE::GaussIntegration gauss_integration_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif