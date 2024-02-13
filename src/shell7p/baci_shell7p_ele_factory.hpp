/*! \file

\brief Factory of Shell elements

\level 1
*/

#ifndef BACI_SHELL7P_ELE_FACTORY_HPP
#define BACI_SHELL7P_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"
#include "baci_lib_element.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  // forward declaration
  class Shell7pEleCalcInterface;
  class Shell7p;

  class Shell7pFactory
  {
   public:
    //! copy constructor
    Shell7pFactory() = default;

    static std::unique_ptr<Shell7pEleCalcInterface> ProvideShell7pCalculationInterface(
        const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech);

   private:
    //! define shell calculation instances dependent on element technology
    template <CORE::FE::CellType distype>
    static std::unique_ptr<Shell7pEleCalcInterface> DefineCalculationInterfaceType(
        const std::set<INPAR::STR::EleTech>& eletech);
  };  // class Shell7pFactory
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif
