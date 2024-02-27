/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-thermo interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_STI_HPP
#define BACI_INPAR_STI_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR
{
  namespace STI
  {
    //! type of coupling between scatra and thermo fields
    enum class CouplingType
    {
      undefined,
      monolithic,
      oneway_scatratothermo,
      oneway_thermotoscatra,
      twoway_scatratothermo,
      twoway_scatratothermo_aitken,
      twoway_scatratothermo_aitken_dofsplit,
      twoway_thermotoscatra,
      twoway_thermotoscatra_aitken
    };

    //! type of scalar transport time integration
    enum class ScaTraTimIntType
    {
      standard,
      elch
    };

    //! set valid parameters for scatra-thermo interaction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    //! set valid conditions for scatra-thermo interaction
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);
  }  // namespace STI
}  // namespace INPAR
BACI_NAMESPACE_CLOSE

#endif  // INPAR_STI_H
