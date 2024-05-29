/*----------------------------------------------------------------------*/
/*! \file

\brief strategy for cloning thermo discretization from scatra discretization


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STI_CLONESTRATEGY_HPP
#define FOUR_C_STI_CLONESTRATEGY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace CORE::Elements
{
  class Element;
}

namespace STI
{
  /*!
  \brief strategy for cloning thermo discretization from scatra discretization

  For a scatra-thermo interaction problem, the thermo discretization is obtained through cloning
  from the scatra discretization. For convenience, we solve the thermo field as a second scatra
  field, hence the thermo discretization is in fact another scatra discretization. Unlike in other
  instances where a scatra discretization is obtained through cloning from another discretization,
  the physical implementation type of the cloned elements is set directly within the strategy and
  not in a subsequent loop after the cloning itself. We can do this because the physical
  implementation type of the cloned elements is unique and known a priori.

  \date 04/15
  */

  class ScatraThermoCloneStrategy
  {
   public:
    //! constructor
    explicit ScatraThermoCloneStrategy() { return; };

    //! destructor
    virtual ~ScatraThermoCloneStrategy() = default;

    //! check material of cloned element
    void check_material_type(const int matid  //! material of cloned element
    );

    //! return map with original names of conditions to be cloned as key values, and final names of
    //! cloned conditions as mapped values
    virtual std::map<std::string, std::string> conditions_to_copy() const;

    //! decide whether element should be cloned or not, and if so, determine type of cloned element
    bool determine_ele_type(
        CORE::Elements::Element* actele,   //! current element on source discretization
        const bool ismyele,                //! ownership flag
        std::vector<std::string>& eletype  //! vector storing types of cloned elements
    );

    //! provide cloned element with element specific data
    void set_element_data(Teuchos::RCP<CORE::Elements::Element>
                              newele,     //! current cloned element on target discretization
        CORE::Elements::Element* oldele,  //! current element on source discretization
        const int matid,                  //! material of cloned element
        const bool isnurbs                //! nurbs flag
    );
  };  // class ScatraThermoCloneStrategy
}  // namespace STI
FOUR_C_NAMESPACE_CLOSE

#endif
