/*----------------------------------------------------------------------*/
/*! \file
\brief Center coordinates of an element
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENT_CENTER_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENT_CENTER_HPP

#include "4C_config.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::Elements
{
  class Element;
}

namespace CORE::FE
{
  // return center coordinates of @p ele. Note that the average of all nodal coordinates is
  // computed.
  std::vector<double> element_center_refe_coords(const CORE::Elements::Element& ele);
}  // namespace CORE::FE

#endif

FOUR_C_NAMESPACE_CLOSE
