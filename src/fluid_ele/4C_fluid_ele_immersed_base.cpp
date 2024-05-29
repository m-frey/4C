/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for the immersed fluid elements


\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_immersed_base.hpp"

#include "4C_fluid_ele_immersed.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::FluidTypeImmersedBase::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDIMMERSED") return Teuchos::rcp(new DRT::ELEMENTS::FluidImmersed(id, owner));

  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            rauch 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidImmersedBase::FluidImmersedBase(int id, int owner) : Fluid(id, owner) {}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       rauch 03/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidImmersedBase::FluidImmersedBase(const DRT::ELEMENTS::FluidImmersedBase& old)
    : Fluid(old)
{
  return;
}

FOUR_C_NAMESPACE_CLOSE
