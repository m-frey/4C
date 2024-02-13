/*--------------------------------------------------------------------------*/
/*! \file

\brief FSI Wrapper for the ALE time integration with internal mesh tying or mesh sliding interface


\level 3

*/
/*--------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_adapter_ale_fsi_msht.hpp"

#include "baci_ale_utils_mapextractor.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFsiMshtWrapper::AleFsiMshtWrapper(Teuchos::RCP<Ale> ale) : AleFsiWrapper(ale)
{
  // create the FSI interface
  fsiinterface_ = Teuchos::rcp(new ALE::UTILS::FsiMapExtractor);
  fsiinterface_->Setup(*Discretization());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::FsiMapExtractor> ADAPTER::AleFsiMshtWrapper::FsiInterface() const
{
  return fsiinterface_;
}

BACI_NAMESPACE_CLOSE
