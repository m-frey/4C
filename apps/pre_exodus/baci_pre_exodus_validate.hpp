/*----------------------------------------------------------------------*/
/*! \file

\brief validate a given .dat-file


\level 1

Validate a given BACI input file (after all preprocessing steps)

*/

#ifndef BACI_PRE_EXODUS_VALIDATE_HPP
#define BACI_PRE_EXODUS_VALIDATE_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <Teuchos_RCP.hpp>

#include <string>

BACI_NAMESPACE_OPEN

namespace EXODUS
{
  // forward declaration
  class Mesh;
  class ElementBlock;

  //! validate a given datfile
  void ValidateInputFile(const Teuchos::RCP<Epetra_Comm> comm, const std::string datfile);

  //! Check Elements for positive Jacobian and otherwise 'rewind' them
  void ValidateMeshElementJacobians(EXODUS::Mesh& mymesh);

  //! Check for positive Jacobian for Element of distype and otherwise 'rewind' them
  void ValidateElementJacobian(
      EXODUS::Mesh& mymesh, const CORE::FE::CellType distype, Teuchos::RCP<EXODUS::ElementBlock>);

  //! Check Elements of distype with full gauss integration rule for positive det at all gps and
  //! return number of negative dets
  int ValidateElementJacobian_fullgp(
      Mesh& mymesh, const CORE::FE::CellType distype, Teuchos::RCP<ElementBlock> eb);

  //! Check one element for positive Jacobi determinant
  bool PositiveEle(const int& eleid, const std::vector<int>& nodes, const EXODUS::Mesh& mymesh,
      const CORE::LINALG::SerialDenseMatrix& deriv);
  int EleSaneSign(
      const std::vector<int>& nodes, const std::map<int, std::vector<double>>& nodecoords);

  //! Rewind one Element
  std::vector<int> RewindEle(std::vector<int> old_nodeids, const CORE::FE::CellType distype);

}  // namespace EXODUS

BACI_NAMESPACE_CLOSE

#endif  // PRE_EXODUS_VALIDATE_H
