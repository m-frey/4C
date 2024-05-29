/*----------------------------------------------------------------------*/
/*! \file

\brief a class to manage an enhanced discretization including varying number
of dofs per node on a fluid discretization for xwall

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DISCRET_XWALL_HPP
#define FOUR_C_LIB_DISCRET_XWALL_HPP

#include "4C_config.hpp"

#include "4C_lib_discret_faces.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class DiscretizationXWall : public DiscretizationFaces
  {
   public:
    /*!
    \brief Standard Constructor

    \param name (in): name of this discretization
    \param comm (in): An epetra comm object associated with this discretization
    */
    DiscretizationXWall(const std::string name, Teuchos::RCP<Epetra_Comm> comm);



    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))

    Additional input nodal dof set: If the node contains more than one set of dofs, which can be
    evaluated, the number of the set needs to be given. Currently only the case for XFEM.

    \param dof (in)         : vector of dof gids (to be filled)
    \param nds (in)         : number of dofset
    \param nodaldofset (in) : number of nodal dofset
    \param node (in)        : the node
    \param element (in)     : the element (optionally)
    */
    void Dof(std::vector<int>& dof, const Node* node, unsigned nds, unsigned nodaldofset,
        const CORE::Elements::Element* element = nullptr) const override
    {
      if (nds > 1) FOUR_C_THROW("xwall discretization can only handle one dofset at the moment");

      FOUR_C_ASSERT(nds < dofsets_.size(), "undefined dof set");
      FOUR_C_ASSERT(havedof_, "no dofs assigned");

      std::vector<int> totaldof;
      dofsets_[nds]->Dof(totaldof, node, nodaldofset);

      if (element == nullptr && element->Shape() == CORE::FE::CellType::hex8)
        FOUR_C_THROW("element required for location vector of hex8 element");

      int size;
      if (element != nullptr)
        size = std::min((int)totaldof.size(), element->NumDofPerNode(*node));
      else
        size = (int)totaldof.size();
      // only take the first dofs that have a meaning for all elements at this node
      for (int i = 0; i < size; i++) dof.push_back(totaldof.at(i));


      return;
    }
  };  // class DiscretizationXFEM
}  // namespace DRT

/// << operator
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationXWall& dis);


FOUR_C_NAMESPACE_CLOSE

#endif
