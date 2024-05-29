/*----------------------------------------------------------------------*/
/*! \file
\brief


\level 2

*--------------------------------------------------------------------------*/

#include "4C_binstrategy_meshfree_bin.hpp"

#include "4C_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
DRT::MESHFREE::MeshfreeBin<ELEMENT>::MeshfreeBin(int id, int owner) : ELEMENT(id, owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
DRT::MESHFREE::MeshfreeBin<ELEMENT>::MeshfreeBin(const DRT::MESHFREE::MeshfreeBin<ELEMENT>& old)
    : ELEMENT(old)
{
  return;
}


/*--------------------------------------------------------------------------*
 | Delete a single node from the element               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
void DRT::MESHFREE::MeshfreeBin<ELEMENT>::DeleteNode(int gid)
{
  for (unsigned int i = 0; i < ELEMENT::nodeid_.size(); i++)
  {
    if (ELEMENT::nodeid_[i] == gid)
    {
      ELEMENT::nodeid_.erase(ELEMENT::nodeid_.begin() + i);
      ELEMENT::node_.erase(ELEMENT::node_.begin() + i);
      return;
    }
  }
  FOUR_C_THROW("Connectivity issues: No node with specified gid to delete in element. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Explicit instantiations                                kronbichler 03/15 |
 *--------------------------------------------------------------------------*/
template class DRT::MESHFREE::MeshfreeBin<CORE::Elements::Element>;
template class DRT::MESHFREE::MeshfreeBin<CORE::Elements::FaceElement>;
template class DRT::MESHFREE::MeshfreeBin<MORTAR::Element>;

FOUR_C_NAMESPACE_CLOSE
