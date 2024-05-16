/*----------------------------------------------------------------------*/
/*! \file

\brief Split conditions into map extrators

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_CONDITION_SELECTOR_HPP
#define FOUR_C_DISCRETIZATION_CONDITION_SELECTOR_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class MultiMapExtractor;
}

namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT

namespace CORE::Conditions
{
  /// Select nodes (and their dofs) that are covered by a condition
  /*!
    We have nodal clouds. Each nodal cloud can have an arbitrary number of
    conditions. And in turn one condition (name) can be assigned to any number
    of nodal clouds. Oftentimes, however, we are not so much concerned with
    nodal clouds. We just want to known if a given node is covered by a
    specific condition. This is what this class is for.

    A stack of ConditionSelector objects is used to build a
    MultiConditionSelector which in turn is used to setup a
    CORE::LINALG::MultiMapExtractor object.

    To put it the other way: A CORE::LINALG::MultiMapExtractor splits a full (row)
    map into non-overlapping parts. If each part corresponds to a Condition, a
    MultiConditionSelector can be used to setup the CORE::LINALG::MultiMapExtractor,
    where each Condition is found by a ConditionSelector object.

    \note The Condition objects know the nodes. The maps contain dofs. The
    ConditionSelector translates between nodes and dofs.

    \see CORE::Conditions::MultiConditionSelector
   */
  class ConditionSelector
  {
   public:
    /// Construct a selector on the given Discretization for the Condition
    /// with the given name
    ConditionSelector(const DRT::Discretization& dis, std::string condname);

    /// construct a selector from a given vector of conditions
    ConditionSelector(const DRT::Discretization& dis,  //!< discretization
        const std::vector<Condition*>& conds           //!< given vector of conditions
    );

    /// Destructor
    virtual ~ConditionSelector() = default;

    /// select all matching dofs of the node and put them into conddofset
    virtual bool SelectDofs(DRT::Node* node, std::set<int>& conddofset);

    /// tell if the node gid is known by any condition of the given name
    virtual bool ContainsNode(int ngid);

    /// tell if the dof of a node from this condition is covered as well
    virtual bool ContainsDof(int dof, int pos) { return true; }

   protected:
    /// Discretization we are looking at
    const DRT::Discretization& Discretization() const { return dis_; }

    /// all conditions that come by the given name
    const std::vector<Condition*>& Conditions() const { return conds_; }

   private:
    /// Discretization
    const DRT::Discretization& dis_;

    /// Conditions
    std::vector<Condition*> conds_;
  };


  /// Select some dofs of the conditioned node
  /*!
    In addition to the assignment of nodes we might want to distinguish
    between different dofs. Sometimes a condition applies only to some dofs of
    each node. This is what the ContainsDof() method is for.

    This selector is used most often. It can be applied e.g. to extract the
    velocity dofs from a fluid node (with velocity and pressure dofs)
   */
  class NDimConditionSelector : public ConditionSelector
  {
   public:
    NDimConditionSelector(
        const DRT::Discretization& dis, std::string condname, int startdim, int enddim)
        : ConditionSelector(dis, std::move(condname)), startdim_(startdim), enddim_(enddim)
    {
    }

    /// Contain a dof number only if the dof nodal position is within the
    /// allowed range.
    bool ContainsDof(int dof, int pos) override { return startdim_ <= pos and pos < enddim_; }

   private:
    int startdim_;
    int enddim_;
  };

  /// a collection of ConditionSelector objects used to create a CORE::LINALG::MultiMapExtractor
  /*!
    Oftentimes the dofs of a field need to be split into a set of disjoint
    maps. The CORE::LINALG::MultiMapExtractor class takes care of these splits. However,
    CORE::LINALG::MultiMapExtractor does not do the splitting itself. This is where
    MultiConditionSelector comes in. Here, different ConditionSelector objects
    are collected and afterwards the CORE::LINALG::MultiMapExtractor is setup according to the
    conditions.

    MultiConditionSelector is a helper class that is needed during the setup
    only.

    \note Each node ends up in one condition only, independent of whether all
    its dofs are selected or not.
   */
  class MultiConditionSelector
  {
   public:
    MultiConditionSelector();

    /// add a new ConditionSelector
    /*!
      \note The order of the selector additions determines the slots within
      CORE::LINALG::MultiMapExtractor.
     */
    void AddSelector(Teuchos::RCP<ConditionSelector> s) { selectors_.push_back(s); }

    /// Do the setup
    void SetupExtractor(const DRT::Discretization& dis, const Epetra_Map& fullmap,
        CORE::LINALG::MultiMapExtractor& extractor);

    /// Activate overlapping
    void SetOverlapping(bool overlapping) { overlapping_ = overlapping; }

   private:
    /// evaluate the ConditionSelector objects
    void SetupCondDofSets(const DRT::Discretization& dis);

    /// condition selectors
    std::vector<Teuchos::RCP<ConditionSelector>> selectors_;

    /// sets of selected dof numbers
    std::vector<std::set<int>> conddofset_;

    /// flag defines if maps are overlapping
    bool overlapping_;
  };

}  // namespace CORE::Conditions

FOUR_C_NAMESPACE_CLOSE

#endif
