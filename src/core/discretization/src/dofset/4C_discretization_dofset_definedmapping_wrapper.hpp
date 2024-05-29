/*----------------------------------------------------------------------*/
/*! \file

\brief A dofset that does not rely on same GID/LID numbers but uses
       a defined node mapping instead (not implemented for element DOFs).

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_DOFSET_DEFINEDMAPPING_WRAPPER_HPP
#define FOUR_C_DISCRETIZATION_DOFSET_DEFINEDMAPPING_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset_base.hpp"
#include "4C_lib_node.hpp"

#include <Epetra_IntVector.h>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::Dofsets
{
  /*!
  \brief A dofset that does not rely on same GID/LID numbers but uses
         a defined node mapping instead

  \warning not implemented for element DOFs

  This dofset was originally written to handle heterogeneous reactions :
  \ref SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy .
  This means reactions between volume-bound scalars and surface-bound scalars in the
  context of scalar transport simulations.
  It is also used in solid-scatra interaction \ref SSI::SSI_Base::SSI_Base .


  <h3>assign_degrees_of_freedom</h3>

  The important work is done in \ref assign_degrees_of_freedom . If you want to use
  or extend this class, please be aware of the documentation of this method. Also
  have in mind the needs of those models and methods in 4C this class was originally
  intended to support (see above).


  <h3>Usage</h3>

  After construction the newly created dofset has to be added to the target discretization.
  The target discretization then has to be filled. During \ref fill_complete() the method
  \ref assign_degrees_of_freedom is called.

  One example for the proper use of this class:

  \code

    Teuchos::RCP<CORE::Dofsets::DofSet> newsourcedofset_for_target_dis =
        Teuchos::rcp(new
  CORE::Dofsets::DofSetMappedProxy(sourcedis->GetDofSetProxy(),sourcedis,"CouplingSourceToTarget",setofcouplingids));

    target_dis->AddDofSet(newsourcedofset_for_target_dis);

    target_dis->fill_complete();

  \endcode

  See also the example '.dat'-file 'immersed_cell_biochemo_mechano_pureProtrusion_h8' in
  the 4C 'Input' folder.


  \date 08/2016
  \author rauch
  \author vuong
  */
  class DofSetDefinedMappingWrapper : public DofSetBase
  {
   public:
    /*! \brief Standard Constructor

    \param dofset       (in) : dofset proxy of source discretization \ref sourcedis
    \param sourcedis    (in) : source discretization for the new dofset
    \param couplingcond (in) : condition to be coupled
    \param condids      (in) : set of condition ids to merge into one dofset    */
    DofSetDefinedMappingWrapper(Teuchos::RCP<DofSetInterface> dofset,
        Teuchos::RCP<const DRT::Discretization> sourcedis, const std::string& couplingcond,
        const std::set<int> condids);

    //! destructor
    ~DofSetDefinedMappingWrapper() override;


    /*! \brief Assign dof numbers using all elements and nodes of the discretization.

           This method is called during \ref fill_complete() on the target
           discretization.


          <h3>Usage</h3>

          An important feature of this dofset proxy is, that portions of source and
          target discretizations, whic hare supposed to be coupled by this dofset do
          not necessarily need to have the same number of nodes, as long as the mapping
          is still unique. We give an example to show what this means:

          Imagine, we have given the following conditions in our .dat-file:

          -----DESIGN SSI COUPLING SOLIDTOSCATRA VOL CONDITIONS
          DVOL   2
          // scatra volume matching to struct volume
          E 2 - 1
          // struct volume matching to scatra volume
          E 1 - 1
          ----DESIGN SSI COUPLING SOLIDTOSCATRA SURF CONDITIONS
          DSURF   2
         // separate surface discretization matching struct volume boundary
         E 1 - 2
         // struct volume boundary
         E 2 - 2

         We have two matching volume discretizations for structure and scatra, and
         we have a scatra surface discretization matching the boundary of the struct
         volume. Thus, we have more scatra nodes, than structure nodes, because we
         have NumNodeScatraVolume + NumNodeScatraSurface scatra nodes, but only
         NumNodeStructVolume (= NumNodeScatraVolume) structure nodes.
         However, since we are matching, each scatra node can be mapped uniquiely to
         a structure node. Note that we could not map each structure node uniquely
         to a scatra node because at the boundary two scatra nodes reside at the same
         location (volume + surface)!
         So in the first case, if we want to add the structure dofset to the scatradis,
         the structure is the source and the scatra is the target.
         We can now create the dofset proxy from the struct dis via

         \code
          Teuchos::RCP<CORE::Dofsets::DofSet> newdofset =
              Teuchos::rcp(new
       CORE::Dofsets::DofSetMappedProxy(structdis->GetDofSetProxy(), structdis,
                                                      "SSICouplingSolidToScatra",
                                                      couplingids));
         \endcode

         The std::set 'couplingids' contains the ids from both, VOL and SURF conditions
         given above. This way, one unique struct dofset for all scatra nodes is created.    */
    int assign_degrees_of_freedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start) override;

    /// reset all internal variables
    void Reset() override;

    //! @name Proxy management
    /// Proxies need to know about changes to the DofSet.

    /// our original DofSet dies
    void Disconnect(DofSetInterface* dofset) override;

    //@}

    /// Returns true if filled
    bool Filled() const override { return filled_ and sourcedofset_->Filled(); };


    //! @name Access methods

    /// Get degree of freedom row map
    const Epetra_Map* dof_row_map() const override { return sourcedofset_->dof_row_map(); };

    /// Get degree of freedom column map
    const Epetra_Map* DofColMap() const override { return sourcedofset_->DofColMap(); };

    /// Get number of dofs for given node
    int NumDof(const DRT::Node* node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        return sourcedofset_->NumDof(node);
      else
        return sourcedofset_->NumDof(sourcenode);
    }

    /// Get number of dofs for given element
    int NumDof(const CORE::Elements::Element*
            element  ///< element, for which you want to know the number of dofs
    ) const override
    {
      // element dofs not yet supported
      return 0;
    }

    /// get number of nodal dofs
    int NumDofPerNode(
        const DRT::Node& node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node.LID());
      if (sourcenode == nullptr)
        return sourcedofset_->NumDofPerNode(node);
      else
        return sourcedofset_->NumDofPerNode(*sourcenode);
    }

    /// Get the gid of all dofs of a node
    std::vector<int> Dof(
        const DRT::Node* node  ///< node, for which you want to know the number of dofs
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        return sourcedofset_->Dof(node);
      else
        return sourcedofset_->Dof(sourcenode);
    }

    /// Get the gid of all dofs of a node
    void Dof(std::vector<int>& dof,  ///< vector of dof gids (to be filled)
        const DRT::Node* node,       ///< node, for which you want the dof positions
        unsigned nodaldofset  ///< number of nodal dof set of the node (currently !=0 only for XFEM)
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        sourcedofset_->Dof(dof, node, nodaldofset);
      else
        sourcedofset_->Dof(dof, sourcenode, nodaldofset);
    };

    /// Get the gid of all dofs of a element
    std::vector<int> Dof(const CORE::Elements::Element* element) const override
    {
      return std::vector<int>();
    }

    /// Get the gid of a dof for given node
    int Dof(const DRT::Node* node, int dof) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        return sourcedofset_->Dof(node, dof);
      else
        return sourcedofset_->Dof(sourcenode, dof);
    }

    /// Get the gid of a dof for given element
    int Dof(const CORE::Elements::Element* element, int dof) const override
    {
      // element dofs not yet supported
      return 0;
    }

    /// Get the gid of all dofs of a node
    void Dof(const DRT::Node* node, std::vector<int>& lm) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        return sourcedofset_->Dof(node, lm);
      else
        return sourcedofset_->Dof(sourcenode, lm);
    }

    /// Get the gid of all dofs of a node
    void Dof(const DRT::Node* node,  ///< node, for which you want the dof positions
        const unsigned startindex,   ///< first index of vector at which will be written to end
        std::vector<int>& lm         ///< already allocated vector to be filled with dof positions
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        sourcedofset_->Dof(node, startindex, lm);
      else
        sourcedofset_->Dof(sourcenode, startindex, lm);
    }

    /// Get the gid of all dofs of a element
    void Dof(const CORE::Elements::Element* element, std::vector<int>& lm) const override
    {
      // element dofs not yet supported
      return;
    }

    /// Get the GIDs of the first DOFs of a node of which the associated element is interested in
    void Dof(const CORE::Elements::Element*
                 element,       ///< element which provides its expected number of DOFs per node
        const DRT::Node* node,  ///< node, for which you want the DOF positions
        std::vector<int>& lm    ///< already allocated vector to be filled with DOF positions
    ) const override
    {
      const DRT::Node* sourcenode = get_source_node(node->LID());
      if (sourcenode == nullptr)
        sourcedofset_->Dof(element, node, lm);
      else
        sourcedofset_->Dof(element, sourcenode, lm);
    }

    /// Get maximum GID of degree of freedom row map
    int MaxAllGID() const override { return sourcedofset_->MaxAllGID(); };

    /// Get minimum GID of degree of freedom row map
    int MinAllGID() const override { return sourcedofset_->MinAllGID(); };

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int MaxGIDinList(const Epetra_Comm& comm) const override
    {
      return sourcedofset_->MaxGIDinList(comm);
    };

    /// are the dof maps already initialized?
    bool Initialized() const override { return sourcedofset_->Initialized(); };

    /// Get Number of Global Elements of degree of freedom row map
    int NumGlobalElements() const override { return sourcedofset_->NumGlobalElements(); };

   private:
    //! get corresponding source node from source discretization
    const DRT::Node* get_source_node(int targetLid) const;

    // dofset
    Teuchos::RCP<DofSetInterface> sourcedofset_;

    //! map containing the mapping of the target node GID to the corresponding source node GID
    //! (value)
    Teuchos::RCP<Epetra_IntVector> targetlidtosourcegidmapping_;

    //! source discretization
    Teuchos::RCP<const DRT::Discretization> sourcedis_;

    //! condition string defining the coupling
    const std::string couplingcond_;

    //! ID of condition the dofset is build from
    const std::set<int> condids_;

    /// filled flag
    bool filled_;
  };  // DofSetDefinedMappingWrapper
}  // namespace CORE::Dofsets



FOUR_C_NAMESPACE_CLOSE

#endif
