/*----------------------------------------------------------------------*/
/*! \file

\brief Projection of state vectors between overlapping meshes

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_XFEM_MESH_PROJECTOR_HPP
#define FOUR_C_XFEM_MESH_PROJECTOR_HPP

#include "4C_config.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::Elements
{
  class Element;
}

namespace CORE::GEO
{
  class SearchTree;
}

namespace XFEM
{
  class MeshProjector
  {
   public:
    //! ctor
    MeshProjector(Teuchos::RCP<const DRT::Discretization> sourcedis,
        Teuchos::RCP<const DRT::Discretization> targetdis, const Teuchos::ParameterList& params,
        Teuchos::RCP<const Epetra_Vector> sourcedisp = Teuchos::null);

    //! set current displacements of source discretization
    void set_source_position_vector(Teuchos::RCP<const Epetra_Vector> sourcedisp = Teuchos::null);

    //! set state vectors - mandatory for interpolation
    void set_source_state_vectors(std::vector<Teuchos::RCP<const Epetra_Vector>> source_statevecs)
    {
      source_statevecs_ = source_statevecs;
    }

    //! main projection routine (pass a map of the target node ids)
    void Project(std::map<int, std::set<int>>&
                     projection_nodeToDof,  //< node-to-dof map of target nodes demanding projection
        std::vector<Teuchos::RCP<Epetra_Vector>>
            target_statevecs,  //< state vectors of target discretization
        Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

    //! projection routine for projection for all nodes of the target discretization
    void project_in_full_target_discretization(
        std::vector<Teuchos::RCP<Epetra_Vector>> target_statevecs,
        Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

    //! write gmsh output for projection details
    void GmshOutput(int step = 0, Teuchos::RCP<const Epetra_Vector> targetdisp = Teuchos::null);

   private:
    /// determine the search radius for the search tree
    template <CORE::FE::CellType distype>
    void find_search_radius();

    //! build a search tree for elements of source discretization
    void setup_search_tree();

    //! for every node search for a covering element from the source discretization
    void find_covering_elements_and_interpolate_values(
        std::vector<CORE::LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<CORE::LINALG::Matrix<8, 1>>& interpolated_vecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values);

    //! compute position of target node w.r.t. source element and interpolate when covered by it
    template <CORE::FE::CellType distype>
    bool check_position_and_project(const CORE::Elements::Element* src_ele,
        const CORE::LINALG::Matrix<3, 1>& node_xyz, CORE::LINALG::Matrix<8, 1>& interpolatedvec);

    //! communicate nodes demanding reconstruction in a Round-Robin pattern
    void communicate_nodes(std::vector<CORE::LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<CORE::LINALG::Matrix<8, 1>>& interpolated_vecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values);

    /// receive a block in the round robin communication pattern
    void receive_block(
        std::vector<char>& rblock, CORE::COMM::Exporter& exporter, MPI_Request& request);

    /// send a block in the round robin communication pattern
    void send_block(
        std::vector<char>& sblock, CORE::COMM::Exporter& exporter, MPI_Request& request);

    /// pack values in the round robin communication pattern
    void pack_values(std::vector<CORE::LINALG::Matrix<3, 1>>& tar_nodepositions,
        std::vector<CORE::LINALG::Matrix<8, 1>>& interpolated_vecs,
        std::vector<int>& projection_targetnodes, std::vector<int>& have_values,
        std::vector<char>& sblock);

    Teuchos::RCP<const DRT::Discretization> sourcedis_;
    Teuchos::RCP<const DRT::Discretization> targetdis_;

    //! search radius factor
    double searchradius_fac_;

    //! 3D seach tree for embedded discretization
    Teuchos::RCP<CORE::GEO::SearchTree> search_tree_;

    //! min. radius needed for the search tree
    double searchradius_;

    //! map of source node to coordinates (including possible displacements)
    std::map<int, CORE::LINALG::Matrix<3, 1>> src_nodepositions_n_;

    //! state vectors from projection source
    std::vector<Teuchos::RCP<const Epetra_Vector>> source_statevecs_;

    //! map between target node id and parent element id
    std::map<int, int> targetnode_to_parent_;
  };
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
