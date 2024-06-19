/*----------------------------------------------------------------------------*/
/** \file

\brief MultiFieldMapExtractor class to handle different discretizations
       with joint interfaces


\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_MULTI_FIELD_MAPEXTRACTOR_HPP
#define FOUR_C_XFEM_MULTI_FIELD_MAPEXTRACTOR_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_xfem_enum_lists.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>

// forward declarations
class Epetra_Map;

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Nodes
{
  class Node;
}

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MultiMapExtractor;
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace Core::LinAlg

namespace XFEM
{
  namespace XFieldField
  {
    class CouplingDofSet;
    class Coupling;
  }  // namespace XFieldField

  /** \class MultiFieldMapExtractor
   *  \brief MultiFieldMapExtractor class to handle different discretizations
   *  with joint interfaces
   *
   *  This class is supposed to be used for problems including one or more
   *  standard field discretizations and one or more xFEM field
   *  discretizations. These discretizations should be generated by splitting
   *  one large discretization, such that the joint interfaces of two (or more)
   *  discretizations are node-matching. The XFEM discretizations are allowed
   *  to be enriched at a number of interface nodes. This is considered during
   *  the creation process of the common unique DoF row map. Furthermore, the
   *  discretizations are allowed to be distributed differently over the
   *  processors (also at the joint interfaces).
   *
   *  The objective of this class is, to create a full map which can be used to
   *  create the necessary state vectors and system matrices and to simplify
   *  the communication between the different discretizations by providing
   *  appropriate extract and insert methods.
   *
   *  Simple 2-dimensional example:
   *
   *  We start with the following discretization. (+) nodes indicate standard
   *  nodes, (o) indicate enriched (or possible enriched) nodes.
   *
   *      +---+---+    nodal GID's from left to right:  0  1  2
   *      |   |   |
   *      o---+---+                                     3  4  5
   *      |   |   |
   *      o---+---+                                     6  7  8
   *      |   |   |
   *      o---+---+                                     9 10 11
   *      |   |   |
   *      o---+---+                                    12 13 14
   *      |   |   |
   *      +---+---+                                    15 16 17
   *
   *  Before you can call this class, you have to split the discretization into
   *  a standard and a xFEM discretization. This can be achieved in different
   *  ways. One possibility is to use the XFEM::UTILS::XFEMDiscretizationBuilder.
   *  Both discretizations are supposed to share still the same nodal global
   *  ID's at their joint interfaces.
   *
   *      +---+---+   0  1  2
   *      |   |   |
   *      +---+---+   3  4  5                   o---+   3  4
   *          |   |                             |   |
   *          +---+      7  8                   o---+   6  7
   *          |   |                 AND         |   |
   *          +---+     10 11                   o---+   9 10
   *          |   |                             |   |
   *      +---+---+  12 13 14                   o---+  12 13
   *      |   |   |
   *      +---+---+  15 16 17
   *
   *      STANDARD                              XFEM
   *
   *  Now this class starts to find all interface nodes between the different
   *  discretizations and creates an auxiliary interface discretization
   *  containing only these shared nodes.
   *
   *               o   +   3  4
   *
   *                   +      7
   *
   *                   +     10
   *
   *               o   +  12 13
   *
   *  The helping discretization is allowed to be distributed differently,
   *  although it is not recommended because of the resulting communication
   *  overhead.
   *
   *  \author hiermeier
   *  \date 09/16 */
  class MultiFieldMapExtractor
  {
    typedef std::vector<Teuchos::RCP<const Core::FE::Discretization>> XDisVec;

    // number of map extractor types
    static constexpr unsigned NUM_MAP_TYPES = 2;

   public:
    /// constructor
    MultiFieldMapExtractor();

    /// destructor
    virtual ~MultiFieldMapExtractor() = default;

    /** \brief Initialize the MultiMapExtractor using a set of filled discretizations
     *
     *  Everything is initialized what is independent of any possible up-coming
     *  redistribution.
     *
     *  \param dis_vec (in): Vector containing different in some way connected
     *                       standard and/or xFEM discretizations.
     *  \param max_num_reserved_dofs_per_node (in): This is the maximal possible number
     *                       of DoF's per enriched node (necessary for the fixed size
     *                       dofset).
     *
     *  \author hiermeier
     *  \date 09/16 */
    void init(const std::vector<Teuchos::RCP<const Core::FE::Discretization>>& dis_vec,
        int max_num_reserved_dofs_per_node);

    /** \brief Setup member variables
     *
     *  Everything is initialized what will change if one of the
     *  discretizations in the discretization vector (see init()) is
     *  redistributed.
     *
     *  \author hiermeier
     *  \date 09/16 */
    virtual void setup();

    /// @name Accessors to the auxiliary interface discretization
    /// @{

    /** \brief access the interface node with given global id
     *
     *  \param gid (in): global id of the interface node
     *
     *  \author hiermeier
     *  \date 10/16 */
    Core::Nodes::Node* gINode(const int& gid) const;

    /// access the interface node row map
    const Epetra_Map* INodeRowMap() const;

    /** \brief get the number of DoF's of the given interface node
     *
     *  Returns the maximum number of DoF's of the jointing discretizations
     *  at the interface node
     *
     *  \param inode (in): pointer to the interface node
     *
     *  \author hiermeier
     *  \date 10/16 */
    int INumDof(const Core::Nodes::Node* inode) const;

    /// get the number of standard DoF's of the discretization
    int INumStandardDof() const;


    int IDof(const Core::Nodes::Node* inode, int dof) const;

    void IDof(const Core::Nodes::Node* inode, std::vector<int>& dofs) const;

    void IDof(std::vector<int>& dof, Core::Nodes::Node* inode, unsigned nodaldofset_id,
        const Core::Elements::Element* element) const;

    /// @}

    const Core::LinAlg::MultiMapExtractor& SlDofMapExtractor(enum FieldName field) const
    {
      return sl_map_extractor(slave_id(field), map_dofs);
    }

    Teuchos::RCP<const Epetra_Map> NodeRowMap(
        enum FieldName field, enum MultiField::BlockType block) const;

    /** \brief return TRUE if the given global node id corresponds to an
     *  interface node
     *
     *  \author hiermeier
     *  \date 09/16 */
    bool IsInterfaceNode(const int& ngid) const;

    /// Access the full maps
    const Teuchos::RCP<const Epetra_Map>& FullMap(enum MapType map_type = map_dofs) const;

    /// @name Extract vector routines
    /// @{
    Teuchos::RCP<Epetra_Vector> ExtractVector(
        const Epetra_Vector& full, enum FieldName field, enum MapType map_type = map_dofs) const;

    Teuchos::RCP<Epetra_MultiVector> ExtractVector(const Epetra_MultiVector& full,
        enum FieldName field, enum MapType map_type = map_dofs) const;

    inline void ExtractVector(Teuchos::RCP<const Epetra_Vector> full, enum FieldName field,
        Teuchos::RCP<Epetra_Vector> partial, enum MapType map_type = map_dofs) const
    {
      ExtractVector(*full, field, *partial, map_type);
    }

    inline void ExtractVector(const Epetra_MultiVector& full, enum FieldName field,
        Epetra_MultiVector& partial, enum MapType map_type = map_dofs) const
    {
      ExtractVector(full, slave_id(field), partial, map_type);
    }

    void ExtractVector(const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial,
        enum MapType map_type = map_dofs) const;

    inline void extract_element_vector(
        const Epetra_MultiVector& full, enum FieldName field, Epetra_MultiVector& partial) const
    {
      extract_element_vector(full, slave_id(field), partial);
    }

    void extract_element_vector(
        const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const;
    /// @}

    /// @name Routines to insert a partial vector into a full vector
    /// @{

    /** \brief Put a partial vector into a full Epetra_Vector
     *
     *  \param partial (in): vector to copy into full vector (Epetra_Vector)
     *  \param field   (in): field name enumerator of the partial vector
     *
     *  \author hiermeier \date 10/16 */
    Teuchos::RCP<Epetra_Vector> InsertVector(
        const Epetra_Vector& partial, enum FieldName field, enum MapType map_type = map_dofs) const;

    /** \brief Put a partial vector into a full vector (Epetra_MultiVector)
     *
     *  \param partial (in): vector to copy into full vector
     *  \param field   (in): field name enumerator of the partial vector
     *
     *  \author hiermeier \date 10/16 */
    Teuchos::RCP<Epetra_MultiVector> InsertVector(const Epetra_MultiVector& partial,
        enum FieldName field, enum MapType map_type = map_dofs) const;

    /** \brief Put a partial vector into a full vector (Epetra_Vector)
     *
     *  \param partial (in): vector to copy into full vector
     *  \param field   (in): field name enumerator of the partial vector
     *  \param full   (out): vector to copy into
     *
     *  \author hiermeier \date 10/16 */
    void InsertVector(Teuchos::RCP<const Epetra_Vector> partial, enum FieldName field,
        Teuchos::RCP<Epetra_Vector> full, enum MapType map_type = map_dofs) const
    {
      InsertVector(*partial, field, *full, map_type);
    }

    /** \brief Put a partial vector into a full vector (Epetra_MultiVector)
     *
     *  \param partial (in): vector to copy into full vector
     *  \param field   (in): field name enumerator of the partial vector
     *  \param full   (out): vector to copy into
     *
     *  \author hiermeier \date 10/16 */
    void InsertVector(const Epetra_MultiVector& partial, enum FieldName field,
        Epetra_MultiVector& full, enum MapType map_type = map_dofs) const
    {
      return InsertVector(partial, slave_id(field), full, map_type);
    }

    /** \brief Put a partial vector into a full vector (Epetra_MultiVector) [derived]
     *
     *  \author hiermeier \date 10/16  */
    void InsertVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full,
        enum MapType map_type = map_dofs) const;

    inline void InsertElementVector(
        const Epetra_MultiVector& partial, enum FieldName field, Epetra_MultiVector& full) const
    {
      InsertElementVector(partial, slave_id(field), full);
    }

    void InsertElementVector(
        const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const;
    /// @}

    /// @name Routines to add a partial vector to the full vector
    /// @{
    /** \brief Add a partial vector to a full vector (Epetra_Vector)
     *
     *  \param partial (in): vector which is added to full vector
     *  \param field   (in): field name enumerator of the partial vector
     *  \param full   (out): sum into this full vector
     *  \param scale   (in): scaling factor for partial vector
     *
     *  \author hiermeier \date 10/16 */
    inline void AddVector(Teuchos::RCP<const Epetra_Vector> partial, enum FieldName field,
        Teuchos::RCP<Epetra_Vector> full, double scale, enum MapType map_type = map_dofs) const
    {
      AddVector(*partial, slave_id(field), *full, scale, map_type);
    }

    /** \brief Add a partial vector to a full vector (Epetra_MultiVector)
     *
     *  \param partial (in): vector which is added to the full vector
     *  \param field   (in): field name enumerator of the partial vector
     *  \param full   (out): sum into this full vector
     *  \param scale   (in): scaling factor for partial vector
     *
     *  \author hiermeier \date 10/16 */
    inline void AddVector(const Epetra_MultiVector& partial, enum FieldName field,
        Epetra_MultiVector& full, double scale, enum MapType map_type = map_dofs) const
    {
      return AddVector(partial, slave_id(field), full, scale, map_type);
    }

    /** \brief Add a partial vector to a full vector (Epetra_MultiVector) [derived]
     *
     *  \author hiermeier \date 10/16 */
    void AddVector(const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full,
        double scale, enum MapType map_type = map_dofs) const;

    inline void AddElementVector(const Epetra_MultiVector& partial, enum FieldName field,
        Epetra_MultiVector& full, double scale) const
    {
      AddElementVector(partial, slave_id(field), full, scale);
    }

    void AddElementVector(
        const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale) const;
    /// @}

    /// @name Add a partial system-matrix to the full matrix
    /// @{
    inline void AddMatrix(const Core::LinAlg::SparseOperator& partial_mat, enum FieldName field,
        Core::LinAlg::SparseOperator& full_mat, double scale)
    {
      AddMatrix(partial_mat, slave_id(field), full_mat, scale);
    }

    void AddMatrix(const Core::LinAlg::SparseOperator& partial_mat, int block,
        Core::LinAlg::SparseOperator& full_mat, double scale);

    void AddMatrix(const Core::LinAlg::BlockSparseMatrixBase& partial_mat, int block,
        Core::LinAlg::SparseMatrix& full_mat, double scale);

    /// @}

    /** \brief return TRUE if the discretization dis_id is a XFEM discretization
     *
     *  \param dis_id (in): entry of the slave discretization vector
     *
     *  \author hiermeier
     *  \date 09/16 */
    bool is_x_fem_dis(enum FieldName field) const { return is_x_fem_dis(slave_id(field)); }

   protected:
    /// check if init() has been called yet
    inline void check_init() const
    {
      if (not isinit_) FOUR_C_THROW("Call init() first!");
    }

    /// check if init() and setup() have been called yet
    inline void check_init_setup() const
    {
      if ((not isinit_) or (not issetup_)) FOUR_C_THROW("Call init() and/or setup() first!");
    }

   private:
    /** \brief return TRUE if the discretization dis_id is a XFEM discretization
     *
     *  \param dis_id (in): entry of the slave discretization vector
     *
     *  \author hiermeier
     *  \date 09/16 */
    bool is_x_fem_dis(int dis_id) const;

    /** \brief  Access the master interface node row map of the interface
     *  between the master interface discretization and the slave discretization
     *  with the ID dis_id
     *
     *  \param dis_id (in): entry of the slave discretization vector
     *
     *  \author hiermeier
     *  \date 09/16 */
    inline const Epetra_Map& master_interface_node_row_map(enum FieldName field) const
    {
      return master_interface_node_row_map(slave_id(field));
    }
    const Epetra_Map& master_interface_node_row_map(unsigned dis_id) const
    {
      check_init();

      if (dis_id >= master_interface_node_maps_.size())
        FOUR_C_THROW(
            "The index %d exceeds the master interface node row map size! "
            "(size = %d)",
            dis_id, master_interface_node_maps_.size());

      if (master_interface_node_maps_[dis_id].is_null())
        FOUR_C_THROW(
            "The master interface node row map %d was not initialized "
            "correctly.",
            dis_id);

      return *(master_interface_node_maps_[dis_id]);
    }

    /** \brief Access the master map extractor
     *
     *  \author hiermeier
     *  \date 09/16 */
    const Core::LinAlg::MultiMapExtractor& ma_map_extractor(enum MapType map_type) const
    {
      if (master_map_extractor_.at(map_type).is_null())
        FOUR_C_THROW("The master dof/node map extractor was not initialized!");

      return *(master_map_extractor_[map_type]);
    }

    /** \brief Access the slave sided node row maps
     *
     *  \param dis_id (in): block id of the desired discretization
     *  \param btype  (in): choose between interface and non-interface nodes
     *
     *  \author hiermeier
     *  \date 10/16 */
    inline const Epetra_Map& slave_node_row_map(
        enum XFEM::FieldName field, enum MultiField::BlockType btype) const
    {
      return slave_node_row_map(slave_id(field), btype);
    }
    const Epetra_Map& slave_node_row_map(unsigned dis_id, enum MultiField::BlockType btype) const;

    const Core::LinAlg::MultiMapExtractor& sl_map_extractor(
        unsigned dis_id, enum MapType map_type) const
    {
      check_init();

      if (dis_id >= slave_map_extractors_.size())
        FOUR_C_THROW(
            "The index %d exceeds the slave map extractor size! "
            "(size = %d)",
            dis_id, slave_map_extractors_.size());

      if (slave_map_extractors_[dis_id].at(map_type).is_null())
        FOUR_C_THROW(
            "The slave dof/node map extractor %d was not initialized "
            "correctly.",
            dis_id);

      return *(slave_map_extractors_[dis_id][map_type]);
    }

    /** \brief Access the interface matrix row transformer for the given field
     *
     *  \author hiermeier \date 10/16 */
    Core::LinAlg::MatrixRowTransform& i_mat_row_transform(enum FieldName field)
    {
      return i_mat_row_transform(slave_id(field));
    }
    Core::LinAlg::MatrixRowTransform& i_mat_row_transform(unsigned dis_id)
    {
      check_init();

      if (dis_id >= interface_matrix_row_transformers_.size())
        FOUR_C_THROW(
            "The index %d exceeds the matrix row transformer size! "
            "(size = %d)",
            dis_id, interface_matrix_row_transformers_.size());

      if (interface_matrix_row_transformers_[dis_id].is_null())
        FOUR_C_THROW(
            "The interface matrix row transformer %d was not initialized "
            "correctly.",
            dis_id);

      return *(interface_matrix_row_transformers_[dis_id]);
    }

    /** \brief Access the interface matrix column transformer for the given field
     *
     *  \author hiermeier \date 10/16 */
    Core::LinAlg::MatrixColTransform& i_mat_col_transform(enum FieldName field)
    {
      return i_mat_col_transform(slave_id(field));
    }
    Core::LinAlg::MatrixColTransform& i_mat_col_transform(unsigned dis_id)
    {
      check_init();

      if (dis_id >= interface_matrix_col_transformers_.size())
        FOUR_C_THROW(
            "The index %d exceeds the matrix column transformer size! "
            "(size = %d)",
            dis_id, interface_matrix_col_transformers_.size());

      if (interface_matrix_col_transformers_[dis_id].is_null())
        FOUR_C_THROW(
            "The interface matrix column transformer %d was not initialized "
            "correctly.",
            dis_id);

      return *(interface_matrix_col_transformers_[dis_id]);
    }

    /** \brief Access the interface matrix row and column transformer for the given field
     *
     *  \author hiermeier \date 10/16 */
    Core::LinAlg::MatrixRowColTransform& i_mat_row_col_transform(enum FieldName field)
    {
      return i_mat_row_col_transform(slave_id(field));
    }
    Core::LinAlg::MatrixRowColTransform& i_mat_row_col_transform(unsigned dis_id)
    {
      check_init();

      if (dis_id >= interface_matrix_row_col_transformers_.size())
        FOUR_C_THROW(
            "The index %d exceeds the matrix row col transformer size! "
            "(size = %d)",
            dis_id, interface_matrix_row_col_transformers_.size());

      if (interface_matrix_row_col_transformers_[dis_id].is_null())
        FOUR_C_THROW(
            "The interface matrix row col transformer %d was not initialized "
            "correctly.",
            dis_id);

      return *(interface_matrix_row_col_transformers_[dis_id]);
    }

    /** \brief Access the interface discretization
     *
     *  \author hiermeier \date 10/16 */
    inline const Core::FE::Discretization& i_discret() const
    {
      check_init();
      return *idiscret_;
    }

    inline const Core::FE::Discretization& sl_discret(enum FieldName field) const
    {
      return sl_discret(slave_id(field));
    }
    const Core::FE::Discretization& sl_discret(unsigned dis_id) const
    {
      check_init();

      if (dis_id >= num_sl_dis())
        FOUR_C_THROW(
            "The index %d exceeds the slave discretization vector size! "
            "(size = %d)",
            dis_id, sl_dis_vec().size());

      if (sl_dis_vec()[dis_id].is_null())
        FOUR_C_THROW(
            "The slave discretization %d was not initialized "
            "correctly.",
            dis_id);

      return *(sl_dis_vec()[dis_id]);
    }

    const XFEM::XFieldField::Coupling& i_coupling(unsigned dis_id) const
    {
      check_init();
      if (dis_id >= interface_couplings_.size())
        FOUR_C_THROW(
            "The index %d exceeds the interface coupling size! "
            "(size = %d)",
            dis_id, interface_couplings_.size());
      if (interface_couplings_[dis_id].is_null())
        FOUR_C_THROW(
            "The interface coupling %d was not initialized "
            "correctly.",
            dis_id);

      return *(interface_couplings_[dis_id]);
    }

    inline const Epetra_Comm& comm() const
    {
      if (comm_.is_null()) FOUR_C_THROW("The Epetra_Comm object has not been initialized!");

      return *comm_;
    }

    inline unsigned num_sl_dis() const { return sl_dis_vec().size(); }

    const std::set<int>& g_interface_node_gid_set() const { return g_interface_node_gid_set_; }

    /** \brief reset class variables at the beginning of each init() and setup() call
     *
     *  \param num_dis_vec (in): number of wrapped discretizations
     *  \param full        (in): TRUE initiates a reset of all class variables */
    void reset(unsigned num_dis_vec) { reset(num_dis_vec, true); }
    void reset(unsigned num_dis_vec, bool full);

    /// get the row node/DoF maps of the wrapped discretizations
    void get_dof_and_node_maps();

    void build_global_interface_node_gid_set();

    void build_master_interface_node_maps(
        const std::vector<std::vector<int>>& my_master_interface_node_gids);

    void build_slave_discret_id_map();

    int slave_id(enum FieldName field) const;

    const std::vector<Teuchos::RCP<const Core::FE::Discretization>>& sl_dis_vec() const
    {
      return slave_discret_vec_;
    }

    void build_slave_dof_map_extractors();

    void build_slave_node_map_extractors();

    void build_master_node_map_extractor();

    void build_master_dof_map_extractor();

    void build_element_map_extractor();

    /** \brief Build the interface coupling DoF set and complete the interface
     *  discretization
     *
     *  \author hiermeier
     *  \date 09/16 */
    void build_interface_coupling_dof_set();

    void build_interface_coupling_adapters();

    void build_interface_matrix_transformers();

   private:
    /// boolean which indicates, that the init() routine has been called
    bool isinit_;

    /// boolean which indicates, that the setup() routine has been called
    bool issetup_;

    int max_num_reserved_dofs_per_node_;

    /// Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm_;

    /// vector containing pointers to all the input discretizations
    std::vector<Teuchos::RCP<const Core::FE::Discretization>> slave_discret_vec_;

    /// mapping between the FieldName enumerator and the slave vector entry number
    std::map<enum FieldName, int> slave_discret_id_map_;

    /** \brief global interface node GID set
     *
     * (containing the same information on all proc's) */
    std::set<int> g_interface_node_gid_set_;

    std::vector<Teuchos::RCP<const Epetra_Map>> master_interface_node_maps_;

    std::vector<std::vector<Teuchos::RCP<Core::LinAlg::MultiMapExtractor>>> slave_map_extractors_;
    std::vector<Teuchos::RCP<Core::LinAlg::MultiMapExtractor>> master_map_extractor_;

    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> element_map_extractor_;

    std::vector<Teuchos::RCP<XFEM::XFieldField::Coupling>> interface_couplings_;

    std::vector<Teuchos::RCP<Core::LinAlg::MatrixRowTransform>> interface_matrix_row_transformers_;
    std::vector<Teuchos::RCP<Core::LinAlg::MatrixColTransform>> interface_matrix_col_transformers_;
    std::vector<Teuchos::RCP<Core::LinAlg::MatrixRowColTransform>>
        interface_matrix_row_col_transformers_;

    std::set<int> xfem_dis_ids_;

    /// interface discretization
    Teuchos::RCP<Core::FE::Discretization> idiscret_;

    /// interface coupling DoF-set
    Teuchos::RCP<XFEM::XFieldField::CouplingDofSet> icoupl_dofset_;

  };  // class MapExtractor
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
