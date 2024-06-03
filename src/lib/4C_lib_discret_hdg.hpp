/*---------------------------------------------------------------------*/
/*! \file

\brief A class to manage an enhanced discretization for hybridizable
     discontinuous Galerkin methods (HDG)

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DISCRET_HDG_HPP
#define FOUR_C_LIB_DISCRET_HDG_HPP

#include "4C_config.hpp"

#include "4C_lib_discret_faces.hpp"
#include "4C_lib_utils_discret.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include <string>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class MapExtractor;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  enum class HDGAction
  {
    project_dirich_field  ///< project dirichlet field
  };

  class DiscretizationHDG : public DiscretizationFaces
  {
   public:
    /*!
    \brief Standard Constructor

    \param name (in): name of this discretization
    \param comm (in): An epetra comm object associated with this discretization
    */
    DiscretizationHDG(const std::string name, Teuchos::RCP<Epetra_Comm> comm);


    /*!
    \brief Complete construction of a discretization  (Filled()==true NOT prerequisite)

    After adding or deleting nodes or elements or redistributing them in parallel,
    or adding/deleting boundary conditions, this method has to be called to (re)construct
    pointer topologies.<br>
    It builds in this order:<br>
    Standard fill_complete of base class
    - row map of nodes
    - column map of nodes
    - row map of elements
    - column map of elements
    - pointers from elements to nodes
    - pointers from nodes to elements
    - assigns degrees of freedoms
    - map of element register classes
    - calls all element register initialize methods
    - build geometries of all Dirichlet and Neumann boundary conditions

    Additional features
    - build internal faces elements
    - build maps and pointers for internal faces

    \param assigndegreesoffreedom (in) : if true, resets existing dofsets and performs
                                         assigning of degrees of freedoms to nodes and
                                         elements.
    \param initelements (in) : if true, build element register classes and call Initialize()
                               on each type of finite element present
    \param doboundaryconditions (in) : if true, build geometry of boundary conditions
                                       present.

    \note In order to receive a fully functional discretization, this method must be called
          with all parameters set to true (the default). The parameters though can be
          used to turn off specific tasks to allow for more flexibility in the
          construction of a discretization, where it is known that this method will
          be called more than once.

    \note Sets Filled()=true
    */
    int fill_complete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true) override;

    /*!
     * this function has the same functionality as the function in the base class,
     * additionally, the degree of the elements is communicated, such that ghosted elements
     * also have full knowledge about the face degrees. This is necessary for discretizations
     * with nonuniform degree distributions and p-adaptivity
     *
     *  schoeder 06/14
     */
    void assign_global_i_ds(const Epetra_Comm& comm,
        const std::map<std::vector<int>, Teuchos::RCP<CORE::Elements::Element>>& elementmap,
        std::map<int, Teuchos::RCP<CORE::Elements::Element>>& finalelements) override;

  };  // class DiscretizationHDG

  namespace UTILS
  {
    /** \brief Specialized Dbc evaluation class for HDG discretizations
     *
     *  \author hiermeier \date 10/16 */
    class DbcHDG : public Dbc
    {
     public:
      /// constructor
      DbcHDG(){};

     protected:
      /** \brief Determine Dirichlet condition
       *
       *  \param cond    (in)  :  The condition object
       *  \param toggle  (out) :  Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids (out) :  Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       *  \author kronbichler \date 06/16 */
      void read_dirichlet_condition(const CORE::UTILS::FunctionManager& function_manager,
          const DRT::Discretization& discret, const CORE::Conditions::Condition& cond, double time,
          DbcInfo& info, const Teuchos::RCP<std::set<int>>* dbcgids,
          int hierarchical_order) const override;
      void read_dirichlet_condition(const CORE::UTILS::FunctionManager& function_manager,
          const DRT::DiscretizationFaces& discret, const CORE::Conditions::Condition& cond,
          double time, DbcInfo& info, const Teuchos::RCP<std::set<int>>* dbcgids,
          int hierarchical_order) const;

      /** \brief Determine Dirichlet condition at given time and apply its
       *         values to a system vector
       *
       *  \param cond            The condition object
       *  \param time            Evaluation time
       *  \param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
       *  \param systemvectord   First time derivative of DBCs
       *  \param systemvectordd  Second time derivative of DBCs
       *  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids         Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       *  \author kronbichler \date 02/08 */
      void do_dirichlet_condition(const CORE::UTILS::FunctionManager& function_manager,
          const DRT::Discretization& discret, const CORE::Conditions::Condition& cond, double time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const override;
      void do_dirichlet_condition(const CORE::UTILS::FunctionManager& function_manager,
          const DRT::DiscretizationFaces& discret, const CORE::Conditions::Condition& cond,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle) const;
    };  // class DbcHDG
  }     // namespace UTILS
}  // namespace DRT

/// << operator
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationHDG& dis);


FOUR_C_NAMESPACE_CLOSE

#endif
