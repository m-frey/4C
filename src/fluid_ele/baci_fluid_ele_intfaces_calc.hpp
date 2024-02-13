/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of fluid internal faces elements


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_INTFACES_CALC_HPP
#define BACI_FLUID_ELE_INTFACES_CALC_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_inpar_xfem.hpp"
#include "baci_mat_material.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Condition;
  class Discretization;
  class DiscretizationFaces;

  namespace ELEMENTS
  {
    class FluidIntFace;
    class FluidEleParameter;
    class FluidEleParameterTimInt;
    class FluidEleParameterIntFace;

    /// Interface base class for FluidIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of FluidIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidIntFaceImpl.
     */
    class FluidIntFaceImplInterface
    {
     public:
      /// Empty constructor
      FluidIntFaceImplInterface() {}

      /// Empty destructor
      virtual ~FluidIntFaceImplInterface() = default;
      //! Assemble internal faces integrals using data from both parent elements
      virtual void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::FluidIntFace* intface,    ///< internal face element
          Teuchos::RCP<MAT::Material>& material,   ///< material associated with the faces
          std::vector<int>& nds_master,            ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,             ///< nodal dofset w.r.t. slave element
          const INPAR::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
          Teuchos::ParameterList& params,          ///< parameter list
          DRT::DiscretizationFaces& discretization,               ///< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int EvaluateInternalFaces(
          DRT::ELEMENTS::FluidIntFace* intface,   ///< internal face element
          Teuchos::RCP<MAT::Material>& material,  ///< material associated with the faces
          Teuchos::ParameterList& params,         ///< parameter list
          DRT::Discretization& discretization,    ///< discretization
          std::vector<int>& patchlm,              ///< patch local map
          std::vector<int>& lm_masterToPatch,     ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) = 0;


      /// Internal implementation class for FluidIntFace elements (the first object is created in
      /// DRT::ELEMENTS::FluidIntFace::Evaluate)
      static FluidIntFaceImplInterface* Impl(const DRT::Element* ele);
    };

    /// Internal FluidIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the FluidIntFace element.

      <h3>Purpose</h3>

      The FluidIntFace element will allocate exactly one object of this class
      for all FluidIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      \author gjb
      \date 08/08
      \author schott
      \date 04/12
    */
    template <CORE::FE::CellType distype>
    class FluidIntFaceImpl : public FluidIntFaceImplInterface
    {
     public:
      /// Singleton access method
      static FluidIntFaceImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Constructor
      FluidIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      void AssembleInternalFacesUsingNeighborData(
          DRT::ELEMENTS::FluidIntFace* intface,    ///< internal face element
          Teuchos::RCP<MAT::Material>& material,   ///< material associated with the faces
          std::vector<int>& nds_master,            ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,             ///< nodal dofset w.r.t. slave element
          const INPAR::XFEM::FaceType& face_type,  ///< which type of face std, ghost, ghost-penalty
          Teuchos::ParameterList& params,          ///< parameter list
          DRT::DiscretizationFaces& discretization,               ///< faces discretization
          Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
          ) override;

      //! Evaluate internal faces
      int EvaluateInternalFaces(DRT::ELEMENTS::FluidIntFace* intface,  ///< internal face element
          Teuchos::RCP<MAT::Material>& material,  ///< material associated with the faces
          Teuchos::ParameterList& params,         ///< parameter list
          DRT::Discretization& discretization,    ///< discretization
          std::vector<int>& patchlm,              ///< patch local map
          std::vector<int>& lm_masterToPatch,     ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,      ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,       ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) override;


     private:
      //! pointer to parameter list for time integration
      DRT::ELEMENTS::FluidEleParameterTimInt* fldparatimint_;
      //! pointer to parameter list for internal faces
      DRT::ELEMENTS::FluidEleParameterIntFace* fldpara_intface_;


    };  // end class FluidIntFaceImpl

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
