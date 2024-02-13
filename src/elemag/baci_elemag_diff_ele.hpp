/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the electromagnetic diffusion element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

<pre>
\level 2

</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ELEMAG_DIFF_ELE_HPP
#define BACI_ELEMAG_DIFF_ELE_HPP

#include "baci_config.hpp"

#include "baci_elemag_ele.hpp"
#include "baci_linalg_serialdensematrix.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class ElemagDiffType : public ElemagType
    {
     public:
      /// Type name
      std::string Name() const override { return "ElemagDiffType"; }

      // Instance
      static ElemagDiffType& Instance();
      /// Create
      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;
      /// Create
      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;
      /// Create
      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      /// Nodal block information
      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      /// Null space computation
      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      /// Element definition
      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      /// Instance
      static ElemagDiffType instance_;
    };


    /*!
    \brief electromagnetic diffusion element
    */
    class ElemagDiff : public Elemag
    {
     public:
      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      ElemagDiff(int id,  ///< A unique global id
          int owner       ///< Owner
      );
      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      ElemagDiff(const ElemagDiff& old);

      /*!
      \brief Deep copy this instance and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and
      the parent_slave element
      */
      Teuchos::RCP<DRT::Element> CreateFaceElement(
          DRT::Element* parent_slave,  //!< parent slave element
          int nnode,                   //!< number of surface nodes
          const int* nodeids,          //!< node ids of surface element
          DRT::Node** nodes,           //!< nodes of surface element
          const int lsurface_master,   //!< local surface number w.r.t master parent element
          const int lsurface_slave,    //!< local surface number w.r.t slave parent element
          const std::vector<int>& localtrafomap  //! local trafo map
          ) override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return ElemagDiffType::Instance().UniqueParObjectId();
      }

      //@}

      //! @name Geometry related methods

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ElemagDiffType::Instance(); }

      /// Element location data
      LocationData lm_;

     private:
      // don't want = operator
      ElemagDiff& operator=(const ElemagDiff& old);

    };  // class ElemagDiff

    /// class ElemagDiffBoundaryType
    class ElemagDiffBoundaryType : public ElemagBoundaryType
    {
     public:
      /// Type name
      std::string Name() const override { return "ElemagDiffBoundaryType"; }
      // Instance
      static ElemagDiffBoundaryType& Instance();
      // Create
      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

     private:
      /// Instance
      static ElemagDiffBoundaryType instance_;
    };

    /// class ElemagDiffBoundary
    class ElemagDiffBoundary : public ElemagBoundary
    {
     public:
      //! @name Constructors and destructors and related methods

      //! number of space dimensions
      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent elemag element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      ElemagDiffBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::ElemagDiff* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      ElemagDiffBoundary(const ElemagDiffBoundary& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return ElemagDiffBoundaryType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override
      {
        return ParentElement()->NumDofPerNode(node);
      }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      /// Return the instance of the element type
      DRT::ElementType& ElementType() const override { return ElemagDiffBoundaryType::Instance(); }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //@}

      /*!
      \brief Return the location vector of this element

      The method computes degrees of freedom this element adresses.
      Degree of freedom ordering is as follows:<br>
      First all degrees of freedom of adjacent nodes are numbered in
      local nodal order, then the element internal degrees of freedom are
      given if present.<br>
      If a derived element has to use a different ordering scheme,
      it is welcome to overload this method as the assembly routines actually
      don't care as long as matrices and vectors evaluated by the element
      match the ordering, which is implicitly assumed.<br>
      Length of the output vector matches number of degrees of freedom
      exactly.<br>
      This version is intended to fill the LocationArray with the dofs
      the element will assemble into. In the standard case these dofs are
      the dofs of the element itself. For some special conditions (e.g.
      the weak dirichlet boundary condtion) a surface element will assemble
      into the dofs of a volume element.<br>

      \note The degrees of freedom returned are not neccessarily only nodal dofs.
            Depending on the element implementation, output might also include
            element dofs.

      \param dis (in)      : the discretization this element belongs to
      \param la (out)      : location data for all dofsets of the discretization
      \param doDirichlet (in): whether to get the Dirichlet flags
      \param condstring (in): Name of condition to be evaluated
      \param condstring (in):  List of parameters for use at element level
      */
      void LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
          const std::string& condstring, Teuchos::ParameterList& params) const override;

     private:
      // don't want = operator
      ElemagDiffBoundary& operator=(const ElemagDiffBoundary& old);

    };  // class ElemagDiffBoundary

    /// class ElemagDiffIntFaceType
    class ElemagDiffIntFaceType : public DRT::ElementType
    {
     public:
      /// Name of the element type
      std::string Name() const override { return "ElemagDiffIntFaceType"; }

      /// Instance
      static ElemagDiffIntFaceType& Instance();

      /// Create
      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      /// Nodal block information
      void NodalBlockInformation(
          Element* dwele, int& numdf, int& dimns, int& nv, int& np) override{};

      /// Null space
      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented");
        return nullspace;
      };

     private:
      /// instance of the class
      static ElemagDiffIntFaceType instance_;
    };

    /// class ElemagDiffIntFace
    class ElemagDiffIntFace : public ElemagIntFace
    {
     public:
      //! @name Constructors and destructors and related methods

      //! number of space dimensions
      /*!
      \brief Standard Constructor

      \param id: A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param master_parent: The master parent elemag element of this surface
      \param slave_parent: The slave parent elemag element of this surface
      \param lsurface_master: the local surface number of this surface w.r.t. the master parent
      element \param lsurface_slave: the local surface number of this surface w.r.t. the slave
      parent element \param localtrafomap: transformation map between the local coordinate systems
      of the face w.r.t the master parent element's face's coordinate system and the slave element's
      face's coordinate system
      */
      ElemagDiffIntFace(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::ElemagDiff* parent_master, DRT::ELEMENTS::ElemagDiff* parent_slave,
          const int lsurface_master, const int lsurface_slave,
          const std::vector<int> localtrafomap);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      ElemagDiffIntFace(const ElemagDiffIntFace& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return ElemagDiffIntFaceType::Instance().UniqueParObjectId();
      }

      //@}

      //! @name Acess methods

      /*!
      \brief create the location vector for patch of master and slave element

      \note All dofs shared by master and slave element are contained only once. Dofs from interface
      nodes are also included.
      */
      void PatchLocationVector(DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& nds_master,        ///< nodal dofset w.r.t master parent element
          std::vector<int>& nds_slave,         ///< nodal dofset w.r.t slave parent element
          std::vector<int>& patchlm,           ///< local map for gdof ids for patch of elements
          std::vector<int>& master_lm,         ///< local map for gdof ids for master element
          std::vector<int>& slave_lm,          ///< local map for gdof ids for slave element
          std::vector<int>& face_lm,           ///< local map for gdof ids for face element
          std::vector<int>& lm_masterToPatch,  ///< local map between lm_master and lm_patch
          std::vector<int>& lm_slaveToPatch,   ///< local map between lm_slave and lm_patch
          std::vector<int>& lm_faceToPatch,    ///< local map between lm_face and lm_patch
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch  ///< local map between slave nodes and nodes in patch
      );
      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ElemagDiffIntFaceType::Instance(); }

      //@}

      /*!
      \brief return the master parent elemag element
      */
      DRT::ELEMENTS::ElemagDiff* ParentMasterElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentMasterElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::ElemagDiff*>(parent) != nullptr,
            "Master element is no elemag_diff element");
        return static_cast<DRT::ELEMENTS::ElemagDiff*>(parent);
      }

      /*!
      \brief return the slave parent elemag element
      */
      DRT::ELEMENTS::ElemagDiff* ParentSlaveElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentSlaveElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::ElemagDiff*>(parent) != nullptr,
            "Slave element is no elemag_diff element");
        return static_cast<DRT::ELEMENTS::ElemagDiff*>(parent);
      }

      //@}

     private:
      // don't want = operator
      ElemagDiffIntFace& operator=(const ElemagDiffIntFace& old);

    };  // class ElemagDiffIntFace

  }  // namespace ELEMENTS
}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_DIFF_ELE_H
