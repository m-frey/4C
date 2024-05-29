/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet4 Element
\level 3
*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_TET4_HPP
#define FOUR_C_SO3_TET4_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_so3_base.hpp"


#define NUMNOD_SOTET4 4       ///< number of nodes
#define NODDOF_SOTET4 3       ///< number of dofs per node
#define NUMDOF_SOTET4 12      ///< total dofs per element
#define NUMGPT_SOTET4 1       ///< total gauss points per element  /****/
#define NUMDIM_SOTET4 3       ///< number of dimensions/****/
#define NUMCOORD_SOTET4 4     ///< number of shape function cooordinates (ksi1-ksi4)
#define NUMNOD_SOTET4_FACE 3  ///< number of nodes on a TET4 face (which is a TRI3)
#define NUMGPT_SOTET4_FACE 1  ///< number of GP on a TET4 face (which is a TRI3)

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class PreStress;

    class SoTet4Type : public CORE::Elements::ElementType
    {
     public:
      std::string Name() const override { return "So_tet4Type"; }

      static SoTet4Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(
          std::string eletype, std::string eledistype, int id, int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(int id, int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void nodal_block_information(
          CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4Type instance_;

      std::string get_element_type_string() const { return "SOLIDT4"; }
    };

    /*!
    \brief A C++ version of the 4-node tet solid element

    */

    class SoTet4 : public SoBase
    {
     public:
      //! @name Friends
      friend class SoTet4Type;


      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      SoTet4(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoTet4(const SoTet4& old);

      // don't want = operator
      SoTet4& operator=(const SoTet4& old) = delete;

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      CORE::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return number of volumes of this element
      */
      int NumVolume() const override { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return 4; }

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return 6; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Surfaces() override;

      virtual std::vector<double> element_center_refe_coords();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoTet4Type::Instance().UniqueParObjectId(); }

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
             (implements pure virtual CORE::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual CORE::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      CORE::Elements::ElementType& ElementType() const override { return SoTet4Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      /*!
      \brief Query names of element data to be visualized using BINIO

      The element fills the provided map with key names of
      visualization data the element wants to visualize AT THE CENTER
      of the element geometry. The values is supposed to be dimension of the
      data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      or 9 (nonsym. tensor)

      Example:
      \code
        // Name of data is 'Owner', dimension is 1 (scalar value)
        names.insert(std::pair<std::string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<std::string,int>("StressesXYZ",6));
      \endcode

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      void VisNames(std::map<std::string, int>& names) override;

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if element recognizes the name
      */
      bool VisData(const std::string& name, std::vector<double>& data) override;

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;


      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate so_tet4 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int evaluate_neumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          CORE::Conditions::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      double EvaluationCost() override
      {
        if (Material()->MaterialType() == CORE::Materials::m_struct_multiscale)
          return 25000.0;
        else
          return 10.0;
      }

      void get_cauchy_n_dir_and_derivatives_at_xi(const CORE::LINALG::Matrix<3, 1>& xi,
          const std::vector<double>& disp, const CORE::LINALG::Matrix<3, 1>& n,
          const CORE::LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
          CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd2,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dn,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
          CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dT,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
          double* d_cauchyndir_dc) override;
      //@}

     protected:
      //! action parameters recognized by so_tet4
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        struct_calc_store_istep,
        struct_calc_recover_istep,
        calc_struct_update_istep,
        calc_struct_reset_istep,     //!< reset elementwise internal variables
                                     //!< during iteration to last converged state
        calc_struct_reset_all,       //!< reset elementwise internal variables
                                     //!< to state in the beginning of the computation
        calc_global_gpstresses_map,  //! basically calc_struct_stress but with assembly of global
                                     //! gpstresses map
        prestress_update,
        calc_struct_energy,
        calc_struct_output_E,
        multi_calc_dens,
        multi_readrestart
      };

      //! volume of the element
      double V_;

      CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void update_jacobian_mapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      /// compute defgrd ypein all gp for given disp
      void def_gradient(const std::vector<double>& disp, CORE::LINALG::SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      /*!
       * \brief Compute the deformation gradient
       *
       * \param defgrd Deformation gradient
       * \param xdisp Displacement vectir for each node (3x4)
       * \param gp Gauss point
       */
      void compute_deformation_gradient(CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4>& defgrd,
          const CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4>& xdisp, int gp) const;


      // internal calculation methods

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void init_jacobian_mapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  // location matrix
          std::vector<double>& disp,                   // current displacements
          std::vector<double>* vel,                    // current velocities
          std::vector<double>* acc,                    // current accelerations
          std::vector<double>& residual,               // current residual displ
          std::vector<double>& dispmat,                // current material displacements
          CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>*
              stiffmatrix,  // element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* massmatrix,  // element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* force,       // element internal force vector
          CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* forceinert,  // element inertial force vector
          CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* force_str,   // element structural force vector
          CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
          CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>*
              eleplstrain,                   // plastic strains at GP
          Teuchos::ParameterList& params,    // algorithmic parameters e.g. time
          INPAR::STR::StressType iostress,   // stress output option
          INPAR::STR::StrainType iostrain,   // strain output option
          INPAR::STR::StrainType ioplstrain  // plastic strain output option
      );

      //! lump mass matrix (bborn 07/08)
      void so_tet4_lumpmass(
          CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* emass);  //!< element mass matrix

      //! remodeling for fibers at the end of time step (st 01/10)
      void so_tet4_remodel(std::vector<int>& lm,          // location matrix
          std::vector<double>& disp,                      // current displacements
          Teuchos::ParameterList& params,                 // algorithmic parameters e.g. time
          const Teuchos::RCP<CORE::MAT::Material>& mat);  // material

      //! Evaluate Tet4 Shapefcts at 1 gausspoint to keep them static
      std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> so_tet4_1gp_shapefcts();
      //! Evaluate Tet4 Derivs at 1 gausspoint to keep them static
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> so_tet4_1gp_derivs();
      //! Evaluate Tet4 Weights at 1 gausspoint to keep them static
      std::vector<double> so_tet4_1gp_weights();

      //! Evaluate Tet4 Shapefcts at 4 gausspoints to keep them static
      std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> so_tet4_4gp_shapefcts();
      //! Evaluate Tet4 Derivs at 4 gausspoints to keep them static
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> so_tet4_4gp_derivs();
      //! Evaluate Tet4 Weights at 4 gausspoints to keep them static
      std::vector<double> so_tet4_4gp_weights();

      //! @name Multi-scale related stuff

      /// Determine a homogenized material density for multi-scale analyses by averaging over the
      /// initial volume
      void sotet4_homog(Teuchos::ParameterList& params);

      /// Read restart on the microscale
      void sotet4_read_restart_multi();

      //@}

      /*!
       * Executes the post setup call for all materials. This method will be called once per element
       * at the first Evaluate call.
       *
       * @param params ParameterList to be passed to the materials
       */
      void material_post_setup(Teuchos::ParameterList& params) override;

     private:
      std::string get_element_type_string() const { return "SOLIDT4"; }
    };  // class So_tet4


  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
