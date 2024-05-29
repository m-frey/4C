/*-----------------------------------------------------------------------*/
/*! \file
\brief header for projector functions

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_PROJECTOR_HPP
#define FOUR_C_MORTAR_PROJECTOR_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_pairedvector.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MORTAR
{
  // forward declarations
  class Node;
  class Element;

  /*!
  \brief A class to perform projections of nodes onto opposing elements

  */

  class Projector
  {
   public:
    Projector(){};

    virtual ~Projector() = default;
    //! @name Access methods

    /// 1. Internal implementation class
    static Projector* Impl(MORTAR::Element& ele);

    /// 2. Internal implementation class
    static Projector* Impl(MORTAR::Element& sele, MORTAR::Element& mele);

    //! @name virtual functions
    virtual bool ProjectNodalNormal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) = 0;

    virtual bool project_element_normal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) = 0;

    virtual bool ProjectGaussPoint2D(
        MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele, double* xi) = 0;

    virtual bool ProjectGaussPoint3D(MORTAR::Element& gpele, const double* gpeta,
        MORTAR::Element& ele, double* xi, double& par) = 0;

    virtual bool project_gauss_point_auxn3_d(const double* globgp, const double* auxn,
        MORTAR::Element& ele, double* xi, double& par) = 0;

    virtual bool project_s_node_by_m_normal(
        MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist) = 0;


    virtual bool project_s_node_by_m_nodal_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) = 0;

    virtual bool project_s_node_by_m_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) = 0;
  };  // class Projector

  /*!
  \brief A class to perform projections of nodes onto opposing elements

  */
  template <CORE::FE::CellType distype>
  class ProjectorCalc : public Projector
  {
   public:
    // constructor
    ProjectorCalc();

    /// Singleton access method
    static ProjectorCalc<distype>* Instance(
        CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

    //! n_: number of element nodes
    static constexpr int n_ = CORE::FE::num_nodes<distype>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = CORE::FE::dim<distype> + 1;

    //! @name 2D and 3D projection methods

    /*!
    \brief Project a node onto an element along the node's normal

    Used to project a slave side node onto an element of the master side

    This method will compute the coordinates of a projection of a node in
    the local coordinate system of an element. The projection point will
    not neccesarily fall inside the element. However, if the projection
    point is far outside the segment's boundaries, problems with the
    internal nonlinear iteration might occur and a warning is issued when
    convergence can not be achieved in a limited number of iterations.

    \param node (in): Slave node to project
    \param ele (in) : Master element to project on
    \param xi (out) : Local coordinates of projection on element

    */
    bool ProjectNodalNormal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) override;

    /*!
    \brief Project a node onto an element along the interpolated
           outward normal field of the element

    Used to project a master side node onto an element of the slave side

    This method will compute the coordinates of a projection of a node in
    the local coordinate system of an element. The projection point will
    not neccesarily fall inside the element. However, if the projection
    point is far outside the segment's boundaries, problems with the
    internal nonlinear iteration might occur and a warning is issued when
    convergence can not be achieved in a limited number of iterations.

    \param node (in): Master node to project
    \param ele (in) : Slave element to project on
    \param xi (out) : Local coordinates of projection on element

    */
    bool project_element_normal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) override;

    /*!
    \brief Project a Gauss point onto an element along GP normal

    Used to project a slave side GP onto an element of the master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param gpele (in): Slave element containing GP to project
    \param gpeta (in): Local coordinates of GP on gpele
    \param ele (in)  : Master element to project on
    \param xi (out)  : Local coordinates of projection on master element

    */
    bool ProjectGaussPoint2D(
        MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele, double* xi) override
    {
      FOUR_C_THROW("Called ele-based projection for segment-based integration!!!");
      return false;
    };

    /*!
    \brief Project a Gauss point onto an element along GP normal (3D)

    Used to project a slave side GP onto an element of the master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param gpele (in): Slave element containing GP to project
    \param gpeta (in): Local coordinates of GP on gpele
    \param ele (in)  : Master element to project on
    \param xi (out)  : Local coordinates of projection on master element
    \param par (out ): Projection parameter alpha

    */
    bool ProjectGaussPoint3D(MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele,
        double* xi, double& par) override
    {
      FOUR_C_THROW("Called ele-based projection for segment-based integration!!!");
      return false;
    };

    /*!
    \brief Project a Gauss point onto an element along AuxPlane normal (3D)

    Used to project an AuxPlane GP onto an element of the slave or master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param globgp(in): Gauss point to project, given in global coords
    \param auxn(in)  : Normal of AuxPlane along which to project
    \param ele (in)  : Slave / master element to project on
    \param xi (out)  : Local coordinates of projection on element
    \param par (out ): Projection parameter alpha

    */
    bool project_gauss_point_auxn3_d(const double* globgp, const double* auxn, MORTAR::Element& ele,
        double* xi, double& par) override;

    // TODO explanation
    bool project_s_node_by_m_normal(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist) override;

    // TODO explanation
    bool project_s_node_by_m_nodal_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) override;

    // TODO explanation
    bool project_s_node_by_m_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) override;

   protected:
    bool project_s_node_by_m_normal3_d(
        MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist);

    bool project_s_node_by_m_normal3_d_lin(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin);

    bool project_s_node_by_m_normal2_d(
        MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist);

    bool project_s_node_by_m_nodal_normal2_d_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin);

    bool project_s_node_by_m_nodal_normal3_d_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin);

    bool project_s_node_by_m_normal2_d_lin(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin);

    /*!
    \brief Evaluate F for nodal normal projection

    */
    double evaluate_f_nodal_normal(MORTAR::Node& node, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate GradF for nodal normal projection

    */
    double evaluate_grad_f_nodal_normal(
        MORTAR::Node& node, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate F for element normal projection

    */
    double evaluate_f_element_normal(MORTAR::Node& node, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate GradF for element normal projection

    */
    double evaluate_grad_f_element_normal(
        MORTAR::Node& node, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate F for AuxPlane Gauss point projection (3D)

    */
    bool evaluate_f_gauss_point_auxn3_d(double* f, const double* globgp, const double* auxn,
        MORTAR::Element& ele, const double* eta, const double& alpha);

    /*!
    \brief Evaluate GradF for AuxPlane Gauss point projection (3D)

    */
    bool evaluate_grad_f_gauss_point_auxn3_d(CORE::LINALG::Matrix<3, 3>& fgrad,
        const double* globgp, const double* auxn, MORTAR::Element& ele, const double* eta,
        const double& alpha);
  };

  /*!
  \brief A class to perform element based projections of nodes onto opposing elements

  */
  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class ProjectorCalcEleBased : public Projector
  {
   public:
    // constructor
    ProjectorCalcEleBased();

    /// Singleton access method
    static ProjectorCalcEleBased<distypeS, distypeM>* Instance(
        CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

    //! ns_: number of slave element nodes
    static constexpr int ns_ = CORE::FE::num_nodes<distypeS>;

    //! nm_: number of master element nodes
    static constexpr int nm_ = CORE::FE::num_nodes<distypeM>;

    //! number of space dimensions ("+1" due to considering only interface elements)
    static constexpr int ndim_ = CORE::FE::dim<distypeS> + 1;

    //! @name 2D and 3D projection methods

    /*!
    \brief Project a node onto an element along the node's normal

    Used to project a slave side node onto an element of the master side

    This method will compute the coordinates of a projection of a node in
    the local coordinate system of an element. The projection point will
    not neccesarily fall inside the element. However, if the projection
    point is far outside the segment's boundaries, problems with the
    internal nonlinear iteration might occur and a warning is issued when
    convergence can not be achieved in a limited number of iterations.

    \param node (in): Slave node to project
    \param ele (in) : Master element to project on
    \param xi (out) : Local coordinates of projection on element

    */
    bool ProjectNodalNormal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) override
    {
      FOUR_C_THROW("Called segment-based projection for element-based integration!!!");
      return false;
    };

    /*!
    \brief Project a node onto an element along the interpolated
           outward normal field of the element

    Used to project a master side node onto an element of the slave side

    This method will compute the coordinates of a projection of a node in
    the local coordinate system of an element. The projection point will
    not neccesarily fall inside the element. However, if the projection
    point is far outside the segment's boundaries, problems with the
    internal nonlinear iteration might occur and a warning is issued when
    convergence can not be achieved in a limited number of iterations.

    \param node (in): Master node to project
    \param ele (in) : Slave element to project on
    \param xi (out) : Local coordinates of projection on element

    */
    bool project_element_normal(MORTAR::Node& node, MORTAR::Element& ele, double* xi) override
    {
      FOUR_C_THROW("Called segment-based projection for element-based integration!!!");
      return false;
    };

    /*!
    \brief Project a Gauss point onto an element along GP normal

    Used to project a slave side GP onto an element of the master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param gpele (in): Slave element containing GP to project
    \param gpeta (in): Local coordinates of GP on gpele
    \param ele (in)  : Master element to project on
    \param xi (out)  : Local coordinates of projection on master element

    */
    bool ProjectGaussPoint2D(
        MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele, double* xi) override;

    /*!
    \brief Project a Gauss point onto an element along GP normal (3D)

    Used to project a slave side GP onto an element of the master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param gpele (in): Slave element containing GP to project
    \param gpeta (in): Local coordinates of GP on gpele
    \param ele (in)  : Master element to project on
    \param xi (out)  : Local coordinates of projection on master element
    \param par (out ): Projection parameter alpha

    */
    bool ProjectGaussPoint3D(MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele,
        double* xi, double& par) override;

    /*!
    \brief Project a Gauss point onto an element along AuxPlane normal (3D)

    Used to project an AuxPlane GP onto an element of the slave or master side

    This method will compute the coordinates of a projection of a Gauss
    point in the local coordinate system of an element.

    \param globgp(in): Gauss point to project, given in global coords
    \param auxn(in)  : Normal of AuxPlane along which to project
    \param ele (in)  : Slave / master element to project on
    \param xi (out)  : Local coordinates of projection on element
    \param par (out ): Projection parameter alpha

    */
    bool project_gauss_point_auxn3_d(const double* globgp, const double* auxn, MORTAR::Element& ele,
        double* xi, double& par) override
    {
      FOUR_C_THROW("Called Aux.-plane projection for element-based integration!!!");
      return false;
    };

    bool project_s_node_by_m_normal(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist) override
    {
      FOUR_C_THROW("ERROR");
      return false;
    };

    bool project_s_node_by_m_nodal_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele,
        double* xi, double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) override
    {
      FOUR_C_THROW("ERROR");
      return false;
    };

    bool project_s_node_by_m_normal_lin(MORTAR::Node& snode, MORTAR::Element& mele, double* xi,
        double* normal, double& dist,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin) override
    {
      FOUR_C_THROW("ERROR");
      return false;
    };

   protected:
    /*!
    \brief Check intersection of projection normal with warped element to project on.

    Create aux-plane for every ele-node with adjacent element edges.
    If all intersection points of the projection normal and the aux-planes
    are located inside the corresponding ele-edges, then the gp lie on the element.
    --> return false for projection check

    \param ele(in)    : element to project on
    \param ngp(in)    : gp-normal to project along
    \param globgp(in) : global gp coordinates

    */
    bool check_projection4_auxplane(MORTAR::Element& ele, double* ngp, double* globgp);

    /*!
    \brief Evaluate F for Gauss point projection

    */
    double evaluate_f_gauss_point2_d(
        const double* gpx, const double* gpn, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate GradF for Gauss point projection

    */
    double evaluate_grad_f_gauss_point2_d(
        const double* gpn, MORTAR::Element& ele, const double* eta);

    /*!
    \brief Evaluate F for Gauss point projection (3D)

    */
    bool evaluate_f_gauss_point3_d(double* f, const double* gpx, const double* gpn,
        MORTAR::Element& ele, const double* eta, const double& alpha);

    /*!
    \brief Evaluate GradF for Gauss point projection (3D)

    */
    bool evaluate_grad_f_gauss_point3_d(CORE::LINALG::Matrix<3, 3>& fgrad, const double* gpx,
        const double* gpn, MORTAR::Element& ele, const double* eta, const double& alpha);
  };
}  // namespace MORTAR


FOUR_C_NAMESPACE_CLOSE

#endif
