/*-----------------------------------------------------------*/
/*! \file

\brief Methods to apply rotationally symmetric periodic boundary
       conditions for fluid problems on element level


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ROTSYM_PERIODICBC_HPP
#define FOUR_C_FLUID_ROTSYM_PERIODICBC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_rotsym_periodicbc_utils.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  /// numdofpernode = number of dofs per node for fluid problem
  template <Core::FE::CellType distype, int numdofpernode,
      Discret::ELEMENTS::Fluid::EnrichmentType enrtype>

  /*!
  \brief  This class manages local transformations(rotation) of velocity fields

          o used for rotationally symmetric boundary conditions
          o used for rotationally symmetric boundary conditions

  \author gjb

   */
  class RotationallySymmetricPeriodicBC
  {
   public:
    /// number of nodes for this element type including virtual nodes
    static constexpr int elenumnode =
        Discret::ELEMENTS::MultipleNumNode<enrtype>::multipleNode * Core::FE::num_nodes<distype>;
    /// number of nodes for this element type (only real nodes)
    static constexpr int elenumnodereal = Core::FE::num_nodes<distype>;

    /// standard constructor
    explicit RotationallySymmetricPeriodicBC() : rotangle_(0.0), rotmat_(true) { return; };

    /// empty destructor
    virtual ~RotationallySymmetricPeriodicBC() = default;


    /// is any rotation needed on the element level?
    bool HasRotSymmPBC()
    {
      if (slavenodelids_.empty())
      {
        return (false);
      }
      else
      {
        return (true);
      }
    }

    /// prepare the class for this element
    void setup(Core::Elements::Element* ele)
    {
      // clean everything
      rotangle_ = 0.0;
      slavenodelids_.clear();
      rotmat_.clear();

      Core::Nodes::Node** nodes = ele->Nodes();
      slavenodelids_.reserve(elenumnodereal);

      for (int inode = 0; inode < elenumnodereal; inode++)
      {
        if (FLD::IsSlaveNodeOfRotSymPBC(nodes[inode], rotangle_))
        {
          // store the lid of this node
          slavenodelids_.push_back(inode);
        }
      }

      // yes, we prepare the rotation matrix (identity matrix when there are no slave node)
      for (int i = 0; i < elenumnode * numdofpernode; i++)
      {
        rotmat_(i, i) = 1.0;
      }
      const double c = cos(rotangle_);
      const double s = sin(rotangle_);
      for (unsigned int slnode = 0; slnode < slavenodelids_.size(); slnode++)
      {
        // get local id of current slave node
        int inode = slavenodelids_[slnode];
        // velocity x- and y-components have to be rotated
        rotmat_(inode * numdofpernode, inode * numdofpernode) = c;
        rotmat_(inode * numdofpernode, inode * numdofpernode + 1) = s;
        rotmat_(inode * numdofpernode + 1, inode * numdofpernode) = -s;
        rotmat_(inode * numdofpernode + 1, inode * numdofpernode + 1) = c;
        // corresponding z- component and pressure entries remain unchanged!
      }

      return;
    }

    /// rotate velocity vector used in element routine if necessary
    void rotate_my_values_if_necessary(std::vector<double>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          int i = slavenodelids_[slnodes];
          double xvalue = myvalues[0 + (i * numdofpernode)];
          double yvalue = myvalues[1 + (i * numdofpernode)];
          myvalues[0 + (i * numdofpernode)] = c * xvalue - s * yvalue;
          myvalues[1 + (i * numdofpernode)] = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate velocity vector used in element routine if necessary
    void rotate_my_values_if_necessary(
        Core::LinAlg::Matrix<numdofpernode - 1, elenumnode>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          const int i = slavenodelids_[slnodes];
          const double xvalue = myvalues(0, i);
          const double yvalue = myvalues(1, i);
          myvalues(0, i) = c * xvalue - s * yvalue;
          myvalues(1, i) = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate velocity vector used in element routine if necessary
    template <int rows, int cols>
    void rotate_my_values_if_necessary(Core::LinAlg::Matrix<rows, cols>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          const int i = slavenodelids_[slnodes];
          const double xvalue = myvalues(0, i);
          const double yvalue = myvalues(1, i);
          myvalues(0, i) = c * xvalue - s * yvalue;
          myvalues(1, i) = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate element matrix and vectors if necessary (first version)
    void rotate_matand_vec_if_necessary(
        Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat1,
        Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat2,
        Core::LinAlg::Matrix<numdofpernode * elenumnode, 1>& elevec1)
    {
      if (HasRotSymmPBC())
      {
        // execute the rotation
        /*
        K_rot = Q* K * Q^T
        b_rot = Q * b
        with Q as defined in the setup() function
         */
        if (elemat1.IsInitialized())  // do not try to access an uninitialized matrix!
        {
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elematold1(
              elemat1);
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(
              true);
          tempmatrix.MultiplyNT(elematold1, rotmat_);
          elemat1.Multiply(rotmat_, tempmatrix);
        }
        if (elemat2.IsInitialized())
        {
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elematold2(
              elemat2);
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(
              true);
          tempmatrix.MultiplyNT(elematold2, rotmat_);
          elemat2.Multiply(rotmat_, tempmatrix);
        }
        if (elevec1.IsInitialized())
        {
          Core::LinAlg::Matrix<numdofpernode * elenumnode, 1> elevec1old(elevec1);
          elevec1.Multiply(rotmat_, elevec1old);
        }
      }
      return;
    }

    /// rotate element matrix and vectors if necessary (second version)
    void rotate_matand_vec_if_necessary(
        Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat1,
        Core::LinAlg::Matrix<numdofpernode * elenumnode, 1>& elevec1)
    {
      if (HasRotSymmPBC())
      {
        // execute the rotation
        /*
        K_rot = Q* K * Q^T
        b_rot = Q * b
        with Q as defined in the setup() function
         */
        if (elemat1.IsInitialized())  // do not try to access an uninitialized matrix!
        {
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elemat1old(
              elemat1);
          Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(
              true);
          tempmatrix.MultiplyNT(elemat1old, rotmat_);
          elemat1.Multiply(rotmat_, tempmatrix);
        }

        if (elevec1.IsInitialized())
        {
          Core::LinAlg::Matrix<numdofpernode * elenumnode, 1> elevec1old(elevec1);
          elevec1.Multiply(rotmat_, elevec1old);
        }
      }
      return;
    }


   protected:
   private:
    //! vector of local slave node ids of the applied periodic surface boundary conditions
    std::vector<int> slavenodelids_;

    //! angle of rotation (slave plane <--> master plane)
    double rotangle_;

    //! rotation matrix
    Core::LinAlg::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> rotmat_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
