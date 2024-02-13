/*----------------------------------------------------------------------*/
/*! \file
\brief A condition of any kind

\level 1


*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef BACI_MAT_PAR_MATERIAL_HPP
#define BACI_MAT_PAR_MATERIAL_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_inpar_material.hpp"
#include "baci_lib_container.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* foward declarations */
namespace MAT
{
  namespace PAR
  {
    class Parameter;
  }
}  // namespace MAT

/*----------------------------------------------------------------------*/
/* declarations */
namespace MAT
{
  namespace PAR
  {
    class ParMaterialType : public CORE::COMM::ParObjectType
    {
     public:
      std::string Name() const override { return "ParMaterialType"; }

      static ParMaterialType& Instance() { return instance_; };

     private:
      static ParMaterialType instance_;
    };

    /// Container for read-in materials
    ///
    /// This object stores the validated material parameters as
    /// DRT::Container.
    ///
    /// \author bborn
    /// \date 02/09
    class Material : public DRT::Container
    {
     public:
      /// @name life span
      //@{

      /// standard constructor
      Material(const int id,                    ///< unique material ID
          const INPAR::MAT::MaterialType type,  ///< type of material
          const std::string name                ///< name of material
      );


      /// Empty Constructor with type condition_none
      Material();


      /// Copy Constructor
      ///
      /// Makes a deep copy of a condition
      Material(const MAT::PAR::Material& old);


      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of this file.
      int UniqueParObjectId() const override
      {
        return ParMaterialType::Instance().UniqueParObjectId();
      }

      /// Pack this class so it can be communicated
      ///
      /// \ref Pack and \ref Unpack are used to communicate this class
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this class
      void Unpack(const std::vector<char>& data) override;

      /// Set pointer to readily allocated 'quick access' material parameters
      ///
      /// This function is called by the material factory MAT::Material::Factory.
      /// To circumvent more than this single major switch of material type to
      /// object, #params_ are allocated externally.
      inline void SetParameter(MAT::PAR::Parameter* matparam)  ///< the pointer
      {
        params_ = Teuchos::rcp(matparam);
      }

      //@}

      /// @name Query methods
      //@{

      /// Return material id
      inline virtual int Id() const { return id_; }

      /// Return material name
      inline virtual std::string Name() const { return name_; }

      /// Print this Condition (std::ostream << is also implemented for DRT::Condition)
      void Print(std::ostream& os) const override;

      /// Return type of condition
      inline virtual INPAR::MAT::MaterialType Type() const { return type_; }

      /// Return communicator
      inline Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

      /// Return quick accessible material parameter data
      ///
      /// These quick access parameters are stored in separate member #params_;
      /// whereas the originally read ones are stored in DRT::Container base
      inline MAT::PAR::Parameter* Parameter() const { return params_.get(); }

      //@}

     protected:
      /// don't want = operator
      Material operator=(const Material& old);

      /// Unique ID of this material, no second material of same ID may exist
      int id_;

      /// Type of this condition
      INPAR::MAT::MaterialType type_;

      /// Name
      std::string name_;

      /// A communicator
      Teuchos::RCP<Epetra_Comm> comm_;

      /// Unwrapped material data for 'quick' access
      Teuchos::RCP<MAT::PAR::Parameter> params_;

    };  // class Material

  }  // namespace PAR

}  // namespace MAT


/// out stream operator
std::ostream& operator<<(std::ostream& os, const MAT::PAR::Material& cond);


BACI_NAMESPACE_CLOSE

#endif  // MAT_PAR_MATERIAL_H
