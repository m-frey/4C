/*----------------------------------------------------------------------*/
/*! \file
\brief A condition of any kind

\level 1

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_MATERIAL_INPUT_BASE_HPP
#define FOUR_C_MATERIAL_INPUT_BASE_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_material_parameter_base.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace CORE::MAT::PAR
{
  // forward declaration
  class Parameter;

  /// Container for read-in materials
  ///
  /// This object stores the validated material parameters as
  /// IO::InputParameterContainer.
  class Material : public IO::InputParameterContainer
  {
   public:
    /// @name life span
    //@{

    /// Standard constructor
    Material(const int id,                        ///< unique material ID
        const CORE::Materials::MaterialType type  ///< type of material
    );

    /// Default constructor without information from the input lines.
    Material() = default;

    /// Set pointer to readily allocated 'quick access' material parameters
    ///
    /// This function is called by the material factory MAT::Factory.
    /// To circumvent more than this single major switch of material type to
    /// object, #params_ are allocated externally.
    inline void SetParameter(CORE::MAT::PAR::Parameter* matparam)
    {
      params_ = Teuchos::rcp(matparam);
    }

    //@}

    /// @name Query methods
    //@{

    /// Return material id
    [[nodiscard]] inline virtual int Id() const { return id_; }

    /// Return type of condition
    [[nodiscard]] inline virtual CORE::Materials::MaterialType Type() const { return type_; }

    /// Return quick accessible material parameter data
    ///
    /// These quick access parameters are stored in separate member #params_;
    /// whereas the originally read ones are stored in IO::InputParameterContainer base
    [[nodiscard]] inline CORE::MAT::PAR::Parameter* Parameter() const { return params_.get(); }

    //@}

    /// don't want = operator
    Material operator=(const Material& old) = delete;
    Material(const CORE::MAT::PAR::Material& old) = delete;

   protected:
    /// Unique ID of this material, no second material of same ID may exist
    int id_{};

    /// Type of this material
    CORE::Materials::MaterialType type_{};

    /// Unwrapped material data for 'quick' access
    Teuchos::RCP<CORE::MAT::PAR::Parameter> params_{};
  };
}  // namespace CORE::MAT::PAR

FOUR_C_NAMESPACE_CLOSE

#endif
