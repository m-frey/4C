/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_SSI_STR_MODEL_EVALUATOR_PARTITIONED_HPP
#define FOUR_C_SSI_STR_MODEL_EVALUATOR_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_ssi_str_model_evaluator_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Structure;
}  // namespace Adapter

namespace SSI
{
  class SSIPart;
}

namespace Solid
{
  namespace ModelEvaluator
  {
    class PartitionedSSI : public BaseSSI
    {
     public:
      //! constructor
      PartitionedSSI(const Teuchos::RCP<const SSI::SSIPart>
              ssi_part  //!< partitioned algorithm for scalar-structure interaction
      );

      void setup() override;

      //! @name Functions which are derived from the base generic class
      //! @{
      [[nodiscard]] Inpar::Solid::ModelType type() const override
      {
        return Inpar::Solid::model_partitioned_coupling;
      }

      bool assemble_force(Core::LinAlg::Vector& f, const double& timefac_np) const override;

      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      void determine_stress_strain() override{};

      void run_pre_compute_x(const Core::LinAlg::Vector& xold, Core::LinAlg::Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override;
      //! @}

     private:
      //! partitioned algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSIPart> ssi_part_;
    };  // class PartitionedSSI

  }  // namespace ModelEvaluator
}  // namespace Solid


FOUR_C_NAMESPACE_CLOSE

#endif
