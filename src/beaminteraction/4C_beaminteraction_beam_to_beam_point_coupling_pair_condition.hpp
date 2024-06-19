/*----------------------------------------------------------------------*/
/*! \file

\brief Class to manage beam-to-beam point couplings.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_POINT_COUPLING_PAIR_CONDITION_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_POINT_COUPLING_PAIR_CONDITION_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_conditions.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BEAMINTERACTION
{
  /**
   * \brief This base class represents a single beam point coupling condition.
   */
  class BeamToBeamPointCouplingCondition : public BeamInteractionConditionBase
  {
   public:
    BeamToBeamPointCouplingCondition(
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_line,
        double positional_penalty_parameter, double rotational_penalty_parameter)
        : BeamInteractionConditionBase(condition_line),
          positional_penalty_parameter_(positional_penalty_parameter),
          rotational_penalty_parameter_(rotational_penalty_parameter)
    {
    }
    /**
     * \brief Check if a combination of beam ids is in this condition.
     */
    bool IdsInCondition(const int id_line, const int id_other) const override;

    /**
     * \brief Clear not reusable data (derived).
     */
    void clear() override;

    /**
     * \brief Create the beam contact pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) override;

    /**
     * \brief Build the ID sets for this condition. The ID sets will be used to check if an element
     * is in this condition.
     */
    void BuildIdSets(const Teuchos::RCP<const Core::FE::Discretization>& discretization) override;


   private:
    /// Penalty parameter used to couple the positional DoFs
    double positional_penalty_parameter_;
    /// Penalty parameter used to couple the rotational DoFs
    double rotational_penalty_parameter_;
    /// Element-local parameter coordinates of the coupling nodes
    std::array<double, 2> local_parameter_coordinates_;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
