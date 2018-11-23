/*!
\file beam_to_solid_volume_meshtying_params.cpp

\brief Data container holding all beam to solid volume meshtying input parameters.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_inpar/inpar_beaminteraction.H"

#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::BeamToSolidVolumeMeshtyingParams()
    : isinit_(false),
      issetup_(false),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID VOLUME MESHTYING");

  // Get parameters form input file.
  {
    // Penalty parameter.
    penalty_parameter_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER");
    if (penalty_parameter_ < 0.0)
      dserror("beam-to-volume-meshtying penalty parameter must not be negative!");

    // Gauss rule for integration along the beam.
    gauss_rule_ = INPAR::BEAMINTERACTION::IntToGaussRule1D(
        beam_to_solid_contact_params_list.get<int>("GAUSS_POINTS"));
  }

  isinit_ = true;
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Setup()
{
  CheckInit();

  // Empty for now.

  issetup_ = true;
}
