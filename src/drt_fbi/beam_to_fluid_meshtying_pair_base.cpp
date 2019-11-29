/*!
\file beam_to_fluid_meshtying_pair_base.cpp

\brief Base meshtying element for meshtying between a 3D beam and a 3D fluid element.

\level 3

\maintainer Nora Hagmeyer
*/


#include "beam_to_fluid_meshtying_pair_base.H"

#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_beaminteraction/beam3contact_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_beaminteraction/beam_to_solid_vtu_output_writer_base.H"
#include "../drt_beaminteraction/beam_to_solid_vtu_output_writer_visualization.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"

#include "../drt_beaminteraction/beam_contact_params.H"
#include "../drt_beaminteraction/beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"
#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"


/**
 *
 */
template <typename beam, typename fluid>
BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::BeamToFluidMeshtyingPairBase()
    : BeamContactPair(), meshtying_is_evaluated_(false), line_to_volume_segments_()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
    std::vector<DRT::Element const*> elements)
{
  // Call Init of base class, the geometry pair will be created and initialized there.
  BeamContactPair::Init(params_ptr, geometry_evaluation_data_ptr, elements);
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam,
    fluid>::Setup()  // todo Muss ich die Positionen nochmal berehcnen obwohl ich sie in der Suche
                     // schon berechnet hab? Vielleicht in extra Funktion reingeben und immer vor
                     // PreEvaluate rufen?
{
  CheckInit();

  // Call setup of base class first.
  BeamContactPair::Setup();

  // Initialize reference nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1posref_(i) = 0.0;

  // Set reference nodal positions (and tangents) for beam element
  for (unsigned int n = 0; n < beam::n_nodes_; ++n)
  {
    const DRT::Node* node = Element1()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele1posref_(3 * beam::n_val_ * n + d) = node->X()[d];

    // tangents
    if (beam::n_val_ == 2)
    {
      LINALG::Matrix<3, 1> tan;
      const DRT::ElementType& eot = Element1()->ElementType();
      if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        dserror("ERROR: Beam3tofluidmeshtying: beam::n_val_=2 detected for beam3 element");
      }
      else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(Element1());
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror(
              "ERROR: Beam3tosolidmeshtying: beam::n_val_=2 detected for beam3r element w/o "
              "Hermite CL");
      }
      else if (eot == DRT::ELEMENTS::Beam3kType::Instance())
      {
        const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(Element1());
        tan = ele->Tref()[n];
      }
      else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
      {
        const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(Element1());
        tan = ele->Tref()[n];
      }
      else
      {
        dserror("ERROR: Beam3tosolidmeshtying: Invalid beam element type");
      }

      for (int d = 0; d < 3; ++d) ele1posref_(3 * beam::n_val_ * n + d + 3) = tan(d, 0);
    }
  }

  // Initialize reference nodal positions for fluid element
  for (unsigned int i = 0; i < fluid::n_dof_; i++) ele2posref_(i) = 0.0;

  // Set reference nodal positions for fluid element
  for (unsigned int n = 0; n < fluid::n_nodes_; ++n)
  {
    const DRT::Node* node = Element2()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele2posref_(3 * n + d) = node->X()[d];
  }

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1pos_(i) = 0.0;

  // Initialize current nodal velocities for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1vel_(i) = 0.0;

  // Initialize current nodal positions for fluid element
  for (unsigned int i = 0; i < fluid::n_dof_; i++) ele2pos_(i) = 0.0;

  // Initialize current nodal velocities for fluid element
  for (unsigned int i = 0; i < fluid::n_dof_; i++) ele2vel_(i) = 0.0;

  issetup_ = true;
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::CreateGeometryPair(
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
{
  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, beam, fluid>(
      geometry_evaluation_data_ptr);
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    CastGeometryPair()->PreEvaluate(ele1poscur_, ele2poscur_, line_to_volume_segments_);
  }
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam,
    fluid>::ResetState(  // todo somehow hand in nodal velocities
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& fluid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < beam::n_dof_; i++)
  {
    ele1pos_(i) = TYPE_BTS_VMT_AD(beam::n_dof_ + fluid::n_dof_, i, beam_centerline_dofvec[i]);
    ele1poscur_(i) = beam_centerline_dofvec[i];
    ele1vel_(i) =
        TYPE_BTS_VMT_AD(beam::n_dof_ + fluid::n_dof_, i, beam_centerline_dofvec[beam::n_dof_ + i]);
  }

  // Fluid element.
  for (unsigned int i = 0; i < fluid::n_dof_; i++)
  {
    ele2pos_(i) =
        TYPE_BTS_VMT_AD(beam::n_dof_ + fluid::n_dof_, beam::n_dof_ + i, fluid_nodal_dofvec[i]);
    ele2poscur_(i) = fluid_nodal_dofvec[i];
    ele2vel_(i) = TYPE_BTS_VMT_AD(
        beam::n_dof_ + fluid::n_dof_, beam::n_dof_ + i, fluid_nodal_dofvec[fluid::n_dof_ + i]);
    /*    if (fluid_nodal_dofvec[fluid::n_dof_ + i] > 0.5)
          printf("Velocity greater than 0.5 while reading in! It's:%f\n",
              fluid_nodal_dofvec[fluid::n_dof_ + i]);
              */
  }
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::Print(std::ostream& out) const
{
  CheckInitSetup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToFluidMeshtyingPair"
      << "\nBeam EleGID:  " << Element1()->Id() << "\nFluid EleGID: " << Element2()->Id();

  out << "\n\nele1 dofvec: " << ele1pos_;
  out << "\nele2 dofvec: " << ele2pos_;
  out << "\nn_segments: " << line_to_volume_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}


/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam,
    fluid>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Only display information if a segment exists for this pair.
  if (line_to_volume_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << Element1()->Id() << ", fluid ID " << Element2()->Id() << ":";
  out << " n_segments = " << line_to_volume_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < line_to_volume_segments_.size();
       index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << line_to_volume_segments_[index_segment].GetEtaA() << ", "
        << line_to_volume_segments_[index_segment].GetEtaB() << "]";
    out << ", Gauss points = "
        << line_to_volume_segments_[index_segment].GetNumberOfProjectionPoints();
    out << "\n";
  }
}

/**
 *
 */
template <typename beam, typename fluid>
void BEAMINTERACTION::BeamToFluidMeshtyingPairBase<beam, fluid>::GetPairVisualization(
    Teuchos::RCP<BeamToSolidVtuOutputWriterBase>
        visualization_writer,  // todo overload outputwriter nicht vergessen!
    const Teuchos::ParameterList& visualization_params) const
{
}

/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
