/*----------------------------------------------------------------------*/
/*! \file

\brief Base element for interactions between a beam and a solid.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_pair_base.H"

#include "baci_beam3_euler_bernoulli.H"
#include "baci_beam3_kirchhoff.H"
#include "baci_beam3_reissner.H"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_base.H"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.H"
#include "baci_geometry_pair_element_functions.H"
#include "baci_geometry_pair_scalar_types.H"

#include <Sacado.hpp>


/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam,
    solid>::BeamToSolidPairBase()
    : BeamContactPair(), line_to_3D_segments_()
{
}


/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam, solid>::Setup()
{
  CheckInit();

  // Call setup of base class first.
  BeamContactPair::Setup();

  // Set reference nodal positions (and tangents) for beam element
  for (unsigned int n = 0; n < beam::n_nodes_; ++n)
  {
    const DRT::Node* node = Element1()->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele1posref_(3 * beam::n_val_ * n + d) = node->X()[d];

    // tangents
    if (beam::n_val_ == 2)
    {
      CORE::LINALG::Matrix<3, 1> tan;
      const DRT::ElementType& eot = Element1()->ElementType();

      if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(Element1());
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror(
              "ERROR: Beam3tosolidmeshtying: beam::n_val_=2 detected for beam3r element w/o "
              "Hermite centerline");
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

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < beam::n_dof_; i++) ele1pos_(i) = 0.0;

  issetup_ = true;
}

/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam,
    solid>::ResetState(const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Beam element.
  for (unsigned int i = 0; i < beam::n_dof_; i++)
    ele1pos_(i) = CORE::FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + solid::n_dof_, i, beam_centerline_dofvec[i]);
}

/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam,
    solid>::SetRestartDisplacement(const std::vector<std::vector<double>>& centerline_restart_vec_)
{
  // Call the parent method.
  BeamContactPair::SetRestartDisplacement(centerline_restart_vec_);
}

/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam, solid>::Print(
    std::ostream& out) const
{
  CheckInitSetup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToSolidPairBase"
      << "\nBeam EleGID:  " << Element1()->Id() << "\nSolid EleGID: " << Element2()->Id();

  out << "\n\nbeam dofvec: " << ele1pos_;
  out << "\nn_segments: " << line_to_3D_segments_.size();
  out << "\n";
  out << "------------------------------------------------------------------------\n";
}

/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam,
    solid>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Only display information if a segment exists for this pair.
  if (line_to_3D_segments_.size() == 0) return;

  // Display the number of segments and segment length.
  out << "beam ID " << Element1()->Id() << ", solid ID " << Element2()->Id() << ":";
  out << " n_segments = " << line_to_3D_segments_.size() << "\n";

  // Loop over segments and display information about them.
  for (unsigned int index_segment = 0; index_segment < line_to_3D_segments_.size(); index_segment++)
  {
    out << "    segment " << index_segment << ": ";
    out << "eta in [" << CORE::FADUTILS::CastToDouble(line_to_3D_segments_[index_segment].GetEtaA())
        << ", " << CORE::FADUTILS::CastToDouble(line_to_3D_segments_[index_segment].GetEtaB())
        << "]";
    out << ", Gauss points = " << line_to_3D_segments_[index_segment].GetNumberOfProjectionPoints();
    out << "\n";
  }
}

/**
 *
 */
template <typename scalar_type, typename segments_scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidPairBase<scalar_type, segments_scalar_type, beam,
    solid>::EvaluateBeamPositionDouble(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>&
                                           integration_point,
    CORE::LINALG::Matrix<3, 1, double>& r_beam, bool reference) const
{
  if (reference)
    GEOMETRYPAIR::EvaluatePosition<beam>(
        integration_point.GetEta(), ele1posref_, r_beam, this->Element1());
  else
    GEOMETRYPAIR::EvaluatePosition<beam>(integration_point.GetEta(),
        CORE::FADUTILS::CastToDouble(ele1pos_), r_beam, this->Element1());
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  // Beam-to-volume pairs
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_hex8>, double,
      t_hermite, t_hex8>;
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_hex20>, double,
      t_hermite, t_hex20>;
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_hex27>, double,
      t_hermite, t_hex27>;
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_tet4>, double,
      t_hermite, t_tet4>;
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_tet10>, double,
      t_hermite, t_tet10>;
  template class BeamToSolidPairBase<line_to_volume_scalar_type<t_hermite, t_nurbs27>, double,
      t_hermite, t_nurbs27>;

  // Beam-to-surface pairs
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_quad4>, double, t_line2,
      t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_quad8>, double, t_line2,
      t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_quad9>, double, t_line2,
      t_quad9>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_tri3>, double, t_line2,
      t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_tri6>, double, t_line2,
      t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_line2, t_nurbs9>, double,
      t_line2, t_nurbs9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_line2, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_line2, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_line2, t_quad9>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_line2, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_line2, t_tri6>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, double, t_line2, t_nurbs9>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_fixed_size<t_line2, t_hex8>,
      double, t_line2, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_fixed_size<t_line2, t_hex20>,
      double, t_line2, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_fixed_size<t_line2, t_hex27>,
      double, t_line2, t_quad9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad9>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>,
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_line2, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_line2, t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_line2, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_line2, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_line2, t_quad9>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>,
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_quad4>, double,
      t_hermite, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_quad8>, double,
      t_hermite, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_quad9>, double,
      t_hermite, t_quad9>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_tri3>, double,
      t_hermite, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_tri6>, double,
      t_hermite, t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_scalar_type<t_hermite, t_nurbs9>, double,
      t_hermite, t_nurbs9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_hermite, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_hermite, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_hermite, t_quad9>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_hermite, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type, double, t_hermite, t_tri6>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, double, t_hermite,
      t_nurbs9>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, double, t_hermite, t_quad4>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, double, t_hermite, t_quad8>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, double, t_hermite, t_quad9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type_1st_order,
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>,
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_hermite, t_tri3>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_hermite, t_tri6>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_hermite, t_quad4>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_hermite, t_quad8>;
  template class BeamToSolidPairBase<line_to_surface_patch_scalar_type,
      line_to_surface_patch_scalar_type, t_hermite, t_quad9>;
  template class BeamToSolidPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>,
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace BEAMINTERACTION
