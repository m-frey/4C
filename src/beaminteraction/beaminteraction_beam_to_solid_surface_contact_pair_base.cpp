/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for contact between a 3D beam and a surface element.

\level 3
*/


#include "beaminteraction_beam_to_solid_surface_contact_pair_base.H"

#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_surface_contact_params.H"
#include "beaminteraction_beam_to_solid_vtu_output_writer_base.H"
#include "beaminteraction_beam_to_solid_vtu_output_writer_visualization.H"
#include "beaminteraction_beam_to_solid_surface_vtk_output_params.H"
#include "beaminteraction_beam_to_solid_utils.H"
#include "beaminteraction_calc_utils.H"
#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_factory.H"
#include "geometry_pair_element_faces.H"
#include "geometry_pair_scalar_types.H"

#include <Epetra_FEVector.h>


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam,
    surface>::BeamToSolidSurfaceContactPairBase()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Clean the segments, as they will be re-evaluated in each iteration.
  this->line_to_3D_segments_.clear();

  // Set the current position of the beam element.
  const int n_patch_dof = face_element_->GetPatchGID().size();
  for (unsigned int i = 0; i < beam::n_dof_; i++)
    this->ele1pos_(i) = FADUTILS::HigherOrderFadValue<scalar_type>::apply(
        beam::n_dof_ + n_patch_dof, i, beam_centerline_dofvec[i]);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  CastGeometryPair()->PreEvaluate(this->ele1pos_, this->face_element_->GetFacePosition(),
      this->line_to_3D_segments_, this->face_element_->GetCurrentNormals());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam,
    surface>::CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>&
        geometry_evaluation_data_ptr)
{
  // Call the method of the base class.
  BeamContactPair::CreateGeometryPair(geometry_evaluation_data_ptr);

  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  this->geometry_pair_ =
      GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<scalar_type, beam, surface>(
          geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::SetFaceElement(
    Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>>(
          face_element, true);

  // The second element in the pair has to be the face element.
  CastGeometryPair()->SetElement2(face_element_->GetDrtFaceElement());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, beam, surface>>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::CastGeometryPair()
    const
{
  return Teuchos::rcp_dynamic_cast<
      GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, beam, surface>>(
      this->geometry_pair_, true);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
std::vector<int>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<scalar_type, beam, surface>::GetPairGID(
    const ::DRT::Discretization& discret) const
{
  // Get the beam centerline GIDs.
  LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  UTILS::GetElementCenterlineGIDIndices(discret, this->Element1(), beam_centerline_gid);

  // Get the patch (in this case just the one face element) GIDs.
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();
  std::vector<int> pair_gid;
  pair_gid.resize(beam::n_dof_ + patch_gid.size());

  // Combine beam and solid GIDs into one vector.
  for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
    pair_gid[i_dof_beam] = beam_centerline_gid(i_dof_beam);
  for (unsigned int i_dof_patch = 0; i_dof_patch < patch_gid.size(); i_dof_patch++)
    pair_gid[beam::n_dof_ + i_dof_patch] = patch_gid[i_dof_patch];

  return pair_gid;
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

}  // namespace BEAMINTERACTION
