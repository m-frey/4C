/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam3_base.cpp

\brief base class for all beam elements

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3_base.H"

#include "../drt_mat/beam_elasthyper.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Sacado.hpp>
#include "../drt_beaminteraction/periodic_boundingbox.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(int id, int owner) :
DRT::Element(id,owner),
interface_ptr_(Teuchos::null),
browndyn_interface_ptr_(Teuchos::null)
{
  // todo: this is a temporary hack, should of course be set from outside
  for (unsigned int i=0; i<1; ++i)
  {
    bspotposxi_.push_back(0.0);
    bspotstatus_[i] = -1;
  }
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(const DRT::ELEMENTS::Beam3Base& old) :
 DRT::Element(old),
 bspotposxi_(old.bspotposxi_),
 bspotstatus_(old.bspotstatus_)
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  // bspotposxi_
  AddtoPack(data,bspotposxi_);
  // bspotstatus_
  AddtoPack(data,bspotstatus_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // bspotposxi_
  ExtractfromPack(position,data,bspotposxi_);
  // bspotstatus_
  ExtractfromPack(position,data,bspotstatus_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ =
        Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>
        (p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> >("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::SetBrownianDynParamsInterfacePtr()
{
  browndyn_interface_ptr_ = interface_ptr_->GetBrownianDynParamInterface();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Beam3Base::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BROWNIANDYN::ParamsInterface> DRT::ELEMENTS::Beam3Base::BrownianDynParamsInterfacePtr() const
{
  return browndyn_interface_ptr_;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetAdditiveDofGIDs(
    const DRT::Discretization& discret,
    const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->PositionDofIndices(dofindices,node);
  this->TangentDofIndices(dofindices,node);
  this->Rotation1DDofIndices(dofindices,node);
  this->TangentLengthDofIndices(dofindices,node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i=0; i<dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0,&node,dofindices[i]));

  return dofgids;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetRotVecDofGIDs(
    const DRT::Discretization& discret,
    const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->RotationVecDofIndices(dofindices,node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i=0; i<dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0,&node,dofindices[i]));

  return dofgids;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3Base::GetCircularCrossSectionRadiusForInteractions() const
{
  return GetBeamMaterial().GetInteractionRadius();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetRefPosAtXi(LINALG::Matrix<3,1>& refpos,
                                             const double& xi) const
{
  const int numclnodes = this->NumCenterlineNodes();
  const int numnodalvalues = this->HermiteCenterlineInterpolation() ? 2 : 1;

  std::vector<double> zerovec;
  zerovec.resize(3*numnodalvalues*numclnodes);

  this->GetPosAtXi(refpos,xi,zerovec);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
const MAT::BeamElastHyperMaterial& DRT::ELEMENTS::Beam3Base::GetBeamMaterial() const
{
  // Todo @grill think about storing the casted pointer as class variable or other solution to
  //      avoid cast in every element evaluation

  const MAT::BeamElastHyperMaterial* beam_material_ptr = NULL;

  // get the material law
  Teuchos::RCP<const MAT::Material> material_ptr = Material();

  switch ( material_ptr->MaterialType() )
  {
    case INPAR::MAT::m_beam_elast_hyper_generic:
    {
      beam_material_ptr =
         static_cast<const MAT::BeamElastHyperMaterial*>( material_ptr.get() );

      if (beam_material_ptr == NULL)
        dserror("cast to beam material class failed!");

      break;
    }
    default:
    {
      dserror("unknown or improper type of material law! expected beam material law!");
      break;
    }
  }

  return *beam_material_ptr;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices(
    LINALG::TMatrix<T,3,3>& CN,
    LINALG::TMatrix<T,3,3>& CM) const
{
  GetBeamMaterial().GetConstitutiveMatrixOfForcesMaterialFrame( CN );
  GetBeamMaterial().GetConstitutiveMatrixOfMomentsMaterialFrame( CM );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor(
    double& mass_inertia_translational,
    LINALG::TMatrix<T,3,3>& J) const
{
  GetTranslationalMassInertiaFactor( mass_inertia_translational );
  GetBeamMaterial().GetMassMomentOfInertiaTensorMaterialFrame( J );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetTranslationalMassInertiaFactor(
    double& mass_inertia_translational) const
{
  mass_inertia_translational = GetBeamMaterial().GetTranslationalMassInertiaFactor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetDampingCoefficients(LINALG::Matrix<3,1>& gamma) const
{
  /* These are coefficients for a straight cylindrical rod taken from
   * Howard, p. 107, table 6.2. The order is as follows:
   * (0) damping of translation parallel to axis,
   * (1) damping of translation orthogonal to axis,
   * (2) damping of rotation around its own axis */

  gamma(0) = 2.0 * PI * BrownianDynParamsInterface().GetViscosity();
  gamma(1) = 4.0 * PI * BrownianDynParamsInterface().GetViscosity();
  gamma(2) = 4.0 * PI * BrownianDynParamsInterface().GetViscosity()
             * std::pow( GetCircularCrossSectionRadiusForInteractions(), 2.0);

  // huge improvement in convergence of non-linear solver in case of artificial factor 4000
//  gamma(2) *= 4000.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<unsigned int ndim, typename T>
void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity(
    Teuchos::ParameterList&      params,  //!<parameter list
    const LINALG::TMatrix<T,ndim,1>&  evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
    LINALG::TMatrix<T,ndim,1>&        velbackground,  //!< velocity of background fluid
    LINALG::TMatrix<T,ndim,ndim>&     velbackgroundgrad) const//!<gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in direction params.get<int>("DBCDISPDIR",0) and orthogonal to z-axis.
   * In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  //velocity at upper boundary of domain
//  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0.0);
  velbackgroundgrad.PutScalar(0.0);

  // fixme @grill: this needs to go somewhere else, outside element level
//  double time = -1.0;
//
//  double shearamplitude = BrownianDynParamsInterface().GetShearAmplitude();
//  int curvenumber = BrownianDynParamsInterface().GetCurveNumber() - 1;
//  int dbcdispdir = BrownianDynParamsInterface().GetDbcDispDir() - 1;
//
//  Teuchos::RCP<std::vector<double> > periodlength = BrownianDynParamsInterface().GetPeriodLength();
//  INPAR::STATMECH::DBCType dbctype = BrownianDynParamsInterface().GetDbcType();
//  bool shearflow = false;
//  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
//     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
//     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
//     dbctype==INPAR::STATMECH::dbctype_affineshear||
//     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
//    shearflow = true;
//
//  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
//  if(periodlength->at(0) > 0.0)
//    if(shearflow &&  curvenumber >=  0 && dbcdispdir >= 0 )
//    {
//      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1))[1];
//
//      //compute background velocity
//      velbackground(dbcdispdir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;
//
//      //compute gradient of background velocity
//      velbackgroundgrad(dbcdispdir,ndim-1) = uppervel / periodlength->at(ndim-1);
//    }
}

/*-----------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of          |
 | periodic boundary conditions; if two nodes within one element are se-       |
 | parated by a periodic boundary, one of them is shifted such that the final  |
 | distance in R^3 is the same as the initial distance in the periodic         |
 | space; the shift affects computation on element level within that           |
 | iteration step, only (no change in global variables performed)              |
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(
    std::vector<double>& disp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& periodic_boundingbox) const
{
  /* get number of degrees of freedom per node; note:
   * the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  // get number of nodes that are used for centerline interpolation
  unsigned int nnodecl = NumCenterlineNodes();

  // loop through all nodes except for the first node which remains
  // fixed as reference node
  for(unsigned int i = 1; i < nnodecl; ++i)
    for( int dim = 0; dim < 3; ++dim )
      periodic_boundingbox->UnShift1D( dim, disp[numdof*i+dim],
          Nodes()[0]->X()[dim] + disp[numdof*0+dim], Nodes()[i]->X()[dim]);
}

//! shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(std::vector<double>& disp) const
{
  this->UnShiftNodePosition(disp, BrownianDynParamsInterface().GetPeriodicBoundingBox() );
}
/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetPosOfBindingSpot(LINALG::Matrix<3,1>& pos,
                                                   std::vector<double>& disp,
                                                   const int& bspotlocn,
                                                   Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& periodic_boundingbox) const
{
  const double xi = bspotposxi_[bspotlocn];
  // get position
  GetPosAtXi(pos,xi,disp);

  // check if pos at xi lies outside the periodic box, if it does, shift it back in
  periodic_boundingbox->Shift3D(pos);

  return;
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetTriadOfBindingSpot(
    LINALG::Matrix<3,3>&                            triad,
    std::vector<double>&                            disp,
    const int&                                      bspotlocn) const
{
  const double xi = bspotposxi_[bspotlocn];
  // get position
  GetTriadAtXi(triad,xi,disp);

}

// explicit template instantiations
template void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices<double>(
    LINALG::TMatrix<double,3,3>& CN,
    LINALG::TMatrix<double,3,3>& CM) const;
template void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices<Sacado::Fad::DFad<double> >(
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3,3>& CN,
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3,3>& CM) const;

template void DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor<double>(
    double&,
    LINALG::TMatrix<double,3,3>&) const;
template void DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor<Sacado::Fad::DFad<double> >(
    double&,
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3,3>&) const;

template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3,double>(
    Teuchos::ParameterList&,
    const LINALG::TMatrix<double,3,1>&,
    LINALG::TMatrix<double,3,1>&,
    LINALG::TMatrix<double,3,3>&) const;
template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3,Sacado::Fad::DFad<double> >(
    Teuchos::ParameterList&,
    const LINALG::TMatrix<Sacado::Fad::DFad<double>,3,1>&,
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3,1>&,
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3,3>&) const;
