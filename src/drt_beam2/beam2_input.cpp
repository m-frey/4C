/*!----------------------------------------------------------------------
\file beam2_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

#include "beam2.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  double shear_correction = 0.0;
  linedef->ExtractDouble("SHEARCORR",shear_correction);
  crosssecshear_ = crosssec_ * shear_correction;

  linedef->ExtractDouble("INERMOM",mominer_);
  return true;
}


#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2::ReadElement()
{
  // read element's nodes; in case of a beam element always line2 shape
  int ierr=0;

  //note: BACI intern type is LINE2, but gid input files work with LIN2
  frchk("LIN2",&ierr);
  // two figures have to be read by frint
  int nnode=2;
  // provide an array of length two in order to store the two figures read
  int nodes[2];
  frint_n("LIN2",nodes,nnode,&ierr);

  // if that does not work try LINE2, in case .dat file was created with pre_exodus
  if (ierr != 1)
  {
    ierr=0;
    frchk("LINE2",&ierr);
    frint_n("LINE2",nodes,nnode,&ierr);
  }

  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  SetMaterial(material);

  // read beam cross section
  crosssec_ = 1.0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");

   // read beam cross section
  double shear_correction = 0.0;
  frdouble("SHEARCORR",&shear_correction,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  crosssecshear_ = crosssec_ * shear_correction;

  // read beam moment of inertia of area
  mominer_ = 1.0;
  frdouble("INERMOM",&mominer_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");


  return true;
} // Beam2::ReadElement()
#endif



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
