/*!----------------------------------------------------------------------
\file so_nurbs27_input.cpp

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "so_nurbs27.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::NURBS::So_nurbs27::ReadElement(const std::string& eletype,
                                                   const std::string& distype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  const int numgp=27;
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(numgp, linedef);

  // read possible gaussian points, obsolete for computation
  std::vector<int> ngp;
  linedef->ExtractIntVector("GP",ngp);
  for (int i=0; i<3; ++i)
    if (ngp[i]!=3)
      dserror("Only version with 3 GP for So_N27 implemented");

  // we expect kintype to be total lagrangian
  kintype_ = INPAR::STR::kinem_nonlinearTotLag;

  // check if material kinematics is compatible to element kinematics
  so3mat->ValidKinematics(kintype_);

  return true;
}
