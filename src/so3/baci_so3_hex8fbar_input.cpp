/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Hex8 element with F-bar modification

\level 1

*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_hex8fbar.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8fbar::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  if (buffer == "linear")
  {
    dserror("Only nonlinear kinematics for SO_HEX8FBAR implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
    dserror("Reading SO_HEX8FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}

BACI_NAMESPACE_CLOSE
