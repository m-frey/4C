/*----------------------------------------------------------------------*/
/*! \file

\brief a common base class for all solid elements

\level 2


 *----------------------------------------------------------------------*/


#include "4C_so3_base.hpp"

#include "4C_mat_so3_material.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoBase::SoBase(int id, int owner)
    : CORE::Elements::Element(id, owner),
      kintype_(INPAR::STR::KinemType::vague),
      interface_ptr_(Teuchos::null),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoBase::SoBase(const DRT::ELEMENTS::SoBase& old)
    : CORE::Elements::Element(old),
      kintype_(old.kintype_),
      interface_ptr_(old.interface_ptr_),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoBase::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data, kintype_);

  // material post setup routine
  AddtoPack(data, static_cast<int>(material_post_setup_));
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoBase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  // material post setup routine
  material_post_setup_ = (ExtractInt(position, data) != 0);
}

/*----------------------------------------------------------------------*
 |  return solid material                                      (public) |
 |                                                           seitz 03/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::So3Material> DRT::ELEMENTS::SoBase::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<MAT::So3Material>(
      CORE::Elements::Element::Material(nummat), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoBase::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<Teuchos::RCP<CORE::Elements::ParamsInterface>>("interface");
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::Elements::ParamsInterface> DRT::ELEMENTS::SoBase::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::ELEMENTS::ParamsInterface& DRT::ELEMENTS::SoBase::str_params_interface()
{
  if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
  return *(Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(interface_ptr_, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoBase::error_handling(const double& det_curr, Teuchos::ParameterList& params,
    const int line_id, const STR::ELEMENTS::EvalErrorFlag flag)
{
  // check, if errors are tolerated or should throw a FOUR_C_THROW
  if (IsParamsInterface() and str_params_interface().IsTolerateErrors())
  {
    str_params_interface().SetEleEvalErrorFlag(flag);
    return;
  }
  else
  {
    // === DEPRECATED (hiermeier, 11/17) ==================================
    bool error_tol = false;
    if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
    if (error_tol)
    {
      params.set<bool>("eval_error", true);
      return;
    }
    // if the errors are not tolerated, throw a FOUR_C_THROW
    else
    {
      if (det_curr == 0.0)
        FOUR_C_THROW("ZERO DETERMINANT DETECTED in line %d", line_id);
      else if (det_curr < 0.0)
        FOUR_C_THROW("NEGATIVE DETERMINANT DETECTED in line %d (value = %.5e)", line_id, det_curr);
    }
  }
}

// Check, whether the material post setup routine was
void DRT::ELEMENTS::SoBase::ensure_material_post_setup(Teuchos::ParameterList& params)
{
  if (!material_post_setup_)
  {
    material_post_setup(params);
  }
}

void DRT::ELEMENTS::SoBase::material_post_setup(Teuchos::ParameterList& params)
{
  // This is the minimal implementation. Advanced materials may need extra implementation here.
  SolidMaterial()->post_setup(params, Id());
  material_post_setup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
