/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional torsion spring element

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_torsion3.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::Torsion3Type DRT::ELEMENTS::Torsion3Type::instance_;

DRT::ELEMENTS::Torsion3Type& DRT::ELEMENTS::Torsion3Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::Torsion3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Torsion3* object = new DRT::ELEMENTS::Torsion3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Torsion3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TORSION3")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Torsion3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Torsion3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Torsion3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Torsion3Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Torsion3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::Torsion3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["TORSION3"];

  defs["LINE3"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE3", 3)
                      .AddNamedInt("MAT")
                      .AddNamedString("BENDINGPOTENTIAL")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(int id, int owner) : CORE::Elements::Element(id, owner)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(const DRT::ELEMENTS::Torsion3& old) : CORE::Elements::Element(old)
{
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Torsion3 and return pointer to it (public)|
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Torsion3::Clone() const
{
  DRT::ELEMENTS::Torsion3* newelement = new DRT::ELEMENTS::Torsion3(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Print(std::ostream& os) const { return; }


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Torsion3::Shape() const { return CORE::FE::CellType::line3; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  Element::Pack(data);
  AddtoPack(data, bendingpotential_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  bendingpotential_ = static_cast<BendingPotential>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Torsion3::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<CORE::Elements::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::Elements::ParamsInterface> DRT::ELEMENTS::Torsion3::ParamsInterfacePtr()
{
  return interface_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
