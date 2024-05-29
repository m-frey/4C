/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_bele3.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::Bele3Type DRT::ELEMENTS::Bele3Type::instance_;


DRT::ELEMENTS::Bele3Type& DRT::ELEMENTS::Bele3Type::Instance() { return instance_; }


CORE::COMM::ParObject* DRT::ELEMENTS::Bele3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Bele3* object = new DRT::ELEMENTS::Bele3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Bele3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  // Search for "BELE3". If found, search for "_"
  // the number after "_" is numdof: so BELE3_4 is a BELE3 element
  // with numdof=4
  std::size_t pos = eletype.rfind("BELE3");
  if (pos != std::string::npos)
  {
    if (eletype.substr(pos + 5, 1) == "_")
    {
      std::istringstream is(eletype.substr(pos + 6, 1));

      int numdof = -1;
      is >> numdof;
      Teuchos::RCP<DRT::ELEMENTS::Bele3> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele3(id, owner));
      ele->set_num_dof_per_node(numdof);
      return ele;
    }
    else
    {
      FOUR_C_THROW("ERROR: Found BELE3 element without specified number of dofs!");
    }
  }

  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Bele3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Bele3Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Bele3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::Bele3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs3 = definitions["BELE3_3"];

  defs3["TRI3"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI3", 3).AddOptionalNamedInt("MAT").Build();

  defs3["TRI6"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI6", 6).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD4"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD8"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddOptionalNamedInt("MAT").Build();

  defs3["QUAD9"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddOptionalNamedInt("MAT").Build();

  std::map<std::string, INPUT::LineDefinition>& defs4 = definitions["BELE3_4"];

  defs4["TRI3"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI3", 3).AddOptionalNamedInt("MAT").Build();

  defs4["TRI6"] =
      INPUT::LineDefinition::Builder().AddIntVector("TRI6", 6).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD4"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD4", 4).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD8"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD8", 8).AddOptionalNamedInt("MAT").Build();

  defs4["QUAD9"] =
      INPUT::LineDefinition::Builder().AddIntVector("QUAD9", 9).AddOptionalNamedInt("MAT").Build();
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Bele3LineType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Bele3Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(int id, int owner)
    : CORE::Elements::Element(id, owner), numdofpernode_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(const DRT::ELEMENTS::Bele3& old)
    : CORE::Elements::Element(old), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Bele3::Clone() const
{
  DRT::ELEMENTS::Bele3* newelement = new DRT::ELEMENTS::Bele3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Bele3::Shape() const
{
  switch (num_node())
  {
    case 3:
      return CORE::FE::CellType::tri3;
    case 4:
      return CORE::FE::CellType::quad4;
    case 6:
      return CORE::FE::CellType::tri6;
    case 8:
      return CORE::FE::CellType::quad8;
    case 9:
      return CORE::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // numdofpernode_
  AddtoPack(data, numdofpernode_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // numdofpernode_
  numdofpernode_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Print(std::ostream& os) const
{
  os << "Bele3_" << numdofpernode_ << " " << CORE::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Bele3::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Bele3Line, Bele3>(CORE::COMM::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Bele3::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}


CORE::FE::GaussRule2D DRT::ELEMENTS::Bele3::get_optimal_gaussrule() const
{
  CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::undefined;
  switch (Shape())
  {
    case CORE::FE::CellType::quad4:
      rule = CORE::FE::GaussRule2D::quad_4point;
      break;
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
      rule = CORE::FE::GaussRule2D::quad_9point;
      break;
    case CORE::FE::CellType::tri3:
      rule = CORE::FE::GaussRule2D::tri_3point;
      break;
    case CORE::FE::CellType::tri6:
      rule = CORE::FE::GaussRule2D::tri_6point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
      break;
  }
  return rule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Bele3::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // check if material is defined
  if (linedef->HaveNamed("MAT"))
  {
    int material = 0;
    // read number of material model
    linedef->ExtractInt("MAT", material);
    SetMaterial(0, MAT::Factory(material));
  }
  return true;
}

FOUR_C_NAMESPACE_CLOSE
