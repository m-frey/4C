/*!
\file scatra_element.cpp
\brief A finite element for simulation transport phenomena

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/matlist_reactions.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/myocard.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "scatra_ele_calc_utils.H"
#include "scatra_ele.H"


DRT::ELEMENTS::TransportType DRT::ELEMENTS::TransportType::instance_;


DRT::ParObject* DRT::ELEMENTS::TransportType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Transport* object =
    new DRT::ELEMENTS::Transport(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportType::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="TRANSP" or eletype=="CONDIF2" or eletype=="CONDIF3" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Transport(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Transport(id,owner));
  return ele;
}


void DRT::ELEMENTS::TransportType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;

  if (DRT::Problem::Instance(0)->ProblemType() == prb_elch)
  {
    if (nv > 1) // only when we have more than 1 dof per node!
    {
      nv -= 1; // ion concentrations
      np = 1;  // electric potential
    }
  }
}

void DRT::ELEMENTS::TransportType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluidDNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::TransportType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["TRANSP"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
   ;

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS2"]
    .AddIntVector("NURBS2",2)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

  defs["NURBS3"]
    .AddIntVector("NURBS3",3)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddOptionalNamedDoubleVector("FIBER1",3)
    ;

}


DRT::ELEMENTS::TransportBoundaryType DRT::ELEMENTS::TransportBoundaryType::instance_;

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportBoundaryType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new TransportBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(int id, int owner) :
DRT::Element(id,owner),
data_(),
numdofpernode_(-1),
distype_(dis_none),
impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(const DRT::ELEMENTS::Transport& old) :
DRT::Element(old),
data_(old.data_),
numdofpernode_(old.numdofpernode_),
distype_(old.distype_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Transport and return pointer to it (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Transport::Clone() const
{
  DRT::ELEMENTS::Transport* newelement = new DRT::ELEMENTS::Transport(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::SetMaterial(int matnum)
{
  // the standard part:
  DRT::Element::SetMaterial(matnum);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  Teuchos::RCP<MAT::Material> mat = Material();
  if(mat->MaterialType() == INPAR::MAT::m_scatra or
     mat->MaterialType() == INPAR::MAT::m_scatra_aniso or
     mat->MaterialType() == INPAR::MAT::m_myocard or
     mat->MaterialType() == INPAR::MAT::m_mixfrac or
     mat->MaterialType() == INPAR::MAT::m_sutherland or
     mat->MaterialType() == INPAR::MAT::m_arrhenius_pv or
     mat->MaterialType() == INPAR::MAT::m_ferech_pv or
     mat->MaterialType() == INPAR::MAT::m_ion or
     mat->MaterialType() == INPAR::MAT::m_biofilm or
     mat->MaterialType() == INPAR::MAT::m_th_fourier_iso or
     mat->MaterialType() == INPAR::MAT::m_thermostvenant or
     mat->MaterialType() == INPAR::MAT::m_yoghurt or
     mat->MaterialType() == INPAR::MAT::m_scatra_growth_scd
     )
    numdofpernode_ = 1; // we only have a single scalar
  else if(mat->MaterialType() == INPAR::MAT::m_electrode)
    numdofpernode_ = 2; // concentration and electric potential
  else if(mat->MaterialType() == INPAR::MAT::m_matlist) // we have a system of scalars
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(mat.get());
    numdofpernode_=actmat->NumMat();

    // for problem type ELCH we have one additional degree of freedom per node
    // for the electric potential
    if (DRT::Problem::Instance()->ProblemType()== prb_elch)
    {
      for (int ii=0; ii<numdofpernode_; ++ii)
      {
        // In the context of ELCH the only valid material combination is m_matlist and m_ion
        if(actmat->MaterialById(actmat->MatID(ii))->MaterialType() != INPAR::MAT::m_ion)
          dserror("In the context of ELCH the material Mat_matlist can be only used in combination with Mat_ion");
      }
      numdofpernode_ += 1;
    }
  }
  else if(mat->MaterialType() == INPAR::MAT::m_matlist_reactions) // we have a system of reactive scalars
    {
      const MAT::MatListReactions* actmat = static_cast<const MAT::MatListReactions*>(mat.get());
      numdofpernode_=actmat->NumMat();

      for (int ii=0; ii<numdofpernode_; ++ii)
      {
        // In the context of reactions the only valid material combination is m_matlist and m_scatra
        if(actmat->MaterialById(actmat->MatID(ii))->MaterialType() != INPAR::MAT::m_scatra)
          dserror("The material Mat_matlist_reaction only supports MAT_scatra as valid main Material");
      }

      int numreac = actmat->NumReac();
      for (int jj=0; jj<numreac; ++jj)
      {
        // In the context of reactions the only valid material combination is m_matlist and m_scatra_reaction
        if(actmat->MaterialById(actmat->ReacID(jj))->MaterialType() != INPAR::MAT::m_scatra_reaction)
          dserror("The material MAT_matlist_reaction only supports MAT_scatra_reaction as valid reaction Material");
      }
  }
  else if(mat->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(mat.get());

    numdofpernode_ = actmat->NumDOF();
  }
  else
    dserror("Transport element got unsupported material type %d", mat->MaterialType());

  return;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::SetMaterial(int matnum,DRT::Element* oldele)
{
  SetMaterial(matnum);

  Teuchos::RCP<MAT::Material> mat = Material();

  if(mat->MaterialType() == INPAR::MAT::m_myocard)
  {
    Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(mat);

    Teuchos::RCP<MAT::ElastHyper> somat = Teuchos::rcp_dynamic_cast<MAT::ElastHyper>(oldele->Material());
    if(somat==Teuchos::null)
      dserror("cast to ElastHyper failed");

    //copy fiber information from solid material to scatra material (for now, only one fiber vector)
    std::vector<LINALG::Matrix<3,1> > fibervecs(0);
    somat->GetFiberVecs(fibervecs);
    actmat->Setup(fibervecs[0]);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Transport element                      (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Transport::Shape() const
{
  return distype_;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  Element::Pack(data);

  // add internal data
  AddtoPack(data,data_);
  AddtoPack(data,numdofpernode_);
  AddtoPack(data,distype_);
  AddtoPack(data,impltype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // extract internal data
  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);
  ExtractfromPack(position,data,numdofpernode_);
  distype_ = static_cast<DiscretizationType>(ExtractInt(position,data));
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position,data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)           gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)        gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)          gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumVolume() const
{
  return DRT::UTILS::getNumberOfElementVolumes(distype_);
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::~Transport()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Print(std::ostream& os) const
{
  os << "Transport element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << DRT::DistypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  std::cout << "Type of scalar transport: " << SCATRA::ImplTypeToString(impltype_) << std::endl;
  std::cout << std::endl;
  std::cout << data_;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                  g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Transport::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumLine() > 1) // 3D and 2D
    return DRT::UTILS::ElementBoundaryFactory<TransportBoundary,Transport>(DRT::UTILS::buildLines,this);
  else
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > lines(1);
    lines[0]= Teuchos::rcp(this, false);
    return lines;
  }

}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Transport::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  if (NumSurface() > 1) // 3D
    return DRT::UTILS::ElementBoundaryFactory<TransportBoundary,Transport>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    // 1D
    dserror("Surfaces() for 1D-Transport element not implemented");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Transport::Volumes()
{
  if (NumVolume() == 1)
  {
    std::vector<Teuchos::RCP<Element> > volumes(1);
    volumes[0]= Teuchos::rcp(this, false);
    return volumes;
  }
  else
  {
    dserror("Volumes() for 1D-/2D-Transport element not implemented");
    return DRT::Element::Volumes();
  }
}


/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                gjb 01/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::VisNames(std::map<std::string,int>& names)
{

  // see whether we have additional data for visualization in our container
  for (int k = 0 ;k<numdofpernode_; k++)
  {
    std::ostringstream temp;
    temp << k;

    // element Peclet number
    std::string name = "Pe_"+temp.str();
    const std::vector<double>* Pe = data_.Get<std::vector<double> >(name);
    if (Pe) names.insert(std::pair<std::string,int>(name,1));

    // element Peclet number (only migration term)
    name = "Pe_mig_"+temp.str();
    const std::vector<double>* Pe_mig = data_.Get<std::vector<double> >(name);
    if (Pe_mig) names.insert(std::pair<std::string,int>(name,1));

    //characteristic element length
    name = "hk_"+temp.str();
    const std::vector<double>* hk = data_.Get<std::vector<double> >(name);
    if (hk) names.insert(std::pair<std::string,int>(name,1));

    // Stabilization parameter at element center
    name = "tau_"+temp.str();
    const std::vector<double>* tau = data_.Get<std::vector<double> >(name);
    if (tau) names.insert(std::pair<std::string,int>(name,1));

  } // loop over transported scalars

  return;
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         gjb 01/09|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport ::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if(DRT::Element::VisData(name,data))
    return true;

  for (int k = 0 ;k<numdofpernode_; k++)
  {
    std::ostringstream temp;
    temp << k;
    if (   (name == "Pe_"+temp.str()    )
        || (name == "Pe_mig_"+temp.str())
        || (name == "hk_"+temp.str()    )
        || (name == "tau_"+temp.str()   )
    )
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      const double value = data_.GetDouble(name);
      data[0] = value;
      return true;
    }
  } // loop over transported scalars

  return false;
}


/*----------------------------------------------------------------------*
 | set implementation type                                   fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport ::SetImplType(const INPAR::SCATRA::ImplType impltype)
{
  // set implementation type
  impltype_ = impltype;

  return;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Transport* parent,
                              const int lbeleid) :
DRT::Element(id,owner)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent,lbeleid);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(const DRT::ELEMENTS::TransportBoundary& old) :
DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::TransportBoundary::Clone() const
{
  DRT::ELEMENTS::TransportBoundary* newelement = new DRT::ELEMENTS::TransportBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::TransportBoundary::Shape() const
{
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Pack(DRT::PackBuffer& data) const
{
  dserror("This TransportBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Unpack(const std::vector<char>& data)
{
  dserror("This TransportBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::~TransportBoundary()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Print(std::ostream& os) const
{
  os << "TransportBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  "<<Shape()<< std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)        gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)    gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::TransportBoundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of TransportBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::TransportBoundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of TransportBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}

