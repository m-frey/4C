/*---------------------------------------------------------------------*/
/*!
\file drt_meshfree_multibin.cpp

\brief element type class of meshfree multibin, creating the same

\maintainer Georg Hammerl

\level 2

*/
/*---------------------------------------------------------------------*/


#include "drt_meshfree_multibin.H"

/// class MeshfreeMultiBinType
DRT::MESHFREE::MeshfreeMultiBinType DRT::MESHFREE::MeshfreeMultiBinType::instance_;

DRT::MESHFREE::MeshfreeMultiBinType& DRT::MESHFREE::MeshfreeMultiBinType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::MESHFREE::MeshfreeMultiBinType::Create( const std::vector<char> & data )
{
  DRT::MESHFREE::MeshfreeMultiBin* object =
    new DRT::MESHFREE::MeshfreeMultiBin(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create( const std::string eletype,
                                                                         const std::string eledistype,
                                                                         const int id,
                                                                         const int owner )
{
  if (eletype=="MESHFREEMULTIBIN")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id,owner));
  return ele;
}



/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(int id, int owner)
:   DRT::MESHFREE::MeshfreeBin<DRT::Element>(id,owner),
    associatedeleid_(INPAR::BINSTRATEGY::enumsize),
    associatedele_(INPAR::BINSTRATEGY::enumsize)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(const DRT::MESHFREE::MeshfreeMultiBin& old)
:   DRT::MESHFREE::MeshfreeBin<DRT::Element>(old),
    associatedeleid_(old.associatedeleid_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::~MeshfreeMultiBin()
{
  return;
}

/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::MESHFREE::MeshfreeMultiBin::Clone() const
{
  DRT::MESHFREE::MeshfreeMultiBin* newele = new DRT::MESHFREE::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 |  << operator                                                 ghamm 04/13 |
 *--------------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeMultiBin& bin)
{
  bin.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                      (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Print(std::ostream& os) const
{
  os << "MeshfreeMultiBin ";
  DRT::Element::Print(os);

  const int nwallele = NumAssociatedEle(INPAR::BINSTRATEGY::Surface);
  const int* walleleids = AssociatedEleIds(INPAR::BINSTRATEGY::Surface);
  if (nwallele > 0)
  {
    os << " Associated wall elements ";
    for (int j=0; j<nwallele; ++j) os << std::setw(10) << walleleids[j] << " ";
  }

  const int nfluidele = NumAssociatedEle(INPAR::BINSTRATEGY::Volume);
  const int* wfluideleids = AssociatedEleIds(INPAR::BINSTRATEGY::Volume);
  if (nfluidele > 0)
  {
    os << " Associated fluid elements ";
    for (int j=0; j<nfluidele; ++j) os << std::setw(10) << wfluideleids[j] << " ";
  }

  const int nbeamele = NumAssociatedEle(INPAR::BINSTRATEGY::Beam);
  const int* wbeameleids = AssociatedEleIds(INPAR::BINSTRATEGY::Beam);
  if (nbeamele > 0)
  {
    os << " Associated beam elements ";
    for (int j=0; j<nbeamele; ++j) os << std::setw(10) << wbeameleids[j] << " ";
  }

  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single element from the bin                (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::DeleteAssociatedEle(INPAR::BINSTRATEGY::BinContent bin_content, int gid)
{
  for (unsigned int i = 0; i<associatedeleid_[bin_content].size(); i++){
    if (associatedeleid_[bin_content][i]==gid){
      associatedeleid_[bin_content].erase(associatedeleid_[bin_content].begin()+i);
      associatedele_[bin_content].erase(associatedele_[bin_content].begin()+i);
      return;
    }
  }
  dserror("Connectivity issues: No element with specified gid to delete in bin. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Delete all wall elements from current bin           (public) ghamm 09/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::RemoveAssociatedEles(INPAR::BINSTRATEGY::BinContent bin_content)
{
  associatedeleid_[bin_content].clear();
  associatedele_[bin_content].clear();

  return;
}

/*--------------------------------------------------------------------------*
 |  Build element pointers                             (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeMultiBin::BuildElePointers(INPAR::BINSTRATEGY::BinContent bin_content, DRT::Element** eles)
{
  associatedele_[bin_content].resize(NumAssociatedEle(bin_content));
  for (int i=0; i<NumAssociatedEle(bin_content); ++i) associatedele_[bin_content][i] = eles[i];
  return true;
}

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Element
  DRT::Element::Pack(data);
  // add vector associatedeleid_
  for(int i=0;i<INPAR::BINSTRATEGY::enumsize;++i)
    AddtoPack(data,associatedeleid_[i]);

  return;
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class DRT::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Element::Unpack(basedata);
  // extract associatedeleid_
  for(int i=0;i<INPAR::BINSTRATEGY::enumsize;++i)
    ExtractfromPack(position,data,associatedeleid_[i]);
  // associatedele_ is NOT communicated
  associatedele_.clear();
  return;
}
