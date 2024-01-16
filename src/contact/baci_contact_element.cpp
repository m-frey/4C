/*-----------------------------------------------------------------------*/
/*! \file
\brief a contact element
\level 2
*/
/*-----------------------------------------------------------------------*/

#include "baci_contact_element.H"

#include "baci_contact_friction_node.H"
#include "baci_contact_node.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"

#include <array>

BACI_NAMESPACE_OPEN
CONTACT::ElementType CONTACT::ElementType::instance_;

CONTACT::ElementType& CONTACT::ElementType::Instance() { return instance_; }

CORE::COMM::ParObject* CONTACT::ElementType::Create(const std::vector<char>& data)
{
  CONTACT::Element* ele =
      new CONTACT::Element(0, 0, CORE::FE::CellType::dis_none, 0, nullptr, false);
  ele->Unpack(data);
  return ele;
}

Teuchos::RCP<DRT::Element> CONTACT::ElementType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Element( id, owner ) );
  return Teuchos::null;
}

void CONTACT::ElementType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

CORE::LINALG::SerialDenseMatrix CONTACT::ElementType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  dserror("method ComputeNullSpace not implemented!");
  return nullspace;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Element::Element(int id, int owner, const CORE::FE::CellType& shape, const int numnode,
    const int* nodeids, const bool isslave, bool isnurbs)
    : MORTAR::Element(id, owner, shape, numnode, nodeids, isslave, isnurbs)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Element::Element(const CONTACT::Element& old) : MORTAR::Element(old)
{
  // empty copy-constructor

  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      mwgee 10/07|
 *----------------------------------------------------------------------*/
DRT::Element* CONTACT::Element::Clone() const
{
  CONTACT::Element* newele = new CONTACT::Element(*this);
  return newele;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::Element& element)
{
  element.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Element::Print(std::ostream& os) const
{
  os << "Contact ";
  MORTAR::Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Element::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class MORTAR::Element
  MORTAR::Element::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Element::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class MORTAR::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MORTAR::Element::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  number of dofs per node (public)                         mwgee 10/07|
 *----------------------------------------------------------------------*/
int CONTACT::Element::NumDofPerNode(const DRT::Node& node) const
{
  const CONTACT::Node* cnode = dynamic_cast<const CONTACT::Node*>(&node);
  if (!cnode) dserror("Node is not a Node");
  return cnode->NumDof();
}

/*----------------------------------------------------------------------*
 |  evaluate element (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
int CONTACT::Element::Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
    std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3)
{
  dserror("CONTACT::Element::Evaluate not implemented!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  Build element normal derivative at node                   popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Element::DerivNormalAtNode(int nid, int& i, CORE::LINALG::SerialDenseMatrix& elens,
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn)
{
  // find this node in my list of nodes and get local numbering
  int lid = GetLocalNodeId(nid);

  // get local coordinates for this node
  double xi[2];
  LocalCoordinatesOfNode(lid, xi);

  // build normal derivative at xi and return it
  DerivNormalAtXi(xi, i, elens, derivn);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal derivative at loc. coord. xi       popp 09/08|
 *----------------------------------------------------------------------*/
void CONTACT::Element::DerivNormalAtXi(double* xi, int& i, CORE::LINALG::SerialDenseMatrix& elens,
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn)
{
  // initialize variables
  const int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("DerivNormalAtXi: Null pointer!");
  CORE::LINALG::SerialDenseVector val(nnodes);
  CORE::LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  double gxi[3];
  double geta[3];

  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes);

  // get local element basis vectors
  Metrics(xi, gxi, geta);

  // derivative weighting matrix for current element
  CORE::LINALG::Matrix<3, 3> W;
  const double lcubeinv = 1.0 / (elens(4, i) * elens(4, i) * elens(4, i));

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      W(j, k) = -lcubeinv * elens(j, i) * elens(k, i);
      if (j == k) W(j, k) += 1 / elens(4, i);
    }
  }

  // now loop over all element nodes for derivatives
  for (int n = 0; n < nnodes; ++n)
  {
    Node* mycnode = dynamic_cast<Node*>(mynodes[n]);
    if (!mycnode) dserror("DerivNormalAtXi: Null pointer!");
    int ndof = mycnode->NumDof();

    // derivative weighting matrix for current node
    static CORE::LINALG::Matrix<3, 3> F;
    F(0, 0) = 0.0;
    F(1, 1) = 0.0;
    F(2, 2) = 0.0;
    F(0, 1) = geta[2] * deriv(n, 0) - gxi[2] * deriv(n, 1);
    F(0, 2) = gxi[1] * deriv(n, 1) - geta[1] * deriv(n, 0);
    F(1, 0) = gxi[2] * deriv(n, 1) - geta[2] * deriv(n, 0);
    F(1, 2) = geta[0] * deriv(n, 0) - gxi[0] * deriv(n, 1);
    F(2, 0) = geta[1] * deriv(n, 0) - gxi[1] * deriv(n, 1);
    F(2, 1) = gxi[0] * deriv(n, 1) - geta[0] * deriv(n, 0);

    // total weighting matrix
    static CORE::LINALG::Matrix<3, 3> WF;
    WF.MultiplyNN(W, F);

    // create directional derivatives
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < ndof; ++k) (derivn[j])[mycnode->Dofs()[k]] += WF(j, k) * NormalFac();
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal of last time step at xi          seitz 05/17 |
 *----------------------------------------------------------------------*/
void CONTACT::Element::OldUnitNormalAtXi(
    const double* xi, CORE::LINALG::Matrix<3, 1>& n_old, CORE::LINALG::Matrix<3, 2>& d_n_old_dxi)
{
  const int nnodes = NumNode();
  CORE::LINALG::SerialDenseVector val(nnodes);
  CORE::LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes);

  n_old.Clear();
  d_n_old_dxi.Clear();

  CORE::LINALG::Matrix<3, 1> tmp_n;
  CORE::LINALG::Matrix<3, 2> tmp_n_deriv;
  for (int i = 0; i < nnodes; ++i)
  {
    Node* cnode = dynamic_cast<CONTACT::Node*>(Nodes()[i]);
    if (!cnode) dserror("this is not a FriNode!");

    for (int d = 0; d < Dim(); ++d)
    {
      if (CORE::LINALG::Matrix<3, 1>(cnode->Data().Normal_old(), true).Norm2() < 0.9)
        dserror("where's my old normal");
      tmp_n(d) += val(i) * cnode->Data().Normal_old()[d];
      for (int x = 0; x < Dim() - 1; ++x)
        tmp_n_deriv(d, x) += deriv(i, x) * cnode->Data().Normal_old()[d];
    }
  }
  const double l = tmp_n.Norm2();
  n_old.Update(1. / l, tmp_n, 0.);

  CORE::LINALG::Matrix<2, 1> dli_dxi;
  dli_dxi.MultiplyTN(-1. / (l * l * l), tmp_n_deriv, tmp_n, 0.);
  d_n_old_dxi.Update(1. / l, tmp_n_deriv, 0.);
  d_n_old_dxi.MultiplyNT(1., tmp_n, dli_dxi, 1.);
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative J,xi of Jacobian determinant          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Element::DJacDXi(
    double* djacdxi, double* xi, const CORE::LINALG::SerialDenseMatrix& secderiv)
{
  // the derivative dJacdXi
  djacdxi[0] = 0.0;
  djacdxi[1] = 0.0;
  CORE::FE::CellType dt = Shape();

  // 2D linear case (2noded line element)
  // 3D linear case (3noded triangular element)
  if (dt == CORE::FE::CellType::line2 || dt == CORE::FE::CellType::tri3)
  {
    // do nothing
  }

  // 2D quadratic case (3noded line element)
  else if (dt == CORE::FE::CellType::line3 || dt == CORE::FE::CellType::nurbs2 ||
           dt == CORE::FE::CellType::nurbs3)
  {
    // get nodal coords for 2nd deriv. evaluation
    CORE::LINALG::SerialDenseMatrix coord(3, NumNode());
    GetNodalCoords(coord);

    // metrics routine gives local basis vectors
    double gxi[3];
    double geta[3];
    Metrics(xi, gxi, geta);

    std::array<double, 3> gsec = {0.0, 0.0, 0.0};
    for (int i = 0; i < NumNode(); ++i)
      for (int k = 0; k < 3; ++k) gsec[k] += secderiv(i, 0) * coord(k, i);

    // the Jacobian itself
    const double jacinv = 1.0 / sqrt(gxi[0] * gxi[0] + gxi[1] * gxi[1] + gxi[2] * gxi[2]);

    // compute dJacdXi (1 component in 2D)
    for (int dim = 0; dim < 3; ++dim) djacdxi[0] += gxi[dim] * gsec[dim] * jacinv;
  }

  // 3D bilinear case    (4noded quadrilateral element)
  // 3D quadratic case   (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::tri6 ||
           dt == CORE::FE::CellType::quad8 || dt == CORE::FE::CellType::quad9 ||
           dt == CORE::FE::CellType::nurbs4 || dt == CORE::FE::CellType::nurbs8 ||
           dt == CORE::FE::CellType::nurbs9)
  {
    // get nodal coords for 2nd deriv. evaluation
    CORE::LINALG::SerialDenseMatrix coord(3, NumNode());
    GetNodalCoords(coord);

    // metrics routine gives local basis vectors
    double gxi[3];
    double geta[3];
    Metrics(xi, gxi, geta);

    // cross product of gxi and geta
    std::array<double, 3> cross = {0.0, 0.0, 0.0};
    cross[0] = gxi[1] * geta[2] - gxi[2] * geta[1];
    cross[1] = gxi[2] * geta[0] - gxi[0] * geta[2];
    cross[2] = gxi[0] * geta[1] - gxi[1] * geta[0];

    // the Jacobian itself
    const double jacinv =
        1.0 / sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

    // 2nd deriv. evaluation
    CORE::LINALG::Matrix<3, 3> gsec(true);
    for (int i = 0; i < NumNode(); ++i)
      for (int k = 0; k < 3; ++k)
        for (int d = 0; d < 3; ++d) gsec(k, d) += secderiv(i, d) * coord(k, i);

    // compute dJacdXi (2 components in 3D)
    djacdxi[0] += jacinv * (cross[2] * geta[1] - cross[1] * geta[2]) * gsec(0, 0);
    djacdxi[0] += jacinv * (cross[0] * geta[2] - cross[2] * geta[0]) * gsec(1, 0);
    djacdxi[0] += jacinv * (cross[1] * geta[0] - cross[0] * geta[1]) * gsec(2, 0);
    djacdxi[0] += jacinv * (cross[1] * gxi[2] - cross[2] * gxi[1]) * gsec(0, 2);
    djacdxi[0] += jacinv * (cross[2] * gxi[0] - cross[0] * gxi[2]) * gsec(1, 2);
    djacdxi[0] += jacinv * (cross[0] * gxi[1] - cross[1] * gxi[0]) * gsec(2, 2);
    djacdxi[1] += jacinv * (cross[2] * geta[1] - cross[1] * geta[2]) * gsec(0, 2);
    djacdxi[1] += jacinv * (cross[0] * geta[2] - cross[2] * geta[0]) * gsec(1, 2);
    djacdxi[1] += jacinv * (cross[1] * geta[0] - cross[0] * geta[1]) * gsec(2, 2);
    djacdxi[1] += jacinv * (cross[1] * gxi[2] - cross[2] * gxi[1]) * gsec(0, 1);
    djacdxi[1] += jacinv * (cross[2] * gxi[0] - cross[0] * gxi[2]) * gsec(1, 1);
    djacdxi[1] += jacinv * (cross[0] * gxi[1] - cross[1] * gxi[0]) * gsec(2, 1);
  }

  // unknown case
  else
    dserror("DJacDXi called for unknown element type!");

  return;
}


void CONTACT::Element::PrepareDderiv(const std::vector<MORTAR::Element*>& meles)
{
  // number of dofs that may appear in the linearization
  int numderiv = 0;
  numderiv += NumNode() * 3 * 12;
  for (unsigned m = 0; m < meles.size(); ++m) numderiv += (meles.at(m))->NumNode() * 3;
  dMatrixDeriv_ = Teuchos::rcp(new CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>(
      numderiv, 0, CORE::LINALG::SerialDenseMatrix(NumNode(), NumNode())));
}

void CONTACT::Element::PrepareMderiv(const std::vector<MORTAR::Element*>& meles, const int m)
{
  // number of dofs that may appear in the linearization
  int numderiv = 0;
  numderiv += NumNode() * 3 * 12;
  for (unsigned i = 0; i < meles.size(); ++i) numderiv += meles[i]->NumNode() * 3;
  mMatrixDeriv_ = Teuchos::rcp(new CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>(
      numderiv, 0, CORE::LINALG::SerialDenseMatrix(NumNode(), meles.at(m)->NumNode())));
}


void CONTACT::Element::AssembleDderivToNodes(bool dual)
{
  if (dMatrixDeriv_ == Teuchos::null)
    dserror("AssembleDderivToNodes called w/o PrepareDderiv first");

  if (dMatrixDeriv_->size() == 0) return;

  for (int j = 0; j < NumNode(); ++j)
  {
    CONTACT::Node* cnode_j = dynamic_cast<CONTACT::Node*>(Nodes()[j]);

    if (!dual)
    {
      for (int k = 0; k < NumNode(); ++k)
      {
        CONTACT::Node* cnode_k = dynamic_cast<CONTACT::Node*>(Nodes()[k]);
        std::map<int, double>& ddmap_jk = cnode_j->Data().GetDerivD()[cnode_k->Id()];

        for (CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>::const_iterator p =
                 dMatrixDeriv_->begin();
             p != dMatrixDeriv_->end(); ++p)
          ddmap_jk[p->first] += (p->second)(j, k);
      }
    }
    else
    {
      std::map<int, double>& ddmap_jj = cnode_j->Data().GetDerivD()[cnode_j->Id()];

      for (CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>::const_iterator p =
               dMatrixDeriv_->begin();
           p != dMatrixDeriv_->end(); ++p)
        ddmap_jj[p->first] += (p->second)(j, j);
    }
  }
  dMatrixDeriv_ = Teuchos::null;
}

void CONTACT::Element::AssembleMderivToNodes(MORTAR::Element& mele)
{
  if (mMatrixDeriv_ == Teuchos::null)
    dserror("AssembleMderivToNodes called w/o PrepareMderiv first");
  if (mMatrixDeriv_->size() == 0) return;

  for (int j = 0; j < NumNode(); ++j)
  {
    CONTACT::Node* cnode_j = dynamic_cast<CONTACT::Node*>(Nodes()[j]);

    for (int k = 0; k < mele.NumNode(); ++k)
    {
      CONTACT::Node* cnode_k = dynamic_cast<CONTACT::Node*>(mele.Nodes()[k]);
      std::map<int, double>& dmmap_jk = cnode_j->Data().GetDerivM()[cnode_k->Id()];

      for (CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>::const_iterator p =
               mMatrixDeriv_->begin();
           p != mMatrixDeriv_->end(); ++p)
        dmmap_jk[p->first] += (p->second)(j, k);
    }
  }
}

BACI_NAMESPACE_CLOSE
