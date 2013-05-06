/*!-----------------------------------------------------------------------*
 \file structporo.cpp

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/



#include <vector>
#include "structporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/so3_material.H"

#include "../drt_lib/drt_utils_factory.H"  // for function Factory in Unpack

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoro::StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  matid_(matdata->GetInt("MATID")),
  bulkmodulus_(matdata->GetDouble("BULKMODULUS")),
  penaltyparameter_(matdata->GetDouble("PENALTYPARAMETER")),
  initporosity_(matdata->GetDouble("INITPOROSITY"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoro(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroType MAT::StructPoroType::instance_;

DRT::ParObject* MAT::StructPoroType::Create(const std::vector<char> & data)
{
  MAT::StructPoro* struct_poro = new MAT::StructPoro();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro() :
  params_(NULL),
  mat_(Teuchos::null),
  porosity_(Teuchos::null),
  surfporosity_(Teuchos::null),
  //dporodt_(Teuchos::null),
  //gradporosity_(Teuchos::null),
  isinitialized_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params) :
  params_(params),
  porosity_(Teuchos::null),
  surfporosity_(Teuchos::null),
  //dporodt_(Teuchos::null),
  //gradporosity_(Teuchos::null),
  isinitialized_(false)
{
  mat_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->matid_));
  if (mat_ == Teuchos::null) dserror("MAT::StructPoro: underlying material should be of type MAT::So3Material");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::PoroSetup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  porosity_ = Teuchos::rcp(new std::vector< double > (numgp,params_->initporosity_));
  surfporosity_ = Teuchos::rcp(new std::map<int, std::vector< double > >);

  /*
  gradporosity_ = Teuchos::rcp(new std::vector< LINALG::Matrix<3,1> > (numgp));
  dporodt_ = Teuchos::rcp(new std::vector<double> (numgp));
   gradJ_ = Teuchos::rcp(new std::vector< LINALG::Matrix<1,3> > (numgp));

  const LINALG::Matrix<3,1> emptyvec(true);
   const LINALG::Matrix<1,3> emptyvecT(true);
  for (int j=0; j<numgp; ++j)
  {
     gradporosity_->at(j) = emptyvec;
    gradJ_->at(j) = emptyvecT;
  }
  */

  isinitialized_=true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Pack(DRT::PackBuffer& data) const
{
  if(not isinitialized_)
    dserror("poro material not initialized. Not a poro element?");

  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);

  // porosity_
  int size=0;
  size = (int)porosity_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*porosity_)[i]);
  }

  // refporosity_
  //AddtoPack(data, refporosity_);

  // surfporosity_ (i think it is not necessary to pack/unpack this...)
  size = (int) surfporosity_->size();
  AddtoPack(data,size);
  // iterator
  std::map<int,std::vector<double> >::const_iterator iter;
  for(iter=surfporosity_->begin();iter!=surfporosity_->end();++iter)
  {
    AddtoPack(data,iter->first);
    AddtoPack(data,iter->second);
  }

  /*
  // gradporosity_
  size = (int)gradporosity_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*gradporosity_)[i]);
  }

  // dporodt_
  size = (int)dporodt_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*dporodt_)[i]);
  }

  // gradJ_
  size = (int)gradJ_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*gradJ_)[i]);
  }
  */

  // Pack data of underlying material
  if (mat_!=Teuchos::null)
    mat_->Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::StructPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // porosity_
  int size = 0;
  ExtractfromPack(position,data,size);
  porosity_=Teuchos::rcp(new std::vector<double >);
  double tmp = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp);
    porosity_->push_back(tmp);
  }

  // refporosity_
  //ExtractfromPack(position,data,refporosity_);

  // surface porosity (i think it is not necessary to pack/unpack this...)
  ExtractfromPack(position,data,size);
  surfporosity_ = Teuchos::rcp(new std::map<int, std::vector< double > >);
  for(int i=0;i<size;i++)
  {
    int dof;
    std::vector<double > value;
    ExtractfromPack(position,data,dof);
    ExtractfromPack(position,data,value);

    //add to map
    surfporosity_->insert(std::pair<int,std::vector<double > >(dof,value));
  }

  /*
  // gradporosity_
  gradporosity_=Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >);
  ExtractfromPack(position,data,size);
  LINALG::Matrix<3,1> tmp2(true);
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp2);
    gradporosity_->push_back(tmp2);
  }

  // dporodt_
  ExtractfromPack(position,data,size);
  dporodt_=Teuchos::rcp(new std::vector<double >);
  tmp = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp);
    dporodt_->push_back(tmp);
  }

  // gradJ_
  gradJ_=Teuchos::rcp(new std::vector<LINALG::Matrix<1,3> >);
  ExtractfromPack(position,data,size);
  LINALG::Matrix<1,3> tmp3(true);
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp3);
    gradJ_->push_back(tmp3);
  }
  */

  // Unpack data of sub material (these lines are copied from drt_element.cpp)
  std::vector<char> datamat;
  ExtractfromPack(position,data,datamat);
  if (datamat.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(datamat);  // Unpack is done here
    MAT::So3Material* mat = dynamic_cast<MAT::So3Material*>(o);
    if (mat==NULL)
      dserror("failed to unpack elastic material");
    mat_ = Teuchos::rcp(mat);
  }
  else mat_ = Teuchos::null;

  isinitialized_=true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( const double& initporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref,
                                       bool save)
{

  const double & bulkmodulus  = params_->bulkmodulus_;
  const double & penalty      = params_->penaltyparameter_;

  const double a = (bulkmodulus / (1 - initporosity) + press - penalty / initporosity) * J;
  const double b = -a + bulkmodulus + penalty;
  const double c = b * b  + 4.0 * penalty * a;
  double d = sqrt(c);

 // const double c = (b / a) * (b / a) + 4.0 * penalty / a;
 // double d = sqrt(c) * a;

  double test = 1 / (2.0 * a) * (-b + d);
  double sign = 1.0;
  if (test >= 1.0 or test < 0.0)
  {
    sign = -1.0;
    d = sign * d;
  }

  const double phi = 1 / (2 * a) * (-b + d);

  if (phi >= 1.0 or phi < 0.0)
    dserror("invalid porosity: %f", porosity);

  const double d_p = J * (-b+2.0*penalty)/d;
  const double d_p_p = ( d * J + d_p * (b - 2.0*penalty) ) / (d * d) * J;
  const double d_J = a/J * ( -b + 2.0*penalty ) / d;
  const double d_J_p = (d_p / J + ( 1-d_p*d_p/(J*J) ) / d *a);
  const double d_J_J = ( a*a/(J*J)-d_J*d_J )/ d;

  //d(porosity) / d(p)
  if(dphi_dp) *dphi_dp = - J * phi/a + (J+d_p)/(2.0*a);

  //d(porosity) / d(J)
  if(dphi_dJ) *dphi_dJ= -phi/J+ 1/(2*J) + d_J / (2.0*a);

  //d(porosity) / d(J)d(pressure)
  if(dphi_dJdp) *dphi_dJdp= -1/J* (*dphi_dp)+ d_J_p/(2*a) - d_J*J/(2.0*a*a);

  //d^2(porosity) / d(J)^2
  if(dphi_dJJ) *dphi_dJJ= phi/(J*J) - (*dphi_dJ)/J - 1/(2.0*J*J) - d_J/(2*a*J) + d_J_J/(2.0*a);

  //d^2(porosity) / d(pressure)^2
  if(dphi_dpp) *dphi_dpp= -J/a* (*dphi_dp) + phi*J*J/(a*a) - J/(2.0*a*a)*(J+d_p) + d_p_p/(2.0*a);

  /*
  dphi_dp = 0.0;
  dphi_dJ = 0.0;
  dphi_dJdp = 0.0;
  dphi_dJJ = 0.0;
  dphi_dpp = 0.0;
  phi = 0.1;
  */

  /*
  double phi = initporosity + (1.0-initporosity) * (1.0-initporosity) / bulkmodulus * press
                + (1.0-initporosity) * (J-1.0);

  dphi_dp = (1.0-initporosity) * (1.0-initporosity) / bulkmodulus;
  dphi_dJ = (1.0-initporosity);
  dphi_dJdp = 0.0;
  dphi_dJJ = 0.0;
  dphi_dpp = 0.0;
  */

  porosity= phi;

  //save porosity
  if(save)
    porosity_->at(gp) = phi;

  if(dphi_dphiref)
  {
    const double dadphiref = J*(bulkmodulus / ((1 - initporosity)*(1 - initporosity)) + penalty / (initporosity*initporosity));
    const double tmp = 2*dadphiref/a * (-b*(a+b)/a - 2*penalty);
    const double dddphiref = sign*(dadphiref * sqrt(c)/a + tmp);

    *dphi_dphiref = ( a * (dadphiref+dddphiref) - dadphiref * (-b + d) )/(2*a*a);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( Teuchos::ParameterList& params,
                                       double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       bool save)
{

  ComputePorosity( params_->initporosity_,
                   press,
                   J,
                   gp,
                   porosity,
                   dphi_dp,
                   dphi_dJ,
                   dphi_dJdp,
                   dphi_dJJ,
                   dphi_dpp,
                   NULL,
                   save);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( Teuchos::ParameterList& params,
                                       double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       bool save)
{

  ComputePorosity( params,
                   press,
                   J,
                   gp,
                   porosity,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   save);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputeSurfPorosity( Teuchos::ParameterList& params,
                                           double     press,
                                           double     J,
                                           const int  surfnum,
                                           int        gp,
                                           double&    porosity,
                                           double*    dphi_dp,
                                           double*    dphi_dJ,
                                           double*    dphi_dJdp,
                                           double*    dphi_dJJ,
                                           double*    dphi_dpp,
                                           bool save)
{
  ComputePorosity(params,
                  press,
                  J,
                  gp,
                  porosity,
                  dphi_dp,
                  dphi_dJ,
                  dphi_dJdp,
                  dphi_dJJ,
                  dphi_dpp,
                  save);

  if(gp==0)  //it's a new iteration, so old values are not needed any more
   ( (*surfporosity_)[surfnum] ).clear();

  ( (*surfporosity_)[surfnum] ).push_back(porosity);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputeSurfPorosity( Teuchos::ParameterList& params,
                                           double     press,
                                           double     J,
                                           const int  surfnum,
                                           int        gp,
                                           double&    porosity,
                                           bool save)
{

  ComputeSurfPorosity( params,
                       press,
                       J,
                       surfnum,
                       gp,
                       porosity,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       save);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::StructPoro::PorosityAv() const
{
  double porosityav = 0.0;

  std::vector<double>::const_iterator m;
  for (m = porosity_->begin(); m != porosity_->end(); ++m)
  {
    porosityav += *m;
  }
  porosityav = porosityav / (porosity_->size());

  return porosityav;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::CouplStress(  const LINALG::Matrix<3,3>& defgrd,
                                    const LINALG::Matrix<3,1>& fluidvel,
                                    const double& press,
                                    LINALG::Matrix<6,1>& couplstress
                                    ) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  LINALG::Matrix<3,3> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // inverse Right Cauchy-Green tensor
  LINALG::Matrix<3,3> C_inv;
  C_inv.Invert(cauchygreen);

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<6,1> C_inv_vec;
  for(int i =0, k=0;i<3; i++)
    for(int j =0;j<3-i; j++,k++)
      C_inv_vec(k)=C_inv(i+j,j);

  for(int i=0; i<6 ; i++)
    couplstress(i)= -1.0*J*press*C_inv_vec(i);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::CouplStress(  const LINALG::Matrix<2,2>& defgrd,
                                    const LINALG::Matrix<2,1>& fluidvel,
                                    const double& press,
                                    LINALG::Matrix<3,1>& couplstress
                                    ) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  LINALG::Matrix<2,2> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // inverse Right Cauchy-Green tensor
  LINALG::Matrix<2,2> C_inv;
  C_inv.Invert(cauchygreen);

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<3,1> C_inv_vec;
  for(int i =0, k=0;i<2; i++)
    for(int j =0;j<2-i; j++,k++)
      C_inv_vec(k)=C_inv(i+j,j);

  for(int i=0; i<3 ; i++)
    couplstress(i)= -1.0*J*press*C_inv_vec(i);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::PorosityGradientAv(LINALG::Matrix<3,1>& porosityav) const
{
  porosityav.Clear();

  std::vector<LINALG::Matrix<3,1> >::const_iterator m;
  for (m = gradporosity_->begin(); m != gradporosity_->end(); ++m)
  {
    porosityav.Update(1.0, *m ,1.0) ;
  }
  double size = gradporosity_->size();

  porosityav.Scale(1.0/size);

  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::PorosityGradientAv(LINALG::Matrix<2,1>& porositygradientav) const
{
  porositygradientav.Clear();

  const int size = gradporosity_->size();
  for(int i=0; i<size;i++)
  {
    porositygradientav(0) += (*gradporosity_)[i](0);
    porositygradientav(1) += (*gradporosity_)[i](1);
  }
  porositygradientav.Scale(1.0/size);
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
double MAT::StructPoro::DPorosityDtAv() const
{
  double av = 0.0;

  std::vector<double>::const_iterator m;
  for (m = dporodt_->begin(); m != dporodt_->end(); ++m)
  {
    av += *m;
  }
  av = av / (dporodt_->size());

  return av;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetDPoroDtAtGP(std::vector<double> dporodt_gp)
{
  //set dporodt values
  std::vector<double>::iterator m = dporodt_gp.begin();
  for (int i = 0; m != dporodt_gp.end(); ++m, ++i)
  {
    double dporodt = *m;
    dporodt_->at(i) = dporodt;
  }
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetGradPorosityAtGP(std::vector<  LINALG::Matrix<3,1> > gradporosity_gp)
{
  std::vector<LINALG::Matrix<3,1> >::iterator m = gradporosity_gp.begin();
  for (int i = 0; m != gradporosity_gp.end(); ++m, ++i)
  {
    LINALG::Matrix<3,1> gradporo  = *m;
    gradporosity_->at(i) = gradporo;
  }
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetGradPorosityAtGP(std::vector<  LINALG::Matrix<2,1> > gradporosity_gp)
{
  for(unsigned int i=0; i<gradporosity_gp.size();i++)
  {
    ((*gradporosity_)[i])(0) = (gradporosity_gp[i])(0);
    ((*gradporosity_)[i])(1) = (gradporosity_gp[i])(1);
  }
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetGradJAtGP(std::vector<  LINALG::Matrix<1,3> > gradJ_gp)
{
  std::vector<LINALG::Matrix<1,3> >::iterator m = gradJ_gp.begin();
  for (int i = 0; m != gradJ_gp.end(); ++m, ++i)
  {
    LINALG::Matrix<1,3> gradJ  = *m;
    (*gradJ_)[i] = gradJ;
  }
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetGradJAtGP(std::vector<  LINALG::Matrix<1,2> > gradJ_gp)
{
  for(unsigned int i=0; i<gradJ_gp.size();i++)
  {
    ((*gradJ_)[i])(0) = (gradJ_gp[i])(0);
    ((*gradJ_)[i])(1) = (gradJ_gp[i])(1);
  }
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::GetGradJAtGP(LINALG::Matrix<1,3>& gradJ, int gp) const
{
  gradJ = (*gradJ_)[gp];
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::GetGradJAtGP(LINALG::Matrix<1,2>& gradJ, int gp) const
{
  gradJ(0) = ((*gradJ_)[gp])(0);
  gradJ(1) = ((*gradJ_)[gp])(1);

  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::GetGradPorosityAtGP(LINALG::Matrix<3,1>& gradporosity, int gp) const
{
  gradporosity = gradporosity_->at(gp);
  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::GetGradPorosityAtGP(LINALG::Matrix<2,1>& gradporosity, int gp) const
{
  gradporosity(0) = ((*gradporosity_)[gp])(0);
  gradporosity(1) = ((*gradporosity_)[gp])(1);

  return;
}
*/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::ConsitutiveDerivatives(Teuchos::ParameterList& params,
                                              double     press,
                                              double     J,
                                              double     porosity,
                                              double*    dW_dp,
                                              double*    dW_dphi,
                                              double*    dW_dJ,
                                              double*    W)
{
  if(porosity == 0.0)
    dserror("porosity equals zero!! Wrong initial porosity?");
  const double & bulkmodulus  = params_->bulkmodulus_;
  const double & penalty      = params_->penaltyparameter_;
  const double & initporosity = params_->initporosity_;

  //some intermediate values
  const double a = bulkmodulus / (1 - initporosity) + press - penalty / initporosity;
  const double b = -1.0*J*a+bulkmodulus+penalty;

  if(W)       *W       = J*a*porosity*porosity + porosity* b - penalty;
  if(dW_dp)   *dW_dp   = -1.0*J*porosity *(1.0-porosity);
  if(dW_dphi) *dW_dphi = 2.0*J*a*porosity + b;
  if(dW_dJ)   *dW_dJ   = a*porosity*porosity - porosity*a;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::VisNames(std::map<string,int>& names)
{
  mat_->VisNames(names);
  std::string porosity = "porosity";
  names[porosity] = 1; // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool MAT::StructPoro::VisData(const string& name, std::vector<double>& data, int numgp)
{
  if (mat_->VisData(name,data,numgp))
    return true;
  if (name=="porosity")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    data[0] = PorosityAv();
    return true;
  }
  return false;
}
