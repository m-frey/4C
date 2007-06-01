/*!----------------------------------------------------------------------
\file wall1_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
}
#include "../drt_lib/dstrc.H"

/*----------------------------------------------------------------------*
 |                                                        mgit 03/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Wall1::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("Wall1::Evaluate");
  DRT::Elements::Wall1::ActionType act = Wall1::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Wall1::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Wall1::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Wall1::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Wall1::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Wall1::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Wall1::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Wall1::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Wall1::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Wall1::calc_struct_update_istep;
  else dserror("Unknown type of action for Wall1");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);
  switch(act)
  {
    case Wall1::calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    case Wall1::calc_struct_nlnstiffmass:
    case Wall1::calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    default:
      dserror("Unknown type of action for Wall1 %d", act);
  }
  return 0;

}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mgit 05/07|
 *----------------------------------------------------------------------*/

int DRT::Elements::Wall1::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::Utils::ExtractMyValues(*disp,mydisp,lm);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // no. of nodes on this surface
  const int iel = NumNode();

 // general arrays
  int       ngauss  = 0;
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  double det;
  // init gaussian points
  W1_DATA w1data;
  w1_integration_points(w1data);

  const int nir = ngp_[0];
  const int nis = ngp_[1];
  const int numdf = 2;
//  vector<double>* thick = data_.Get<vector<double> >("thick");
//  if (!thick) dserror("Cannot find vector of nodal thickness");

  vector<double> funct(iel);
  Epetra_SerialDenseMatrix deriv(2,iel);

  double xrefe[2][MAXNOD_WALL1];
  double xcure[2][MAXNOD_WALL1];

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {

    xrefe[0][k] = Nodes()[k]->X()[0];
    xrefe[1][k] = Nodes()[k]->X()[1];

    xcure[0][k] = xrefe[0][k] + mydisp[k*numdf+0];
    xcure[1][k] = xrefe[1][k] + mydisp[k*numdf+1];

  }


  // get values and switches from the condition
  vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  vector<double>* val   = condition.Get<vector<double> >("val");
  /*=================================================== integration loops */
  for (int lr=0; lr<nir; ++lr)
  {
    /*================================== gaussian point and weight at it */
    const double e1   = w1data.xgrr[lr];
    double facr = w1data.wgtr[lr];
      for (int ls=0; ls<nis; ++ls)
    {
      const double e2   = w1data.xgss[ls];
      double facs = w1data.wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      w1_shapefunctions(funct,deriv,e1,e2,iel,1);
      /*--------------------------------------- compute jacobian Matrix */
      w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);
      /*------------------------------------ integration factor  -------*/
     double fac=0;
     fac = facr * facs * det;

    // load vector ar
    double ar[2];
    // loop the dofs of a node
    // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]
    for (int i=0; i<2; ++i)
   {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;

    }

    // add load components
    for (int node=0; node<NumNode(); ++node)
      for (int dof=0; dof<2; ++dof)
         elevec1[node*2+dof] += funct[node] *ar[dof];

      ngauss++;
    } // for (int ls=0; ls<nis; ++ls)
  } // for (int lr=0; lr<nir; ++lr)


cout << elevec1;
return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_nlnstiffmass(vector<int>&               lm,
                                            vector<double>&           disp,
                                            vector<double>&           residual,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            struct _MATERIAL*         material)
{
  DSTraceHelper dst("Wall1::w1_nlnstiffmass");
  const int numnode = NumNode();
  const int numdf   = 2;
  int       ngauss  = 0;
  const int nd      = numnode*numdf;


   // general arrays
  vector<double>           funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix boplin;
  boplin.Shape(4,2*numnode);
  Epetra_SerialDenseVector F;
  F.Size(4);
  Epetra_SerialDenseVector strain;
  strain.Size(4);
  double det;
  double xrefe[2][MAXNOD_WALL1];
  double xcure[2][MAXNOD_WALL1];
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);

  // gaussian points
  W1_DATA w1data;
  w1_integration_points(w1data);

  // ------------------------------------ check calculation of mass matrix
  int imass=0;
  double density=0.0;
  if (massmatrix)
  {
    imass=1;
    w1_getdensity(material,&density);
  }

  const int nir = ngp_[0];
  const int nis = ngp_[1];
  const int iel = numnode;

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {

    xrefe[0][k] = Nodes()[k]->X()[0];
    xrefe[1][k] = Nodes()[k]->X()[1];

    xcure[0][k] = xrefe[0][k] + disp[k*numdf+0];
    xcure[1][k] = xrefe[1][k] + disp[k*numdf+1];

  }

  /*=================================================== integration loops */
  for (int lr=0; lr<nir; ++lr)
  {
    /*================================== gaussian point and weight at it */
    const double e1   = w1data.xgrr[lr];
    double facr = w1data.wgtr[lr];
      for (int ls=0; ls<nis; ++ls)
    {
      const double e2   = w1data.xgss[ls];
      double facs = w1data.wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      w1_shapefunctions(funct,deriv,e1,e2,iel,1);
      /*--------------------------------------- compute jacobian Matrix */
      w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);
      /*------------------------------------ integration factor  -------*/
     double fac=0;
     fac = facr * facs * det * thickness_;

      /*------------------------------compute mass matrix if imass-----*/
      if (imass)
      {
       double facm = fac * density;
       for (int a=0; a<iel; a++)
       {
        for (int b=0; b<iel; b++)
        {
         (*massmatrix)(2*a,2*b)     += facm * funct[a] * funct[b]; /* a,b even */
         (*massmatrix)(2*a+1,2*b+1) += facm * funct[a] * funct[b]; /* a,b odd  */
        }
       }
      }
     /*----------------------------------- calculate operator Blin  ---*/
     w1_boplin(boplin,deriv,xjm,det,iel);

     /*----------------- calculate defgrad F, Green-Lagrange-strain ---*/
     w1_defgrad(F,strain,xrefe,xcure,boplin,iel);

     /*-calculate defgrad F in matrix notation and Blin in curent conf.*/
     w1_boplin_cure(b_cure,boplin,F,numeps,nd);
     /*------------------------------------------ call material law ---*/
     w1_call_matgeononl(strain,stress,C,numeps,material);
     /*---------------------- geometric part of stiffness matrix kg ---*/
     w1_kg(*stiffmatrix,boplin,stress,fac,nd,numeps);
     /*------------------ elastic+displacement stiffness matrix keu ---*/
     w1_keu(*stiffmatrix,b_cure,C,fac,nd,numeps);
     /*--------------- nodal forces fi from integration of stresses ---*/
       if (force) w1_fint(stress,b_cure,*force,fac,nd);

      ngauss++;
    } // for (int ls=0; ls<nis; ++ls)
  } // for (int lr=0; lr<nir; ++lr)

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element integration points (private)        mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_integration_points(struct _W1_DATA& data)
{
  DSTraceHelper dst("Wall1::w1_integration_points");

  const int numnode = NumNode();

  const double invsqrtthree = 1./sqrt(3.);
  const double sqrtthreeinvfive = sqrt(3./5.);
  const double wgt  = 5.0/9.0;
  const double wgt0 = 8.0/9.0;


  // quad elements
  if (numnode==4 || numnode==8 || numnode==9)
  {
    switch(ngp_[0]) // r direction
    {
      case 1:
        data.xgrr[0] = 0.0;
        data.xgrr[1] = 0.0;
        data.xgrr[2] = 0.0;
        data.wgtr[0] = 2.0;
        data.wgtr[1] = 0.0;
        data.wgtr[2] = 0.0;
      break;
      case 2:
        data.xgrr[0] = -invsqrtthree;
        data.xgrr[1] =  invsqrtthree;
        data.xgrr[2] =  0.0;
        data.wgtr[0] =  1.0;
        data.wgtr[1] =  1.0;
        data.wgtr[2] =  0.0;
      break;
      case 3:
        data.xgrr[0] = -sqrtthreeinvfive;
        data.xgrr[1] =  0.0;
        data.xgrr[2] =  sqrtthreeinvfive;
        data.wgtr[0] =  wgt;
        data.wgtr[1] =  wgt0;
        data.wgtr[2] =  wgt;
      break;
      default:
        dserror("Unknown no. of gaussian points in r-direction");
      break;
    } // switch(ngp_[0]) // r direction

    switch(ngp_[1]) // s direction
    {
      case 1:
        data.xgss[0] = 0.0;
        data.xgss[1] = 0.0;
        data.xgss[2] = 0.0;
        data.wgts[0] = 2.0;
        data.wgts[1] = 0.0;
        data.wgts[2] = 0.0;
      break;
      case 2:
        data.xgss[0] = -invsqrtthree;
        data.xgss[1] =  invsqrtthree;
        data.xgss[2] =  0.0;
        data.wgts[0] =  1.0;
        data.wgts[1] =  1.0;
        data.wgts[2] =  0.0;
      break;
      case 3:
        data.xgss[0] = -sqrtthreeinvfive;
        data.xgss[1] =  0.0;
        data.xgss[2] =  sqrtthreeinvfive;
        data.wgts[0] =  wgt;
        data.wgts[1] =  wgt0;
        data.wgts[2] =  wgt;
      break;
      default:
        dserror("Unknown no. of gaussian points in s-direction");
      break;
    } // switch(ngp_[0]) // s direction

  } // if (numnode==4 || numnode==8 || numnode==9)

//  else if (numnode==3 || numnode==6) // triangle elements
//  {
//    switch(ngptri_)
//    {
//      case 1:
//      {
//        const double third = 1.0/3.0;
//        data.xgrr[0] =  third;
//        data.xgrr[1] =  0.0;
//        data.xgrr[2] =  0.0;
//        data.xgss[0] =  third;
//        data.xgss[1] =  0.0;
//        data.xgss[2] =  0.0;
//        data.wgtr[0] =  0.5;
//        data.wgtr[1] =  0.0;
//        data.wgtr[2] =  0.0;
//        data.wgts[0] =  0.5;
//        data.wgts[1] =  0.0;
//        data.wgts[2] =  0.0;
//      }
//      break;
//      case 3:
//      {
//        const double wgt = 1.0/6.0;
//        data.xgrr[0] =  0.5;
//        data.xgrr[1] =  0.5;
//        data.xgrr[2] =  0.0;
//        data.xgss[0] =  0.0;
//        data.xgss[1] =  0.5;
//        data.xgss[2] =  0.5;
//        data.wgtr[0] =  wgt;
//        data.wgtr[1] =  wgt;
//        data.wgtr[2] =  wgt;
//        data.wgts[0] =  wgt;
//        data.wgts[1] =  wgt;
//        data.wgts[2] =  wgt;
//      }
//      break;
//      default:
//        dserror("Unknown no. of gaussian points for triangle");
//      break;
//    }
//  } // else if (numnode==3 || numnode==6)

  return;
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives (private)                mgit 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_shapefunctions(
                             vector<double>& funct,
                             Epetra_SerialDenseMatrix& deriv,
                             const double r, const double s, const int numnode,
                             const int doderiv) const
{
  DSTraceHelper dst("Wall1::w1_shapefunctions");

  const double q12 = 0.5;
  const double q14 = 0.25;
  const double rr = r*r;
  const double ss = s*s;
  const double rp = 1.0+r;
  const double rm = 1.0-r;
  const double sp = 1.0+s;
  const double sm = 1.0-s;
  const double r2 = 1.0-rr;
  const double s2 = 1.0-ss;
  int i;
  int ii;

  //const double t;

  switch(numnode)
  {
    case 4:
    {
      funct[0] = q14*rp*sp;
      funct[1] = q14*rm*sp;
      funct[2] = q14*rm*sm;
      funct[3] = q14*rp*sm;
      if (doderiv)
      {
        deriv(0,0)= q14*sp;
        deriv(0,1)=-q14*sp;
        deriv(0,2)=-q14*sm;
        deriv(0,3)= q14*sm;
        deriv(1,0)= q14*rp;
        deriv(1,1)= q14*rm;
        deriv(1,2)=-q14*rm;
        deriv(1,3)=-q14*rp;
      }
      return;
    }
    break;
    case 8:
    {
      funct[0] = q14*rp*sp;
      funct[1] = q14*rm*sp;
      funct[2] = q14*rm*sm;
      funct[3] = q14*rp*sm;
      funct[4] = q12*r2*sp;
      funct[5] = q12*rm*s2;
      funct[6] = q12*r2*sm;
      funct[7] = q12*rp*s2;
      funct[0] = funct[0] - q12*(funct[4] + funct[7]);
      if (doderiv)
      {
         deriv(0,0)= q14*sp;
         deriv(0,1)=-q14*sp;
         deriv(0,2)=-q14*sm;
         deriv(0,3)= q14*sm;
         deriv(1,0)= q14*rp;
         deriv(1,1)= q14*rm;
         deriv(1,2)=-q14*rm;
         deriv(1,3)=-q14*rp;
         deriv(0,4)=-1*r*sp;
         deriv(0,5)=-q12*  s2;
         deriv(0,6)=-1*r*sm;
         deriv(0,7)= q12*  s2;
         deriv(1,4)= q12*r2  ;
         deriv(1,5)=-1*rm*s;
         deriv(1,6)=-q12*r2  ;
         deriv(1,7)=-1*rp*s;

         deriv[0][0]=deriv[0][0] - q12*(deriv[0][4] + deriv[0][7]);
         deriv[1][0]=deriv[1][0] - q12*(deriv[1][4] + deriv[1][7]);
      }
      for (i=1; i<=3; i++)
      {
         ii=i + 3;
         funct[i]=funct[i] - q12*(funct[ii] + funct[ii+1]);
         if (doderiv)              /*--- check for derivative evaluation ---*/
         {
             deriv(0,i)=deriv(0,i) - q12*(deriv(0,ii) + deriv(0,ii+1));
             deriv(1,i)=deriv(1,i) - q12*(deriv(1,ii) + deriv(1,ii+1));
         }
      }
      return;
    }
    break;
    case 9:
    {
      const double rh  = q12*r;
      const double sh  = q12*s;
      const double rs  = rh*sh;
      const double rhp = r+q12;
      const double rhm = r-q12;
      const double shp = s+q12;
      const double shm = s-q12;
      funct[0] = rs*rp*sp;
      funct[1] =-rs*rm*sp;
      funct[2] = rs*rm*sm;
      funct[3] =-rs*rp*sm;
      funct[4] = sh*sp*r2;
      funct[5] =-rh*rm*s2;
      funct[6] =-sh*sm*r2;
      funct[7] = rh*rp*s2;
      funct[8] = r2*s2;
      if (doderiv==1)
      {
         deriv(0,0)= rhp*sh*sp;
         deriv(0,1)= rhm*sh*sp;
         deriv(0,2)=-rhm*sh*sm;
         deriv(0,3)=-rhp*sh*sm;
         deriv(0,4)=-2.0*r*sh*sp;
         deriv(0,5)= rhm*s2;
         deriv(0,6)= 2.0*r*sh*sm;
         deriv(0,7)= rhp*s2;
         deriv(0,8)=-2.0*r*s2;
         deriv(1,0)= shp*rh*rp;
         deriv(1,1)=-shp*rh*rm;
         deriv(1,2)=-shm*rh*rm;
         deriv(1,3)= shm*rh*rp;
         deriv(1,4)= shp*r2;
         deriv(1,5)= 2.0*s*rh*rm;
         deriv(1,6)= shm*r2;
         deriv(1,7)=-2.0*s*rh*rp;
         deriv(1,8)=-2.0*s*r2;
      }
      return;
   }
   break;
//   case 3:
//   {
//     funct[0]=1-r-s;
//     funct[1]=r;
//     funct[2]=s;
//     if (doderiv==1)
//     {
//       deriv(0,0)=  -1.0;
//       deriv(0,1)=  1.0;
//       deriv(0,2)=  0.0;
//       deriv(1,0)=  -1.0;
//       deriv(1,1)=  0.0;
//       deriv(1,2)=  1.0;
//     }
//   }
//   break;
//    case 6:
//     funct[0]=(1-2*r-2*s)*(1-r-s);
//     funct[1]=2*r*r-r;
//     funct[2]=2*s*s-s;
//     funct[3]=4*(r-r*r-r*s);
//     funct[4]=4*r*s;
//     funct[5]=4*(s-s*s-s*r);
//     if (doderiv==1)
//     {
//       deriv(0,0)= -3.0+4.0*r+4.0*s;
//       deriv(0,1)= 4.0*r-1.0;
//       deriv(0,2)= 0.0;
//       deriv(0,3)= 4.0*(1-2.0*r-s);
//       deriv(0,4)= 4.0*s;
//       deriv(0,5)= -4.0*s;
//       deriv(1,0)= -3.0+4.0*r+4.0*s;
//       deriv(1,1)= 0.0;
//       deriv(1,2)= 4.0*s-1.0;
//       deriv(1,3)= -4.0*r;
//       deriv(1,4)= 4.0*r;
//       deriv(1,5)= 4.0*(1.0-2.0*s-r);
//     }
//   break;

///*------------------------------------------------- triangular elements */
//case tri3: /* LINEAR shape functions and their natural derivatives -----*/
///*----------------------------------------------------------------------*/
//   funct(0)=ONE-r-s;
//   funct[1]=r;
//   funct[2]=s;
//
//   if(option==1) /* --> first derivative evaluation */
//   {
//      deriv[0][0]=-ONE;
//      deriv[1][0]=-ONE;
//      deriv[0][1]= ONE;
//      deriv[1][1]=ZERO;
//      deriv[0][2]=ZERO;
//      deriv[1][2]= ONE;
//   } /* endif (option==1) */
//break;
///*-------------------------------------------------------------------------*/
//case tri6: /* Quadratic shape functions and their natural derivatives -----*/
//    t = ONE-r-s;
//
//   funct[0] = t*(TWO*t-ONE);
//    funct[1] = r*(TWO*r-ONE);
//    funct[2] = s*(TWO*s-ONE);
//    funct[3] = FOUR*r*t;
//    funct[4] = FOUR*r*s;
//    funct[5] = FOUR*s*t;
//
//    if (option == 1) /* --> first derivative evaluation */
//    {
//        /* first natural derivative of funct[0] with respect to r */
//        deriv[0][0] = -FOUR*t + ONE;
//        /* first natural derivative of funct[0] with respect to s */
//        deriv[1][0] = -FOUR*t + ONE;
//        deriv[0][1] = FOUR*r - ONE;
//        deriv[1][1] = ZERO;
//        deriv[0][2] = ZERO;
//        deriv[1][2] = FOUR*s - ONE;
//        deriv[0][3] = FOUR*t - FOUR*r;
//        deriv[1][3] = -FOUR*r;
//        deriv[0][4] = FOUR*s;
//        deriv[1][4] = FOUR*r;
//        deriv[0][5] = -FOUR*s;
//        deriv[1][5] = FOUR*t - FOUR*s;
//    } /* end if (option==1) */
//break;
default:
   dserror("Unknown no. of nodes %d to wall1 element",numnode);
break;
} /* end of switch typ */
/*----------------------------------------------------------------------*/

  return;

} /* DRT::Elements::Wall1::w1_shapefunctions */

/*----------------------------------------------------------------------*
 |  jacobian matrix (private)                                  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_jacobianmatrix(double xrefe[2][MAXNOD_WALL1],
                          const Epetra_SerialDenseMatrix& deriv,
                          Epetra_SerialDenseMatrix& xjm,
			  double* det,
                          const int iel)
{

   memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double));

   for (int k=0; k<iel; k++)
   {
        xjm(0,0) += deriv(0,k) * xrefe[0][k];
        xjm(0,1) += deriv(0,k) * xrefe[1][k];
        xjm(1,0) += deriv(1,k) * xrefe[0][k];
        xjm(1,1) += deriv(1,k) * xrefe[1][k];
   }

/*------------------------------------------ determinant of jacobian ---*/
     *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

      if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
/*----------------------------------------------------------------------*/

   return;
} // DRT::Elements::Wall1::w1_jacobianmatrix

/*----------------------------------------------------------------------*
 |  Matrix boplin in reference configuration (private)         mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_boplin(Epetra_SerialDenseMatrix& boplin,
                           Epetra_SerialDenseMatrix& deriv,
                           Epetra_SerialDenseMatrix& xjm,
                           double& det,
                           const int iel)
{

  int inode;
  int dnode;
  double dum;
  double xji[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  dum = 1.0/det;
  xji[0][0] = xjm(1,1)* dum;
  xji[0][1] =-xjm(0,1)* dum;
  xji[1][0] =-xjm(1,0)* dum;
  xji[1][1] = xjm(0,0)* dum;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,x    0   |
       |   0    Nk,y |
       | Nk,y    0   |
       |  0     Nk,x |
  */
  for (inode=0; inode<iel; inode++)
  {
    dnode = inode*2;

    boplin(0,dnode+0) = deriv(0,inode)*xji[0][0] + deriv(1,inode)*xji[0][1];
    boplin(1,dnode+1) = deriv(0,inode)*xji[1][0] + deriv(1,inode)*xji[1][1];
    boplin(2,dnode+0) = boplin(1,dnode+1);
    boplin(3,dnode+1) = boplin(0,dnode+0);
  } /* end of loop over nodes */
 return;
}

/* DRT::Elements::Wall1::w1_boplin */

/*----------------------------------------------------------------------*
 | Deformation gradient F and Green-Langrange strain (private)  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_defgrad(Epetra_SerialDenseVector& F,
                           Epetra_SerialDenseVector& strain,
                           const double xrefe[][MAXNOD_WALL1],
                           const double xcure[][MAXNOD_WALL1],
                           Epetra_SerialDenseMatrix& boplin,
                           const int iel)
{
  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,x  |
        |  1 + Uy,y  |
        |      Ux,y  |
        |      Uy,x  |
  */

  memset(F.A(),0,F.N()*F.M()*sizeof(double));

  F[0]=1;
  F[1]=1;
  for (int inode=0; inode<iel; inode++)
  {
     F[0] += boplin(0,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);
     F[1] += boplin(1,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);
     F[2] += boplin(2,2*inode)   * (xcure[0][inode] - xrefe[0][inode]);
     F[3] += boplin(3,2*inode+1) * (xcure[1][inode] - xrefe[1][inode]);
  } /* end of loop over nodes */
  /*-----------------------calculate Green-Lagrange strain ---------------*/
  strain[0]=0.5 * (F[0] * F[0] + F[3] * F[3] - 1.0);
  strain[1]=0.5 * (F[2] * F[2] + F[1] * F[1] - 1.0);
  strain[2]=0.5 * (F[0] * F[2] + F[3] * F[1]);
  strain[3]=strain[2];

 return;
}

/* DRT::Elements::Wall1::w1_defgrad */


/*----------------------------------------------------------------------*
 | Deformation gradient F in matrix notation and B in
 reference configuration (private)                             mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_boplin_cure(Epetra_SerialDenseMatrix& b_cure,
                                          Epetra_SerialDenseMatrix& boplin,
                                          Epetra_SerialDenseVector& F,
                                          int numeps,
                                          int nd)
{


     Epetra_SerialDenseMatrix Fmatrix;
     Fmatrix.Shape(4,4);


  /*---------------------------write Vector F as a matrix Fmatrix*/

     Fmatrix(0,0) = F[0];
     Fmatrix(0,2) = 0.5 * F[2];
     Fmatrix(0,3) = 0.5 * F[2];
     Fmatrix(1,1) = F[1];
     Fmatrix(1,2) = 0.5 * F[3];
     Fmatrix(1,3) = 0.5 * F[3];
     Fmatrix(2,1) = F[2];
     Fmatrix(2,2) = 0.5 * F[0];
     Fmatrix(2,3) = 0.5 * F[0];
     Fmatrix(3,0) = F[3];
     Fmatrix(3,2) = 0.5 * F[1];
     Fmatrix(3,3) = 0.5 * F[1];
    /*-------------------------------------------------int_b_cure operator*/
      memset(b_cure.A(),0,b_cure.N()*b_cure.M()*sizeof(double));
      for(int i=0; i<numeps; i++)
        for(int j=0; j<nd; j++)
          for(int k=0; k<numeps; k++)
          b_cure(i,j) += Fmatrix(k,i)*boplin(k,j);
    /*----------------------------------------------------------------*/

  return;
}

/* DRT::Elements::Wall1::w1_boplin_cure */

/*----------------------------------------------------------------------*
 | Constitutive matrix C and stresses (private)                mgit 05/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_call_matgeononl(Epetra_SerialDenseVector& strain,
                                              Epetra_SerialDenseMatrix& stress,
                                              Epetra_SerialDenseMatrix& C,
                                              const int numeps,
                                              struct _MATERIAL* material)

{
  /*--------------------------- call material law -> get tangent modulus--*/
    switch(material->mattyp)
    {
    case m_stvenant:/*--------------------------------- linear elastic ---*/
    {
      double ym = material->m.stvenant->youngs;
      double pv = material->m.stvenant->possionratio;


  /*-------------- some comments, so that even fluid people are able to
     understand this quickly :-)
     the "strain" vector looks like:

         | EPS_xx |
         | EPS_yy |
         | EPS_xy |
         | EPS_yx |

  */
  /*---------------------------------material-tangente-- plane stress ---*/
    switch(wtype_)
    {
    case plane_stress:
      {
      double e1=ym/(1. - pv*pv);
      double e2=pv*e1;
      double e3=e1*(1. - pv)/2.;

      C(0,0)=e1;
      C(0,1)=e2;
      C(0,2)=0.;
      C(0,3)=0.;

      C(1,0)=e2;
      C(1,1)=e1;
      C(1,2)=0.;
      C(1,3)=0.;

      C(2,0)=0.;
      C(2,1)=0.;
      C(2,2)=e3;
      C(2,3)=e3;

      C(3,0)=0.;
      C(3,1)=0.;
      C(3,2)=e3;
      C(3,3)=e3;
      }
      break;
      default:
     /*----------- material-tangente - plane strain, rotational symmetry ---*/
      {
      double c1=ym/(1.0+pv);
      double b1=c1*pv/(1.0-2.0*pv);
      double a1=b1+c1;

      C(0,0)=a1;
      C(0,1)=b1;
      C(0,2)=0.;
      C(0,3)=0.;

      C(1,0)=b1;
      C(1,1)=a1;
      C(1,2)=0.;
      C(1,3)=0.;

      C(2,0)=0.;
      C(2,1)=0.;
      C(2,2)=c1/2.;
      C(2,3)=c1/2.;

      C(3,0)=0.;
      C(3,1)=0.;
      C(3,2)=c1/2;
      C(3,3)=c1/2;
      }
     }
  /*-------------------------- evaluate 2.PK-stresses -------------------*/
  /*------------------ Summenschleife -> += (2.PK stored as vecor) ------*/

  Epetra_SerialDenseVector svector;
  svector.Size(4);

  for (int k=0; k<3; k++)
  {
        for (int i=0; i<numeps; i++)
    {
       svector[k] += C[k][i] * strain[i];
    }
  }
  /*------------------ 2.PK stored as matrix -----------------------------*/
  stress(0,0)=svector[0];
  stress(0,2)=svector[2];
  stress(1,1)=svector[1];
  stress(1,3)=svector[2];
  stress(2,0)=svector[2];
  stress(2,2)=svector[1];
  stress(3,1)=svector[2];
  stress(3,3)=svector[0];

  }

    break;
    default:
    dserror(" unknown type of material law");
    }

  return;
}

/* DRT::Elements::Wall1::w1_call_matgeononl */


/*----------------------------------------------------------------------*
| geometric stiffness part (total lagrange)                   mgit 05/07|
*----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_kg(Epetra_SerialDenseMatrix& estif,
                                 Epetra_SerialDenseMatrix& boplin,
                                 Epetra_SerialDenseMatrix& stress,
                                 double fac,
                                 int nd,
                                 int numeps)
{
  /*---------------------------------------------- perform B^T * SIGMA * B*/
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
      for(int r=0; r<numeps; r++)
         for(int m=0; m<numeps; m++)
            estif(i,j) += boplin(r,i)*stress(r,m)*boplin(m,j)*fac;

  return;
}

/* DRT::Elements::Wall1::w1_kg */

/*----------------------------------------------------------------------*
| elastic and initial displacement stiffness (total lagrange)  mgit 05/07                   mgit 05/07|
*----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_keu(Epetra_SerialDenseMatrix& estif,
                                  Epetra_SerialDenseMatrix& b_cure,
                                  Epetra_SerialDenseMatrix& C,
                                  double fac,
                                  int nd,
                                  int numeps)
{
  /*------------- perform B_cure^T * D * B_cure, whereas B_cure = F^T * B */
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
        for(int k=0; k<numeps; k++)
           for(int m=0; m<numeps; m++)

            estif(i,j) +=  b_cure(k,i)*C(k,m)*b_cure(m,j)*fac;

  return;
}

/* DRT::Elements::Wall1::w1_keu */

/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) mgit 05/07  |
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_fint(Epetra_SerialDenseMatrix& stress,
                                   Epetra_SerialDenseMatrix& b_cure,
                                   Epetra_SerialDenseVector& intforce,
                                   double fac,
                                   int nd)

{
  Epetra_SerialDenseVector st;
  st.Size(4);

  st[0] = fac * stress(0,0);
  st[1] = fac * stress(1,1);
  st[2] = fac * stress(0,2);
  st[3] = fac * stress(0,2);

  for(int i=0; i<nd; i++)
    for(int j=0; j<4; j++)
      intforce[i] += b_cure(j,i)*st[j];

  return;
}

/* DRT::Elements::Wall1::w1_fint */






#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
