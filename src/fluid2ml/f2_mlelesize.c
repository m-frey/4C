/*!----------------------------------------------------------------------
\file
\brief Calculate submesh element size and streamlength for fluid2

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief routine to calculate submesh element size for fluid2

<pre>                                                       gravem 07/03

In this routine, the characteristic submesh element length is calculated
and the routine for the calculation of the stabilization parameter or 
the subgrid viscosity (depending on the respective flag) is called. 

</pre>
\param  *ele       ELEMENT	     (i)   actual element
\param  *data      FLUID_DATA	     (i)
\param  *dynvar    FLUID_DYN_CALC    (i/o)
\param  *mlvar     FLUID_DYN_ML      (i)
\param  *funct     DOUBLE 	     (-)   l-s shape functions
\param **deriv     DOUBLE 	     (-)   deriv. of l-s shape funcs
\param **deriv2    DOUBLE 	     (-)   2nd deriv. of l-s sh. funcs
\param  *smfunct   DOUBLE 	     (-)   submesh shape functions
\param **smderiv   DOUBLE 	     (-)   deriv. of submesh shape fun.
\param **smderiv2  DOUBLE 	     (-)   2nd deriv. of sm sh. fun.
\param **derxy     DOUBLE 	     (-)   global l-s sh. fun. deriv.
\param **xjm       DOUBLE 	     (-)   jacobian matrix
\param **evel      DOUBLE 	     (i)   l-s element velocities
\param  *velint    DOUBLE 	     (-)   l-s vel. at integr. point
\param **vderxy    DOUBLE 	     (-)   l-s vel. deriv. at int. point
\param **smxyze    DOUBLE 	     (i)   submesh element coordinates
\param **smxyzep   DOUBLE 	     (i)   sm ele. coord. on parent dom.
\param **cutp      DOUBLE 	     (-)   cutting points
\return void             

------------------------------------------------------------------------*/
void f2_smelesize(ELEMENT         *ele,    
		  FLUID_DATA	  *data, 
		  FLUID_DYN_CALC  *dynvar,
		  FLUID_DYN_ML    *mlvar,
	          DOUBLE	  *funct,  
	          DOUBLE	 **deriv,  
	          DOUBLE	 **deriv2,		
	          DOUBLE	  *smfunct,  
	          DOUBLE	 **smderiv,  
	          DOUBLE	 **smderiv2,		
	          DOUBLE	 **derxy,  
		  DOUBLE	 **xjm,    
		  DOUBLE	 **evel,		 
		  DOUBLE	  *velint, 
	          DOUBLE	 **vderxy,  
	          DOUBLE	 **smxyze,  
	          DOUBLE	 **smxyzep,  
		  DOUBLE	 **cutp)
{

INT     i;              /* simply a counter	        		*/
INT     ntyp,nsmtyp;    /* l-s and sm element type (TRI or QUAD)        */
INT     actmat;         /* number of actual material		        */
INT     iel,smiel;      /* l-s and sm number of nodes of actual element */
DOUBLE  visc;           /* fluid viscosity                              */
DOUBLE  area;           /* element area                                 */
DOUBLE  det;            /* determinant of jacobian                      */
DOUBLE  strle;          /* streamlength                                 */
DOUBLE  e1,e2;          /* natural coordinates of inegration point      */
DOUBLE  fac,facr,facs;  /* factors                                      */
DOUBLE  dia,dia1,dia2;  /* values used for calculation of element size  */
DOUBLE  dx,dy;          /* values used for calculation of element size  */
DOUBLE  gcoor[2],coor[2];/* global and local coordinates                */
DIS_TYP typ,smtyp;

#ifdef DEBUG 
dstrc_enter("f2_smelesize");
#endif		

/*---------------------------------------------------------- initialize */
ntyp   = ele->e.f2->ntyp;
nsmtyp = mlvar->submesh.ntyp;
typ    = ele->distyp;
smtyp  = mlvar->submesh.typ;
iel    = ele->numnp;
smiel  = mlvar->submesh.numen;
area  = ZERO;
strle = ZERO;
 
actmat=ele->mat-1;
visc = mat[actmat].m.fluid->viscosity;

/*------ get values of integration parameters, shape functions and their
         derivatives for submesh element -------------------------------*/
switch(nsmtyp)
{
case 1:    /* --> quad - element */
   e1	= data->qxg[0][0];
   facr = data->qwgt[0][0];
   e2	= data->qxg[0][0];
   facs = data->qwgt[0][0];
   f2_rec(smfunct,smderiv,smderiv2,e1,e2,smtyp,2);
break;
case 2:       /* --> tri - element */		   
   e1	= data->txgr[0][0];
   facr = data->twgt[0][0];
   e2	= data->txgs[0][0];
   facs = ONE;
   f2_tri(smfunct,smderiv,smderiv2,e1,e2,smtyp,2);
break;
default:
   dserror("nsmtyp unknown!\n");      
} /*end switch(nsmtyp) */

if (mlvar->smesize<4)     /* compute area */
{
/* -------------------------------------------- compute jacobian matrix */      
  f2_jaco3(smxyze,smfunct,smderiv,xjm,&det,smiel,ele);
  fac=facr*facs*det;
  area += fac;
}/* endif (mlvar->smesize<4) */
  
if (mlvar->smesize==4)    /* compute diagonal based diameter */
{
  switch(nsmtyp)
  {
  case 1:
    dx = smxyze[0][0] - smxyze[0][2];
    dy = smxyze[1][0] - smxyze[1][2];
    dia1 = sqrt(dx*dx+dy*dy);
    dx = smxyze[0][1] - smxyze[0][3];
    dy = smxyze[1][1] - smxyze[1][3];
    dia2 = sqrt(dx*dx+dy*dy);
/*------ dia=sqrt(2)*area/(1/2*(dia1+dia2))=sqrt(8)*area/(dia1+dia2) ---*/
    dia = sqrt(EIGHT)*area/(dia1+dia2); 
  break;
  case 2:    /* get global coordinate of element center */
    f2_gcoor2(smfunct,smxyze,smiel,gcoor);
    dia = ZERO;
    for (i=0;i<3;i++)
    {
       dx = gcoor[0] - smxyze[0][i];
       dy = gcoor[1] - smxyze[1][i];
       dia += dx*dx + dy*dy;
    }
    dia = FOUR*area/sqrt(THREE*dia);
  break;
  default:
    dserror("ntyp unknown!\n");
  } /* end switch(nsmtyp) */
} /* endif (mlvar->smesize==4) */

if (mlvar->smesize==5)    /* compute streamlength based on l-s velocity */
{
  f2_gcoor2(smfunct,smxyzep,smiel,coor);
  switch(ntyp)
  {
  case 1:    /* --> quad - element */
    f2_rec(funct,deriv,deriv2,coor[0],coor[1],typ,2);
  break;
  case 2:	/* --> tri - element */ 	     
    f2_tri(funct,deriv,deriv2,coor[0],coor[1],typ,2);
  break;
  default:
    dserror("ntyp unknown!\n");      
  } /*end switch(ntyp) */
  f2_veli(velint,funct,evel,iel);
  f2_gcoor2(smfunct,smxyze,smiel,gcoor);
  f2_smstrlen(&strle,velint,smxyze,gcoor,cutp,ntyp);	      
} /* endif (mlvar->smesize==5) */

/*----------------------------------- set characteristic element length */
if      (mlvar->smesize==1) ele->e.f2->smcml = sqrt(area);
else if (mlvar->smesize==2) ele->e.f2->smcml = TWO*sqrt(area/PI);
else if (mlvar->smesize==3) ele->e.f2->smcml = sqrt(TWO*area/PI);
else if (mlvar->smesize==4) ele->e.f2->smcml = dia;
else if (mlvar->smesize==5) ele->e.f2->smcml = strle;  

if (mlvar->smesize<5)   /* compute l-s velocity */
{
  f2_gcoor2(smfunct,smxyzep,smiel,coor);
  switch(ntyp)
  {
  case 1:    /* --> quad - element */
    f2_rec(funct,deriv,deriv2,coor[0],coor[1],typ,2);
  break;
  case 2:	/* --> tri - element */ 	     
    f2_tri(funct,deriv,deriv2,coor[0],coor[1],typ,2);
  break;
  default:
    dserror("ntyp unknown!\n");      
  } /*end switch(ntyp) */
  f2_veli(velint,funct,evel,iel);
}  

/*----------------------------------- calculate stabilization parameter */               
if (mlvar->smstabi>0) f2_smstabpar(ele,dynvar,mlvar,velint,visc,smiel,ntyp); 

/*--------------------------------------------------- subgrid viscosity */               
if (mlvar->smsgvi==1 || mlvar->smsgvi==2)
{ 
  f2_jaco(funct,deriv,xjm,&det,ele,iel);
/*-------------------------------------------- compute global derivates */
  f2_gder(derxy,deriv,xjm,det,iel);
/*----------------------- get velocity derivatives at integration point */
  f2_vder(vderxy,derxy,evel,iel);
/*----------------------------------------- calculate subgrid viscosity */               
  f2_smsgvisc(ele,mlvar,velint,vderxy,visc,smiel,ntyp);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_smelesize */

/*!---------------------------------------------------------------------                                         
\brief routine to calculate streamlength for submesh for fluid2

<pre>                                                       gravem 07/03

In this routine, the submesh element streamlength for calculation of 
stabilization parameter is calculated.
		     
</pre>
\param  *strle     DOUBLE   (o)    streamlength
\param  *velint    DOUBLE   (i)    l-s velocities at integr. point
\param **smxyze    DOUBLE   (i)    submesh element coordinates
\param  *gcoor     DOUBLE   (i)    global coord. of int. point
\param **cutp      DOUBLE   (-)    cutting points
\param   ntyp	   INT      (i)    flag for element type
\return void                                               
\sa f2_calelesize()                               

------------------------------------------------------------------------*/
void f2_smstrlen(DOUBLE   *strle,     
		 DOUBLE   *velint,   
		 DOUBLE  **smxyze,	
                 DOUBLE   *gcoor,    
		 DOUBLE  **cutp,	     
		 int	   ntyp)
{
INT     nodcut=-1;
INT     nodmax;
INT     inod;
DOUBLE dl,dx,dy,dxh,dyh;
DOUBLE x1,y1,x2,y2;
DOUBLE dsub,dval;

#ifdef DEBUG 
dstrc_enter("f2_smstrlen");
#endif

dval = FABS(velint[0])+FABS(velint[1]);
if (dval == ZERO)  /* no flow at this point - take some arbitr. measure for streamlength */
{
   dx = smxyze[0][2] - smxyze[0][0];
   dy = smxyze[1][2] - smxyze[1][0];
   goto calc2;   
} /* enidf (dval == ZERO) */

/*----------------------------------------------------------------------*
   streamlength is calculated via cutting points of velocity vector
   with straight boundaries                                             
*/
switch(ntyp)
{
case 1: /* max number of nodes for quad: 4 --> C-numbering nodmax = 4-1 */
   nodmax = 3;
break;
case 2:  /* max number of nodes for tri: 3 --> C-numbering nodmax = 3-1 */
   nodmax = 2;
break;
default:
   dserror("ntyp unknown!\n");   
} /* end switch(ntyp) */        
 /*------------------------------------------------- get cutting points */
for (inod=0;inod<nodmax;inod++)
{
   dxh = smxyze[0][inod+1] - smxyze[0][inod];
   dyh = smxyze[1][inod+1] - smxyze[1][inod];
   dsub = dxh*velint[1]-dyh*velint[0];
   if (dsub==ZERO)  /* check for parallel vectors */
      continue;
   dl = ((smxyze[1][inod]-gcoor[1])*velint[0] -	\
	 (smxyze[0][inod]-gcoor[0])*velint[1])/dsub;
   if (dl>=ZERO && dl<=ONE)
   {
      nodcut++;
      cutp[0][nodcut]=smxyze[0][inod]+dl*dxh;
      cutp[1][nodcut]=smxyze[1][inod]+dl*dyh;
      if (nodcut==1)
	 goto calc1;
   } /* endif (dl>=ZERO && dl<=ONE) */
} /* end loop over inod */
/*------------------------------------------------- the last boundary */
dxh = smxyze[0][0]-smxyze[0][nodmax];
dyh = smxyze[1][0]-smxyze[1][nodmax];

dsub = dxh*velint[1] - dyh*velint[0];
if (dsub==ZERO)
   dserror("Couldn't find two cutting points!\n");
dl = ((smxyze[1][nodmax]-gcoor[1])*velint[0] -	\
      (smxyze[0][nodmax]-gcoor[0])*velint[1])/dsub;
if (dl>=ZERO && dl <= ONE)
{
   nodcut++;
   cutp[0][nodcut]=smxyze[0][nodmax]+dl*dxh;
   cutp[1][nodcut]=smxyze[1][nodmax]+dl*dyh; 
   if(nodcut==1)
      goto calc1;
} /* endif  (dl>=ZERO && dl <= ONE) */

dserror("Couldn't find two cutting points!\n");

calc1:
dx = cutp[0][1]-cutp[0][0];
dy = cutp[1][1]-cutp[1][0];

calc2:
*strle = sqrt(dx*dx+dy*dy);

/*----------------------------------------------------------------------*/		
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_smstrlen */		  
#endif
