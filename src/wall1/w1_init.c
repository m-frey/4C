#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | initialize the element                                    al 6/01    |
 *----------------------------------------------------------------------*/
void w1init(PARTITION *actpart,MATERIAL    *mat )
{
int          i,j,k,ncm;
int          size_i, size_j;
ELEMENT     *actele;
NODE        *actnode;
W1_DATA      data;

ARRAY    funct_a_h;  /* shape functions */    
double  *funct_h;     
ARRAY    deriv_a_h;  /* derivatives of shape functions */   
double **deriv_h;     
ARRAY    xjm_a_h;    /* jacobian matrix */     
double **xjm_h;         

#ifdef DEBUG 
dstrc_enter("w1init");
#endif
/*----------------------------------------------------------------------*/
funct_h     = amdef("funct_h"  ,&funct_a_h,MAXNOD_WALL1,1 ,"DV");       
deriv_h     = amdef("deriv_h"  ,&deriv_a_h,2,MAXNOD_WALL1 ,"DA");       
xjm_h       = amdef("xjm_h"    ,&xjm_a_h  ,2,2            ,"DA");           
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numele; i++)
{
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_wall1) continue;
  /*---------------------------------------- init integration points ---*/
  w1intg(actele,&data,0);
  
  /*----------------------------------------- init stress structures ---*/
  size_i = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
  size_j = actele->numnp;
  
  actele->e.w1->stress = (W1_ELE_STRESS*)calloc(1,sizeof(W1_ELE_STRESS));
  
  actele->e.w1->stress[0].gprr = (double*)calloc(size_i,sizeof(double));
  actele->e.w1->stress[0].gpss = (double*)calloc(size_i,sizeof(double));
  actele->e.w1->stress[0].gprs = (double*)calloc(size_i,sizeof(double));
  
  actele->e.w1->stress[0].fps  = (double*)calloc(size_i,sizeof(double));
  actele->e.w1->stress[0].sps  = (double*)calloc(size_i,sizeof(double));
  actele->e.w1->stress[0].aps  = (double*)calloc(size_i,sizeof(double));
  
  actele->e.w1->stress[0].ferr = (double*)calloc(size_j,sizeof(double));
  actele->e.w1->stress[0].fess = (double*)calloc(size_j,sizeof(double));
  actele->e.w1->stress[0].fers = (double*)calloc(size_j,sizeof(double));
  /*--------------------------------------------- init working array ---*/
  if(mat[actele->mat-1].mattyp == m_pl_mises || 
     mat[actele->mat-1].mattyp == m_pl_dp || 
     mat[actele->mat-1].mattyp == m_pl_epc )
  {
    size_i = 1;
    actele->e.w1->elewa = (W1_ELE_WA*)CALLOC(size_i,sizeof(W1_ELE_WA));
    if (actele->e.w1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    } 
    
    size_j = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
    actele->e.w1->elewa[0].ipwa = 
                               (W1_IP_WA*)CALLOC(size_j,sizeof(W1_IP_WA));
    if (actele->e.w1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    } 
    for (k=0; k<size_j; k++)
    {
      actele->e.w1->elewa[0].ipwa[k].epstn = 0.;
      actele->e.w1->elewa[0].ipwa[k].yip   = -1;
      actele->e.w1->elewa[0].ipwa[k].qn = (double*)calloc(4,sizeof(double));

      if(mat[actele->mat-1].mattyp == m_pl_epc )
      {
      actele->e.w1->elewa[0].ipwa[k].sigc = (double*)calloc(4,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].grad = (double*)calloc(4,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].dlam = (double*)calloc(2,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].dlam[0] = 0.;
      actele->e.w1->elewa[0].ipwa[k].dlam[1] = 0.;
      actele->e.w1->elewa[0].ipwa[k].sigi = (double*)calloc(4,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].epsi = (double*)calloc(4,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].di   = (double*)calloc(4,sizeof(double));
      }
      
      ncm = mat[actele->mat-1].m.pl_epc->maxreb;
      
      actele->e.w1->elewa[0].ipwa[k].rsig   = (double*)calloc(ncm,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].reps   = (double*)calloc(ncm,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].repstn = (double*)calloc(ncm,sizeof(double));
      actele->e.w1->elewa[0].ipwa[k].ryip   = (int*)calloc(ncm,sizeof(int));
      for (j=0; j<ncm; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].rsig[j]   = 0.;
        actele->e.w1->elewa[0].ipwa[k].reps[j]   = 0.;
        actele->e.w1->elewa[0].ipwa[k].repstn[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].ryip[j]   = -1;
      }
      
      
      for (j=0; j<4; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].sig[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].eps[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].qn[ j] = 0.;
        if(mat[actele->mat-1].mattyp == m_pl_epc )
        {
        actele->e.w1->elewa[0].ipwa[k].sigc[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].grad[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].sigi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].epsi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].di[  j] = 0.;
        }
      }
    }
  /*------------------------------------- calculate element diameter ---*/
    if(mat[actele->mat-1].mattyp == m_pl_mises && 
     (fabs(0.0001 - mat[actele->mat-1].m.pl_mises->GF) > 0.0001) )
    {
       w1cdia(actele, &data, funct_h, deriv_h, xjm_h);
    }
    else if(mat[actele->mat-1].mattyp == m_pl_epc)
    {
       w1cdia(actele, &data, funct_h, deriv_h, xjm_h);
    }
   /*-------------------------------------------------------------------*/
  }
}
/*----------------------------------------------------------------------*/
amdel(&funct_a_h);       
amdel(&deriv_a_h);       
amdel(&xjm_a_h  );           
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1init */
 
