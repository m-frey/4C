/*!----------------------------------------------------------------------
\file
\brief setting dirichlet conditions for fluid

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;

/*!--------------------------------------------------------------------- 
\brief routine to initialise the dirichlet boundary conditions

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are initialised:
- real pressure from input-file is transformed to kinematic pressure:
  actgnode->dirich->dirich_val.a.dv[predof] /= dens

- if the initial field is a zero-field, then set the dirichlet values in
  the solution history of the node:
    'actnode->sol.a.da[0][j] = initval*acttimefac' (--> output)
    'actnode->sol_increment.a.da[1][j] = initval*acttimefac'
                                 |        |         |
                            time (n)      |         |               
                          initial value from input  |               
                                       factor from timecurve (T=0.0)

</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)    

\return void                                            

------------------------------------------------------------------------*/
void fluid_initdirich(  FIELD          *actfield, 
                        FLUID_DYNAMIC  *fdyn
		     )
{
INT        i,j;
INT        numnp_total;               /* total number of fluid nodes    */
INT        numele_total;              /* total number of fluid elements */
INT        predof;	              /* number of pressure dof	        */
INT        numdf;	              /* number of fluid dofs	        */
INT        actcurve;	              /* actual timecurve  	        */
INT        numveldof;
DOUBLE     dens;	              /* density			*/
DOUBLE     timefac[MAXTIMECURVE];     /* factors from time-curve        */
DOUBLE     T=0.0;	              /* starting time		        */
DOUBLE     acttimefac;                /* actual factor from timecurve   */
DOUBLE     initval;	              /* intial dirichlet value	        */
GNODE     *actgnode;	              /* actual GNODE		        */
NODE      *actnode;	              /* actual NODE		        */
ELEMENT   *actele;	              /* actual ELEMENT		        */

INT counter=0;

#ifdef DEBUG 
dstrc_enter("fluid_initdirich");
#endif  

numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele; 
numdf        = fdyn->numdf; 
predof       = numdf-1;
numveldof    = numdf-1;

/*-------------------------- since different materials are not allowed
              one can work with the material parameters of any element */
actele = &(actfield->dis[0].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

/*------------------------------------------ check dirichlet conditions */
for (i=0;i<numnp_total;i++) /* loop all nodes */
{
   actnode  = &(actfield->dis[0].node[i]); 
   actgnode = actnode->gnode; 
   if (actgnode->dirich==NULL)
      continue;
   if (actgnode->dirich->dirich_type==dirich_FSI)
      counter++;
   if (actgnode->dirich->dirich_type==dirich_none)
   { 
      for (j=0;j<actnode->numdf;j++) /* loop all dofs */    
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
            actcurve = actgnode->dirich->curve.a.iv[j];
	    if(actcurve>numcurve)
	       dserror("Load curve: actual curve > number defined curves\n");   
      } /* end of loop over all dofs */
      /* transform real pressure from input to kinematic pressure ---*/
      if (actgnode->dirich->dirich_onoff.a.iv[predof]!=0)      
          actgnode->dirich->dirich_val.a.dv[predof] /= dens; 
   }
} /* end of loop over all nodes */

if (counter>0 && par.myrank==0)
{
printf("\n");
printf("          | FIELD FLUID     | number of nodes coupled with structure: %d \n",counter);
printf("\n");
}

/*---------- set dirichlet conditions at time (0) for zero intial field */
if (fdyn->init==0)
{
/*------------------------------------------ get values from time curve */
   for (actcurve=0;actcurve<numcurve;actcurve++)
   {
     dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
   }/* end loop over active timecurves */   
/*------------------------------------------------- loop over all nodes */
   for (i=0;i<numnp_total;i++)
   {
      actnode  = &(actfield->dis[0].node[i]); 
      actgnode = actnode->gnode;      
      if (actgnode->dirich==NULL)
         continue;
      switch(actgnode->dirich->dirich_type)
      {
      case dirich_none:
         for (j=0;j<actnode->numdf;j++) /* loop all dofs */
         {
            if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
               continue;
            actcurve = actgnode->dirich->curve.a.iv[j]-1;
            if (actcurve<0)
               acttimefac = ONE;
            else
               acttimefac = timefac[actcurve];
            initval  = actgnode->dirich->dirich_val.a.dv[j];               
            actnode->sol_increment.a.da[1][j] = initval*acttimefac;
	    actnode->sol.a.da[0][j] = initval*acttimefac;
         } /* end loop over dofs */
      break;
      case dirich_FSI: /* FSI --> dirichvalues = grid velocity!!! */
	 for (j=0;j<numveldof;j++) /* loop vel-dofs */
	 {
	    initval = actnode->sol_increment.a.da[4][j];  
	    actnode->sol_increment.a.da[1][j] = initval; 
	    actnode->sol.a.da[0][j] = initval;
	 }
      break;
      default:
         dserror("dirch_type unknown!\n");
      } /* end switch */
   } /*end loop over nodes */   
} /* endif fdyn->init */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_initdirich*/

/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions at time <time>

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->time>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[pos][j] = initval*acttimefac'
                                  |         |         |
                            time (n+1)      |         |               
                          initial value from input    |               
                                         factor from timecurve
</pre>
\param *actfield FIELD		(i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC	(i)  
\param  pos	 INT		(i)	position, where to write dbc

\return void     

------------------------------------------------------------------------*/
void fluid_setdirich(   FIELD           *actfield, 
                        FLUID_DYNAMIC   *fdyn,
			INT		 pos
	            )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        actcurve;	             /* actual timecurve		*/
INT        numveldof;
DOUBLE     timefac[MAXTIMECURVE];    /* factors from time-curve         */
DOUBLE     T;		             /* actual time		        */
DOUBLE     acttimefac;               /* actual factor from timecurve    */
DOUBLE     initval;	             /* intial dirichlet value	        */
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG 
dstrc_enter("fluid_setdirich");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
T            = fdyn->time;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
} /* end loop over active timecurves */ 

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++) 
{
   actnode  = &(actfield->dis[0].node[i]); 
   actgnode = actnode->gnode;      
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
         actcurve = actgnode->dirich->curve.a.iv[j]-1;
         if (actcurve<0)
            acttimefac = ONE;
         else
            acttimefac = timefac[actcurve];
         initval  = actgnode->dirich->dirich_val.a.dv[j];               
         actnode->sol_increment.a.da[pos][j] = initval*acttimefac;	 
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */     
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[pos][j]
	=actnode->sol_increment.a.da[4][j];
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich*/


/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions for fluid to determine
Relaxation parameter

<pre>                                                             cf 08/03

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set for the fluid solution needed to determine the Relaxation
parameter via steepest descent method.

the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[7][j] = 0.0 
                                 |          | 
                fluid sol. for RelaxParam   | 
					at Dirichlet boundaries
  AND:
    'actnode->sol_increment.a.da[7][j] = actnode->sol_increment.a.da[4][j]
                                 |                  | 
                fluid sol. for RelaxParam           | 
					at fsi coupling interface,
					      grid velocity
</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)  

\return void     

------------------------------------------------------------------------*/
void fluid_setdirich_sd(
                        FIELD           *actfield, 
                        FLUID_DYNAMIC   *fdyn
	               )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        numveldof;
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG 
dstrc_enter("fluid_setdirich_sd");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++) 
{
   actnode  = &(actfield->dis[0].node[i]); 
   actgnode = actnode->gnode;      
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
         actnode->sol_increment.a.da[7][j] = 0.0;
         actnode->sol_increment.a.da[6][j] = 0.0;	 	
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */     
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
         actnode->sol_increment.a.da[7][j]
	=actnode->sol_increment.a.da[4][j];
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich_sd*/




/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions of acceleration

<pre>                                                             cf 09/03

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set for the fluid acceleration. These conditions depend 
upon the dbc for the velocity set in the input file according to the 
time stepping sceme. 
The actual implementation serves the generalised-alpha scheme.

sol_increment[pos_to][i] = fac1 * sol_increment[pos_from1][i]
                         + fac2 * sol_increment[pos_from2][i]
                         + fac3 * sol_increment[pos_from3][i]


</pre>
\param *actfield	FIELD		(i)	actual field (fluid)   
\param *fdyn	 	FLUID_DYNAMIC	(i)	fluid dynamic
\param  pos_to		INT		(i)	pos in sol_increment to write to
\param  pos1_from	INT		(i)	1st pos to read from
\param  pos2_from	INT		(i)	2nd pos to read from
\param  pos3_from	INT		(i)	3rd pos to read from
\param	fac1		DOUBLE		(i)	factor of value at pos1_from
\param	fac2		DOUBLE		(i)	factor of value at pos2_from
\param	fac3		DOUBLE		(i)	factor of value at pos3_from

\return void     

------------------------------------------------------------------------*/
void fluid_setdirich_acc(
                         FIELD		*actfield, 
                         FLUID_DYNAMIC	*fdyn,
			 INT		 pos_to,
			 INT		 pos1_from,
			 INT		 pos2_from,
			 INT		 pos3_from,
			 DOUBLE		 fac1,
			 DOUBLE		 fac2,
			 DOUBLE		 fac3
	                )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	             /* number of fluid dofs    	*/
INT        numveldof;
GNODE     *actgnode;	             /* actual GNODE		        */
NODE      *actnode;                  /* actual NODE                     */

#ifdef DEBUG 
dstrc_enter("fluid_setdirich_acc");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numdf        = fdyn->numdf;
numveldof    = numdf-1;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++) 
{
   actnode  = &(actfield->dis[0].node[i]); 
   actgnode = actnode->gnode;      
   if (actgnode->dirich==NULL)
         continue;
   switch(actgnode->dirich->dirich_type)
   {
   case dirich_none:
      for (j=0;j<actnode->numdf;j++) /* loop dofs */
      {
         if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
	       actnode->sol_increment.a.da[pos_to][j] 
	       = fac1 * actnode->sol_increment.a.da[pos1_from][j]
               + fac2 * actnode->sol_increment.a.da[pos2_from][j]
               + fac3 * actnode->sol_increment.a.da[pos3_from][j];
      } /* end loop over dofs */
   break;
   case dirich_FSI: /* dirichvalues = grid velocity!!! */     
      for (j=0;j<numveldof;j++)  /* loop vel-dofs */
	 dserror("generalised alpha with FSI not yet implemented");
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich_acc*/



/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector

<pre>                                                         genk 04/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[3][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= estif[i][j] * dirich[j];
\endcode			     

\param  *actele    ELEMENT   (i)   actual element	  
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces
\param   is_relax  INT       (i)   flag if it is for relaxation param

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_caldirich(
                        ELEMENT         *actele,  
		        DOUBLE          *dforces, 
                        DOUBLE         **estif,   
		        INT             *hasdirich,
			INT              is_relax
		    )     
{

INT         i,j;
INT         nrow;
INT         numdf;                      /* number of fluid dofs         */
INT         nd=0;                      
DOUBLE      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele  */ 
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */

#ifdef DEBUG 
dstrc_enter("fluid_caldirich");
#endif  

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;   
   if (actgnode->dirich==NULL) 
      continue;
   else
      *hasdirich=1;
      break;
} /* end loop over nodes */					  

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                           written to the nodes (sol_increment[3][j]) */
nrow=0;
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];   
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[nrow+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      if (is_relax) /* calculation for relax.-Param. reads
                       dbc from (sol_increment[7][j]) */
         dirich[nrow+j] = actnode->sol_increment.a.da[7][j];
      else 
         dirich[nrow+j] = actnode->sol_increment.a.da[3][j];
   } /* end loop over dofs */
   nrow+=numdf;
} /* end loop over nodes */
dsassert(nrow==nd,"failure during calculation of dirich forces\n");
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_caldirich*/ 

#endif
/*! @} (documentation module close)*/
