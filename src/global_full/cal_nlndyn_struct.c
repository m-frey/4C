#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];


double acttime;

/*----------------------------------------------------------------------*
 |  routine to control nonlinear dynamic structural analysis m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_nln_structural() 
{
int             i,j,k;                  /* simply a counter */
int             numeq;              /* number of equations on this proc */
int             numeq_total;        /* total number of equations */
int             init;               /* flag for solver_control call */
int             itnum;              /* counter for NR-Iterations */
int             convergence;        /* convergence flag */
int             mod_disp,mod_stress;
int             mod_res_write;
int             restart;
double          maxtime;
double          t0_res,t1_res;

double          dt;
int             nstep;

double          t0,t1;

double          dmax;               /* infinity norm of residual displacements */

int             stiff_array;        /* indice of the active system sparse matrix */
int             mass_array;         /* indice of the active system sparse matrix */
int             damp_array;         /* indice of the active system sparse matrix */
int             num_array;          /* indice of global stiffness matrices */
int             actcurve;           /* indice of active time curve */

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structure cal_action enum */
STRUCT_DYNAMIC *sdyn;               /* pointer to structural dynamic input data */

DIST_VECTOR    *vel;                /* total velocities */              
DIST_VECTOR    *acc;                /* total accelerations */
DIST_VECTOR    *fie;                /* internal forces and working array */
DIST_VECTOR    *dispi;              /* distributed vector to hold incremental displacments */ 
DIST_VECTOR    *work;               /* working vectors */

ARRAY           intforce_a;         /* redundant vector of full length for internal forces */
double         *intforce;
ARRAY           dirich_a;           /* redundant vector of full length for dirichlet-part of rhs */
double         *dirich;             
double          dirichfacs[10];     /* factors needed for dirichlet-part of rhs */
 
STRUCT_DYN_CALC dynvar;             /* variables to perform dynamic structural simulation */              

CONTAINER       container;          /* contains variables defined in container.h */
container.isdyn = 1;                /* dynamic calculation */

#ifdef DEBUG 
dstrc_enter("dyn_nln_structural");
#endif
/*----------------------------------------------------------------------*/
restart = genprob.restart;
/*--------------------------------------------------- set some pointers */
actfield           = &(field[0]);
actsolv            = &(solv[0]);
actpart            = &(partition[0]);
action             = &(calc_action[0]);
sdyn               =   alldyn[0].sdyn;
container.fieldtyp = actfield->fieldtyp;
/*----------------------------------------- check for explicit dynamics */
if (sdyn->Typ == centr_diff)
{
   dyn_nln_stru_expl();
   goto end_expl;
}
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
/* if we are not parallel, we have to allocate an alibi intra-communicator structure */
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*-------------------------------- init the variables in dynvar to zero */
dynvar.rldfac = 0.0;
dynvar.rnorm  = 0.0;
dynvar.epot   = 0.0;
dynvar.eout   = 0.0;
dynvar.etot   = 0.0;
dynvar.ekin   = 0.0;
dynvar.dinorm = 0.0;
dynvar.dnorm  = 0.0;
for (i=0; i<20; i++) dynvar.constants[i] = 0.0;
acttime=0.0;
/*------------------------------------ check presence of damping matrix */
/*                 and set indice of stiffness and mass sparse matrices */
   stiff_array = 0;
   mass_array  = 1;

if (sdyn->damp==1) 
{
   damp_array  = 2;
   actsolv->nsysarray=3;
}
else
{
   damp_array  =-1;
   actsolv->nsysarray=2;
}

/*--------------- stiff_array already exists, so copy the mask of it to */
/*------------------------------- mass_array (and damp_array if needed) */
/* reallocate the vector of sparse matrices and the vector of there types */
/* formerly lenght 1, now lenght 2 or 3 dependent on presence of damp_array */
actsolv->sysarray_typ = 
(SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

actsolv->sysarray = 
(SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

/*-copy the matrices sparsity mask from stiff_array to mass_array (and to damp_array) */
solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[stiff_array]),
                            &(actsolv->sysarray[stiff_array]),
                            &(actsolv->sysarray_typ[mass_array]),
                            &(actsolv->sysarray[mass_array]));
if (damp_array>0)
{
solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[stiff_array]),
                            &(actsolv->sysarray[stiff_array]),
                            &(actsolv->sysarray_typ[damp_array]),
                            &(actsolv->sysarray[damp_array]));
}

/*------------------------------- init the dist sparse matrices to zero */
for (i=0; i<actsolv->nsysarray; i++)
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[i]),
                 &(actsolv->sysarray_typ[i])
                );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(actsolv->sysarray[stiff_array],
                   actsolv->sysarray_typ[stiff_array],
                   &numeq,
                   &numeq_total
                   );


/*---------------------------------------allocate 4 dist. vectors 'rhs' */
/*  these hold original load vector, load vector at time t and t-dt and */
/*                                             interpolated load vector */
actsolv->nrhs = 4;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));

/*-------------------- there are 2 solution vector to hold total displ.*/
/*                                  one at time t and one at time t-dt */
actsolv->nsol= 2;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));

/*-------------- there is one vector to hold incremental displacements */
solserv_create_vec(&dispi,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(dispi[i]));

/*-------------------------------------------- allocate one vector vel */
solserv_create_vec(&vel,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(vel[i]));

/*-------------------------------------------- allocate one vector acc */
solserv_create_vec(&acc,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(acc[i]));

/*-------------- allocate one redundant vector intforce of full lenght */
/* this is used by the element routines to assemble the internal forces*/
intforce = amdef("intforce",&intforce_a,numeq_total,1,"DV");
/*----------- create a vector of full length for dirichlet part of rhs */
dirich = amdef("dirich",&dirich_a,numeq_total,1,"DV");
/*----------------------------------------- allocate 3 DIST_VECTOR fie */
/*                    to hold internal forces at t, t-dt and inbetween */ 
solserv_create_vec(&fie,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(fie[i]));

/*--------------------------------------allocate three working vectors */
/*   By optimizing this routine one could live with one or two working */
/*    vectors, I needed three to make things straight-forward and easy */
solserv_create_vec(&work,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(work[i]));
/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to 
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init phase
      is cheap in computation (can be costly in memory)
      There will be no solver call on mass or damping array.
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[0]),
               init);
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[mass_array]),
               &(actsolv->sysarray[mass_array]),
               &work[0],
               &work[1],
               init);
if (damp_array>0)
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[damp_array]),
               &(actsolv->sysarray[damp_array]),
               &work[0],
               &work[1],
               init);
/*----------------- init the assembly for stiffness and for mass matrix */
/*                                           (damping is not assembled) */
init_assembly(actpart,actsolv,actintra,actfield,stiff_array,0);
init_assembly(actpart,actsolv,actintra,actfield,mass_array,0);

/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action,&container);

/*----------------------- call elements to calculate stiffness and mass */
*action = calc_struct_nlnstiffmass;
container.dvec          = NULL;
container.dirich        = NULL;
container.global_numeq  = 0;
container.dirichfacs    = NULL;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);

/*-------------------------------------------- calculate damping matrix */
if (damp_array>0)
{
   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   sdyn->k_damp);

   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[mass_array]),
                   &(actsolv->sysarray[mass_array]),
                   sdyn->m_damp);  
}
/*-------------------------------------- create the original rhs vector */
/*-------------------------- the approbiate action is set inside calrhs */
/*---------------------- this vector holds loads due to external forces */
container.kstep = 0;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[2]),action,&container);
/*------------------------------------------------- copy the rhs vector */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[3]));
/*----------------------- init the time curve applied to the loads here */
/*-------------- this control routine at the moment always uses curve 0 */
/*-------------------------------------------------- init the timecurve */
actcurve = 0;
dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);

/* put a zero the the place 7 in node->sol to init the velocities and accels */
/* of prescribed displacements */
solserv_putdirich_to_dof(actfield,0,0.0,8);

/*---------------------------------- get factor at a certain time t=0.0 */
dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));

/*-------------------------------------- multiply load vector by rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[2]),dynvar.rldfac);

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,dynvar.rldfac,0);


/*----- also put prescribed displacements to the nodes in field sol at  */
/*                                  place 3 separate from the free dofs */
solserv_putdirich_to_dof(actfield,0,dynvar.rldfac,3);

/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[2]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);

sdyn->step = -1;
sdyn->time = 0.0;

/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains(actfield);
}
/*------------------------------------------------------- printout head */
if (par.myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);
/*----------------------------------------------------------------------*/
/*                     START LOOP OVER ALL STEPS                        */
/*----------------------------------------------------------------------*/
/*
   rhs[3]    original load vector
   rhs[2]             load vector at time t-dt
   rhs[1]             load vector at time t
   rhs[0]    interpolated load vector and working array

   fie[2]    internal forces at step t
   fie[1]    internal forces at step t-dt
   fie[0]    interpolated internal forces and working array

   dispi[0]  displacement increment from t-dt to t

   sol[0]    total displacements at time t-dt
   sol[1]    total displacements at time t
   
   vel[0]    velocities    at t-dt
   acc[0]    accelerations at t-dt

   work[2]   working vector for sums and matrix-vector products 
   work[1]   working vector for sums and matrix-vector products 
   work[0]   working vector for sums and matrix-vector products 
   work[0]   is used to hold residual displacements in corrector 
             iteration
             
   in the nodes, displacements are kept in node[].sol[0][0..numdf-1]
                 velocities    are kept in node[].sol[1][0..numdf-1]
                 accelerations are kept in node[].sol[2][0..numdf-1]

Values of the different vectors from above in one loop:
  /......no change in this step
  =,+=...evaluation in this step

   vector	Predictor - Start     Precictor - End	  Corrector - Start	Corrector - End			Update - End

   rhs[3]  	/{=orig. load vect.} 	/	 		/		/				/
   rhs[2]       /{=rhs(t-dt)}     	/			/		/				=rhs[1]{=rhs(t)}
   rhs[1]       =rhs(t)      		/			/		/				/
   rhs[0]    	/{=rhs(t-2dt)}		=feff_p			/		=feff_c				=rhs[2]{=rhs(t-dt)}

   fie[2]    	/			/			=fint(t)	/				/
   fie[1]    	=fint(t-dt)		/			/		/				/
   fie[0]    	/			/			/		=(1-alpha_f)*fie[2]+alpha_f*fie[1]	/

   dispi[0]  	=0			=Keff^-1*feff-p		/		+=work[0]			/

   sol[0]    	/{=d(t-dt)}		/			/		/				=sol[1]{=d(t)}
   sol[1]    	{=d(t-dt)}		=sol[0]+dispi[0]{=d(t)}	/		=sol[0]+dispi[0]		/
   
   vel[0]    	/{=v(t-dt)}		/			/		/				=v(t)
   acc[0]    	/{=a(t-dt)}		/			/		/				=a(t)	

   work[2]    	/{=v(t-2dt)}		/			/		/				=v(t-dt)
   work[1]    	/{=a(t-2dt)}		/			/		/				=a(t-dt)
   work[0]    	/ 			/			/		=Keff^-1*feff-c			=M*vel[0]

*/
timeloop:
t0 = ds_cputime();
/*------------------------------------------------- write memory report */
if (par.myrank==0) dsmemreport();
/*--------------------------------------------------- check for restart */
if (restart)
{
   t0_res = ds_cputime();
   /*-------------- save the stepsize as it will be overwritten in sdyn */
   dt    = sdyn->dt;
   /*------ save the number of steps, as it will be overwritten in sdyn */
   nstep = sdyn->nstep;
   maxtime = sdyn->maxtime;
   /*------------- save the restart interval, as it will be overwritten */
   mod_res_write = sdyn->res_write_evry;
   /*----------------------------------- the step to read in is restart */
   restart_read_nlnstructdyn(restart,
                             sdyn,
                             &dynvar,
                             actfield,
                             actpart,
                             actintra,
                             action,
                             actsolv->nrhs, actsolv->rhs,
                             actsolv->nsol, actsolv->sol,
                             1            , dispi       ,
                             1            , vel         ,
                             1            , acc         ,
                             3            , fie         ,
                             3            , work        ,
                             &intforce_a,
                             &dirich_a,
                             &container);     /* contains variables defined in container.h */
   /*-------------------------------------- put the dt to the structure */
   sdyn->dt = dt;
   /*--------------------------------------- put nstep to the structure */
   sdyn->nstep = nstep;
   sdyn->maxtime = maxtime;
   /*-------------------------------- put restart interval to structure */
   sdyn->res_write_evry = mod_res_write;
   /*------------------------------------------- switch the restart off */
   restart=0;
   /*----------------------------------------------------- measure time */
   t1_res = ds_cputime();
   fprintf(allfiles.out_err,"TIME for restart reading is %f sec\n",t1_res-t0_res);
}
/*--------------------------------------------- increment step and time */
sdyn->step++;
/*------------------- modifications to time steps siye can be done here */
/*------------------------------------------------ set new absolue time */
sdyn->time += sdyn->dt;
/*--- put time to global variable for time-dependent load distributions */
acttime = sdyn->time;
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);

/*---------------------- set incremental displacements dispi[0] to zero */
solserv_zero_vec(&dispi[0]);

/*------------------------- set residual displacements in nodes to zero */
solserv_result_resid(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*----------------------------------------------------------------------*/
/*                     PREDICTOR                                        */
/*----------------------------------------------------------------------*/
/*---------------------- this vector holds loads due to external forces */
solserv_zero_vec(&(actsolv->rhs[1]));
container.kstep = 0;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[1]),action,&container);
/*------------------------------------------------ get factor at time t */
dyn_facfromcurve(actcurve,sdyn->time,&(dynvar.rldfac));

/*------------------------ multiply rhs[1] by actual load factor rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,dynvar.rldfac,0);

/* put the prescribed scaled displacements to the nodes in field sol at */
/*                                  place 4 separate from the free dofs */
solserv_putdirich_to_dof(actfield,0,dynvar.rldfac,4);

/*-------- put presdisplacements(t) - presdisplacements(t-dt) in place 5 */
solserv_adddirich(actfield,0,3,4,5,-1.0,1.0);

/*----- set factors needed for prescribed displacement terms on rhs eff */
/*
dirichfacs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))         
dirichfacs[1] =  (1.0-alpham)*(1.0/beta)/dt                 
dirichfacs[2] =  (1.0-alpham)/(2*beta) - 1                  
dirichfacs[3] = -(1.0-alphaf)*(gamma/beta)/dt               
dirichfacs[4] =  (1.0-alphaf)*gamma/beta - 1                
dirichfacs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)            
dirichfacs[6] = -(1.0-alphaf) or 0                          
dirichfacs[7] =  raleigh damping factor for mass            
dirichfacs[8] =  raleigh damping factor for stiffness       
dirichfacs[9] =  dt     
see phd theses Mok page 165: generalized alfa time integration with prescribed displ.                                    
*/
dirichfacs[0] = -dynvar.constants[0];
dirichfacs[1] =  dynvar.constants[1];
dirichfacs[2] =  dynvar.constants[2];
dirichfacs[3] = -dynvar.constants[3];
dirichfacs[4] =  dynvar.constants[4];
dirichfacs[5] =  dynvar.constants[5];
dirichfacs[6] = -dynvar.constants[6]; 
dirichfacs[9] =  sdyn->dt; 
if (damp_array>0) {
   dirichfacs[7] =  sdyn->m_damp;
   dirichfacs[8] =  sdyn->k_damp;}
else {
   dirichfacs[7] =  0.0;
   dirichfacs[8] =  0.0;}

/*- calculate tangential stiffness/mass and internal forces at time t-dt */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),&(actsolv->sysarray_typ[mass_array]));
amzero(&dirich_a);
amzero(&intforce_a);
*action = calc_struct_nlnstiffmass;
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);
/*---------------------------- store positive internal forces on fie[1] */
solserv_zero_vec(&fie[1]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(fie[1]),intforce,1.0);

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/*---------- subtract internal forces from interpolated external forces */
solserv_add_vec(&(fie[1]),&(actsolv->rhs[0]),-1.0);

/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);

/*--------------------- create effective load vector (rhs[0]-fie[2])eff */
/*
  Peff = rhs[0] - fie[0] 
         + M*(-a1*dispi[0]+a2*vel[0]+a3*acc[0]) 
         + D*(-a4*dispi[0]+a5*vel[0]+a6*acc[0]) (if present)
    
    a1 =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
    a2 = ((1.0-alpham) * (1.0/beta)/(DSQR(dt)))*dt
    a3 =  (1.0-alpham) / (2.0*beta) - 1.0
    a4 =  (1.0-alphaf) * ((gamma/beta)/dt)
    a5 = ((1.0-alphaf) * ((gamma/beta)/dt))*dt - 1.0
    a6 =  (gamma/beta)/2.0 - 1.0) * dt * (1.0-alphaf)
*/
/*----------------------------------------------------------------------*/
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
              mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
/*
  keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
*/  
/*----------------------------------------------------------------------*/
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
              damp_array);
/*------------- call for solution of system dispi[0] = Keff^-1 * rhs[0] */
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[0]),
               init);
/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);
/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/* here put incremental prescribed displacements from sol[5] to sol_increment[0] ? */
/*----------------------------------------------------------------------*/
/*                     PERFORM EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/
itnum=0;
iterloop:
/*----- set factors needed for prescribed displacement terms on rhs eff */
dirichfacs[0] = -dynvar.constants[0];
dirichfacs[1] =  dynvar.constants[1];
dirichfacs[2] =  dynvar.constants[2];
dirichfacs[3] = -dynvar.constants[3];
dirichfacs[4] =  dynvar.constants[4];
dirichfacs[5] =  dynvar.constants[5];
dirichfacs[6] =  0.0; 
dirichfacs[9] =  sdyn->dt; 
if (damp_array>0) {
   dirichfacs[7] =  sdyn->m_damp;
   dirichfacs[8] =  sdyn->k_damp;}
else {
   dirichfacs[7] =  0.0;
   dirichfacs[8] =  0.0;}
/* zero the stiffness matrix and vector for internal forces and dirichlet forces */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),&(actsolv->sysarray_typ[mass_array]));
amzero(&intforce_a);
amzero(&dirich_a);
/* call element routines for calculation of tangential stiffness and intforce */
*action = calc_struct_nlnstiffmass;
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);
/*---------------------------- store positive internal forces on fie[2] */
solserv_zero_vec(&fie[2]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(fie[2]),intforce,1.0);

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/* interpolate internal forces fie[0] = (1-alfaf)fie[2] + alphaf*fie[1] */
solserv_copy_vec(&fie[2],&fie[0]);
solserv_scalarprod_vec(&fie[0],(1.0-sdyn->alpha_f));
solserv_add_vec(&fie[1],&fie[0],sdyn->alpha_f);

/*-- subtract interpolated internal forces from interp. external forces */
solserv_add_vec(&fie[0],&(actsolv->rhs[0]),-1.0);

/*------------------ add dirichlet forces from prescribed displacements */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);
/*--------------------- create effective load vector (rhs[0]-fie[0])eff */
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
              mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
              damp_array);

/*---------------------------------------- solve keff * rsd[0] = rhs[0] */
/* solve for residual displacements to correct incremental displacements*/
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(work[0]),
               &(actsolv->rhs[0]),
               init);

/*-------------------------- return residual displacements to the nodes */
solserv_result_resid(actfield,actintra,&work[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*-- update the incremental displacements by the residual displacements */
solserv_add_vec(&work[0],&dispi[0],1.0);

/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);
/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*----------------------------------------------- check for convergence */
convergence = 0;
dmax        = 0.0;
solserv_vecnorm_euclid(actintra,&(work[0]),&(dynvar.dinorm));
solserv_vecnorm_euclid(actintra,&(dispi[0]),&(dynvar.dnorm));
solserv_vecnorm_Linf(actintra,&(work[0]),&dmax);
if (dynvar.dinorm < sdyn->toldisp ||
    dynvar.dnorm  < EPS14 ||
    (dynvar.dinorm < EPS14 && dmax < EPS12) )
{
   convergence = 1;
}    
else
{
   itnum++;
   if (itnum==sdyn->maxiter) dserror("No convergence in maxiter steps");
   goto iterloop;
}
/*----------------------------------------------------------------------*/
/*                      END OF EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/
/*----------- make temporary copy of actsolv->rhs[2] to actsolv->rhs[0] */
/*                                   (load at t-dt)                     */
/* because in  dyn_nlnstructupd actsolv->rhs[2] is overwritten but is   */
/* still needed to compute energies in dynnle                           */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
/*------------------ update displacements, velocities and accelerations */
dyn_nlnstructupd(actfield,
                 &dynvar,sdyn,actsolv,
                 &(actsolv->sol[0]),   /* total displacements at time t-dt */
                 &(actsolv->sol[1]),   /* total displacements at time t    */
                 &(actsolv->rhs[1]),   /* load vector         at time t    */
                 &(actsolv->rhs[2]),   /* load vector         at time t-dt */
                 &vel[0],              /* velocities          at time t    */
                 &acc[0],              /* accelerations       at time t    */
                 &work[0],             /* working arrays                   */
                 &work[1],             /* working arrays                   */
                 &work[2]);            /* working arrays                   */

/* 
in the nodes the results are stored the following way: 
place 0 holds total displacements of free dofs at time t 
place 1 holds velocities at time t
place 2 holds accels at time t
place 3 holds prescribed displacements at time t-dt 
place 4 holds prescribed displacements at time t
place 5 holds place 4 - place 3
place 6 holds the  velocities of prescribed dofs
place 7 holds the  accels of prescribed dofs
place 8 is working space
*/
/*-------------------------------------- return velocities to the nodes */
solserv_result_total(actfield,actintra, &vel[0],1,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*---------------------------------------velocities for prescribed dofs */
solserv_adddirich(actfield,0,6,0,1,1.0,0.0);
/*------------------------------------------ return accel. to the nodes */
solserv_result_total(actfield,actintra, &acc[0],2,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*-------------------------------------------accel. for prescribed dofs */
solserv_adddirich(actfield,0,7,0,2,1.0,0.0);
/*------------------------------------------ make all types of energies */
dynnle(&dynvar,sdyn,actintra,actsolv,&dispi[0],&fie[1],&fie[2],
       &(actsolv->rhs[1]),&(actsolv->rhs[0]),&work[0]);
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);
dynvar.etot = dynvar.epot + dynvar.ekin;

/*------------------------------- check whether to write results or not */
mod_disp      = sdyn->step % sdyn->updevry_disp;
mod_stress    = sdyn->step % sdyn->updevry_stress;
/*------------------------------- check whether to write restart or not */
mod_res_write = sdyn->step % sdyn->res_write_evry;
/*------------------------------------------ perform stress calculation */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress_file==1 || ioflags.struct_stress_gid==1)
{
   *action = calc_struct_stress;
   container.dvec          = NULL;
   container.dirich        = NULL;
   container.global_numeq  = 0;
   container.dirichfacs    = NULL;
   container.kstep         = 0;
   calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   container.kstep = 0;
   calreduce(actfield,actpart,actintra,action,&container);
}
/*-------------------------------------------- print out results to out */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress_file==1 && ioflags.struct_disp_file==1)
{
  out_sol(actfield,actpart,actintra,sdyn->step,0);
}
/*--------------------------------------------- printout results to gid */
if (par.myrank==0) 
{
   if (mod_disp==0)
   if (ioflags.struct_disp_gid==1)
   {
      out_gid_sol("displacement",actfield,actintra,sdyn->step,0);
      out_gid_sol("velocities",actfield,actintra,sdyn->step,1);
      out_gid_sol("accelerations",actfield,actintra,sdyn->step,2);
   }
   if (mod_stress==0)
   if (ioflags.struct_stress_gid==1)
   {
   out_gid_sol("stress"      ,actfield,actintra,sdyn->step,0);
   }
}
/*-------------------------------------- write restart data to pss file */
if (mod_res_write==0)
restart_write_nlnstructdyn(sdyn,
                           &dynvar,
                           actfield,
                           actpart,
                           actintra,
                           action,
                           actsolv->nrhs, actsolv->rhs,
                           actsolv->nsol, actsolv->sol,
                           1            , dispi       ,
                           1            , vel         ,
                           1            , acc         ,
                           3            , fie         ,
                           3            , work        ,
                           &intforce_a,
                           &dirich_a,
                           &container);     /* contains variables defined in container.h */
/*----------------------------------------------------- print time step */
if (par.myrank==0) dyn_nlnstruct_outstep(&dynvar,sdyn,itnum);

/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
fprintf(allfiles.out_err,"TIME for step %d is %f sec\n",sdyn->step,t1-t0);
/*-------------------------------------- check time and number of steps */
if (sdyn->step < sdyn->nstep-1 && sdyn->time <= sdyn->maxtime)
goto timeloop;
/*----------------------------------------------------------------------*/
end:
/*--------------------------------------------------- cleaning up phase */
amdel(&intforce_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&dispi,1);
solserv_del_vec(&vel,1);
solserv_del_vec(&acc,1);
solserv_del_vec(&fie,3);
solserv_del_vec(&work,3);
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
end_expl:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */
