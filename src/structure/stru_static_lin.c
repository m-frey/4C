/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];


/*----------------------------------------------------------------------*
 |  routine to control static execution                  m.gee 6/01     |
 *----------------------------------------------------------------------*/
void calsta()
{
#ifdef DEBUG
dstrc_enter("calsta");
#endif
/*----------------------------------------------------------------------*/

if (statvar->linear==1 && statvar->nonlinear==1)
   dserror("linear and nonlinear static analysis on");

if (statvar->linear==1)
{
   stalin();
}
if (statvar->nonlinear==1)
{
   stanln();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of calsta */




/*----------------------------------------------------------------------*
 |  routine to control linear static structural analysis    m.gee 6/01  |
 *----------------------------------------------------------------------*/
void stalin()
{
INT           i;                /* a counter */
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

static ARRAY   dirich_a;
static DOUBLE *dirich;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */

CONTAINER     container;        /* contains variables defined in container.h */

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

#ifdef BINIO
BIN_OUT_FIELD out_context;
#endif

INT           disnum = 0;

container.isdyn   = 0;           /* static calculation */
container.kintyp  = 0;           /* kintyp  = 0: geo_lin*/
container.disnum  = disnum;           /* disnum  = 0: only one discretisation*/

#ifdef DEBUG
dstrc_enter("stalin");
#endif
/*----------------------------------------------------------------------*/
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
container.fieldtyp  = actfield->fieldtyp;
/*----------------------------------------------------------------------*/
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
#ifdef PARALLEL
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/*    intracommunicator (in case of linear statics, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=2;
actsolv->nsol=2;
solserv_create_vec(&(actsolv->rhs),2,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),2,numeq_total,numeq,"DV");
dirich = amdef("dirich",&dirich_a,numeq_total,1,"DV");
amzero(&dirich_a);
/*------------------------------ init the created dist. vectors to zero */
for (i=0; i<actsolv->nrhs; i++)
   solserv_zero_vec(&(actsolv->rhs[i]));
for (i=0; i<actsolv->nsol; i++)
   solserv_zero_vec(&(actsolv->sol[i]));
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[actsysarray]),
                  &(actsolv->rhs[actsysarray]),
                    init
                 );
/*--------------------------------- init the dist sparse matrix to zero */
/*               NOTE: Has to be called after solver_control(init=1) */
solserv_zero_mat(
                    actintra,
                    &(actsolv->sysarray[actsysarray]),
                    &(actsolv->sysarray_typ[actsysarray])
                   );
/*----------------------------- init the assembly for ONE sparse matrix */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action,&container);

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   &(actsolv->sysarray_typ[actsysarray]), &(actsolv->sysarray[actsysarray]),
                   actfield, actpart, actintra, 0);
#endif

/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0 && ioflags.output_gid==1)
   out_gid_msh();

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,0,0,0.0);

/*------call element routines to calculate & assemble stiffness matrice */
*action = calc_struct_linstiff;
container.dvec         = NULL;
container.dirich       = dirich;
container.global_numeq = numeq_total;
container.kstep        = 0;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
/*----------------------------------- call rhs-routines to assemble rhs */
/*------------------------------------ set action before call of calrhs */
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,actsysarray,
       &(actsolv->rhs[actsysarray]),action,&container);
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
	     &(actsolv->rhs[actsysarray]),dirich,-1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[actsysarray]),
                  &(actsolv->rhs[actsysarray]),
                    init
                 );
/*-------------------------allreduce the result and put it to the nodes */
solserv_result_total(
                     actfield,
                     disnum,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*------------------------------------------ perform stress calculation */
if (ioflags.struct_stress==1)
{
   *action = calc_struct_stress;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.kstep        = 0;
   calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   container.kstep = 0;
   calreduce(actfield,actpart,disnum,actintra,action,&container);
}


/* printout results to out */
if (ioflags.output_out==1 && ioflags.struct_disp==1)
{
  out_sol(actfield,actpart,disnum,actintra,0,0);
}


/* printout results to binary file */
#ifdef BINIO
if (ioflags.output_bin==1)
{
  if (ioflags.struct_disp==1) {
    out_results(&out_context, 0, 0, 0, OUTPUT_DISPLACEMENT);

#ifdef D_AXISHELL
    out_results(&out_context, 0, 0, 0, OUTPUT_THICKNESS);
    out_results(&out_context, 0, 0, 0, OUTPUT_AXI_LOADS);
#endif
  }

  if (ioflags.struct_stress==1) {
    out_results(&out_context, 0, 0, 0, OUTPUT_STRESS);
  }
}
#endif


/* printout results to gid */
if (ioflags.output_gid==1 && ioflags.struct_disp==1 && par.myrank==0)
{
   out_gid_sol("displacement",actfield,disnum,actintra,0,0,ZERO);
   out_gid_domains(actfield, disnum);
#ifdef D_AXISHELL
   out_gid_sol("thickness",actfield,disnum,actintra,0,0,ZERO);
   out_gid_sol("axi_loads",actfield,disnum,actintra,0,0,ZERO);
#endif
}

/* printout stress to gid */
if (ioflags.output_gid==1 && ioflags.struct_stress==1 && par.myrank==0)
{
   out_gid_sol("stress"      ,actfield,disnum,actintra,0,0,ZERO);
}


end:

#ifdef BINIO
destroy_bin_out_field(&out_context);
#endif

#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of stalin */
