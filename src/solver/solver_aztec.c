#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |  control solver lib AZTEC                             m.gee 9/01     |
 *----------------------------------------------------------------------*/
void solver_az_msr( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _AZ_ARRAY_MSR   *msr_array,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     )
{
#ifdef AZTEC_PACKAGE
int         i;
int         dim;
int         reuse;
int         azname;
AZVAR      *azvar;

double     *dfrom,*dto;

double     *tmpsol;
ARRAY       tmpsol_a;

double     *tmprhs;
ARRAY       tmprhs_a;

#ifdef DEBUG 
dstrc_enter("solver_az_msr");
#endif
/*----------------------------------------------------------------------*/
azvar = actsolv->azvar;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   /*----------- make processor configuration, dependend on parallelism */
   #ifdef PARALLEL
   AZ_set_proc_config(msr_array->proc_config, (MPI_AZComm)(actintra->MPI_INTRA_COMM));
   #else
   AZ_set_proc_config(msr_array->proc_config, AZ_NOT_MPI);
   #endif
   /*------------------------------------- set default value to options */
   AZ_defaults(msr_array->options,msr_array->params);
   /*-------------------------------------- perform check of msr matrix */
   #ifdef DEBUG 
   AZ_check_msr(
                &(msr_array->bindx.a.iv[0]),
                  msr_array->numeq, 
                  msr_array->N_external,
                  AZ_GLOBAL, 
                  msr_array->proc_config
               );
   #endif
   /*----------------------- set options and params from the input file */
   switch(azvar->azsolvertyp)/*--------------------------- set solver */
   {
   case azsolv_CG:
      msr_array->options[AZ_solver] = AZ_cg;
   break;
   case azsolv_GMRES:
      msr_array->options[AZ_solver] = AZ_gmres;
      msr_array->options[AZ_kspace] = azvar->azsub;
   break;
   case azsolv_CGS:
      msr_array->options[AZ_solver] = AZ_cgs;
   break;
   case azsolv_BiCGSTAB:
      msr_array->options[AZ_solver] = AZ_bicgstab;
   break;
   case azsolv_LU:
      msr_array->options[AZ_solver] = AZ_lu;
   break;
   case azsolv_TFQMR:
      msr_array->options[AZ_solver] = AZ_tfqmr;
   break;
   default:
      dserror("No correct solver for Aztec");
   }
   switch(azvar->azprectyp)/*--------------------- set preconditioner */
   {
   case azprec_none:
      msr_array->options[AZ_precond]         = AZ_none;
   break;
   case azprec_ILUT:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_ilut;
      msr_array->params[AZ_ilut_fill]        = azvar->azfill;
   break;
   case azprec_ILU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_ilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_Jacobi:
      msr_array->options[AZ_precond]         = AZ_Jacobi;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_Neumann:
      msr_array->options[AZ_precond]         = AZ_Neumann;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_Least_Squares:
      msr_array->options[AZ_precond]         = AZ_ls;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_SymmGaussSeidel:
      msr_array->options[AZ_precond]         = AZ_sym_GS;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_LU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_lu;
      msr_array->params[AZ_drop]             = azvar->azdrop;
   break;
   case azprec_RILU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_rilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_BILU:
      dserror("Block Preconditioning Bilu cannot be used in MSR format"); 
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_bilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_ICC:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_icc;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   default:
      dserror("No correct preconditioner for Aztec");
   }
   /*---------------------------------------------- set rest of options */
   msr_array->options[AZ_max_iter] = azvar->aziter;
   msr_array->options[AZ_overlap]  = 0;
   msr_array->options[AZ_poly_ord] = azvar->azpoly;
   msr_array->options[AZ_output]   = AZ_none;/*AZ_all;AZ_warnings;AZ_last;300;*/
   msr_array->options[AZ_conv]     = AZ_r0;
   msr_array->params[AZ_tol]       = azvar->aztol;
   msr_array->params[AZ_drop]      = azvar->azdrop;
   /*--------- make backup copy of bindx, as it is permuted in solution */
   am_alloc_copy(&(msr_array->bindx),&(msr_array->bindx_backup));
   /*-------------------------------------- allocate backup copy of val */
   amdef("val_back",&(msr_array->val_backup),msr_array->val.fdim,1,"DV");
   /*----------------------------- set NULL-pointers for Amat and Aprec */
   msr_array->Amat  = NULL;
   msr_array->Aprec = NULL;
   msr_array->ncall=0;
   /* set flag, that this matrix has been initialized and is ready for solve */   
   msr_array->is_init=1;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
/*--------------------------------------------- check the reuse feature */
/* NOTE: This is not multifield yet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/* 
   check for type of field, each field uses an own aztec internal storage 
   name called azname, so if using aztec with 2 or more different fields, 
   the reuse properties, which are stored internal in aztec (that means the 
   preconditioner) don't become mixed up 
   
   structure uses azname=1
   fluid     uses azname=2
   ale       uses azname=3    
*/
switch (actsolv->fieldtyp)
{
case structure: azname=1; break;
case fluid:     azname=2; break;
case ale:       azname=3; break;
default:
   dserror("Unknown type of field");
break;
}
/*--------------------------------------------- check the reuse feature */
if (msr_array->is_factored==0) amcopy(&(msr_array->val),&(msr_array->val_backup));
else                           amcopy(&(msr_array->val_backup),&(msr_array->val));                           
/*----------------------- transform matrix to processor local numbering */
AZ_transform(
             msr_array->proc_config,
             &(msr_array->external),
             msr_array->bindx.a.iv,
             msr_array->val.a.dv,
             msr_array->update.a.iv,
             &(msr_array->update_index),
             &(msr_array->extern_index),
             &(msr_array->data_org),
             msr_array->numeq,
             NULL,
             NULL,
             NULL,
             NULL,
             AZ_MSR_MATRIX
            );
/*------------------------------------ create Aztec structure AZ_MATRIX */
msr_array->Amat = AZ_matrix_create(
                                   msr_array->data_org[AZ_N_internal]+
                                   msr_array->data_org[AZ_N_border]
                                  );
/*---------------------------------attach dmsr-matrix to this structure */
AZ_set_MSR(
           msr_array->Amat, 
           msr_array->bindx.a.iv, 
           msr_array->val.a.dv, 
           msr_array->data_org, 
           0, 
           NULL, 
           AZ_LOCAL
          );
/*--------------------- save number of external components on this proc */            
msr_array->N_external = msr_array->data_org[AZ_N_external];
/*-------------------------------------------------- reorder rhs-vector */            
tmprhs = am_alloc_copy(&(rhs->vec),&(tmprhs_a));
AZ_reorder_vec(
               tmprhs,
               msr_array->data_org,
               msr_array->update_index,
               NULL
              );
/*--------------------------- reorder initial guess and solution-vector */            
AZ_reorder_vec(
               sol->vec.a.dv,
               msr_array->data_org,
               msr_array->update_index,
               NULL
              );
/*----- allocate temporary solution vector large enough for N_externals */
tmpsol = amdef("tmpsol",&tmpsol_a,(msr_array->numeq+msr_array->N_external),1,"DV");
/*--------------------- copy initial guess to temporary solution vector */
/*             (this looks a bit strange, but it's supposed to be fast) */
dfrom = sol->vec.a.dv;
dto   = tmpsol_a.a.dv;
dim   = sol->vec.fdim;
for (i=0; i<dim; i++) *(dto++) = *(dfrom++);
dto   = &(tmpsol_a.a.dv[dim]);
dim   = tmpsol_a.fdim - sol->vec.fdim;
for (i=0; i<dim; i++) *(dto++) = 0.0;
/*--------------------------------------------- check the reuse feature */
msr_array->data_org[AZ_name]=azname;
/*---------------------------------------------------------- first call */
if (msr_array->ncall==0)
{
   msr_array->options[AZ_pre_calc]  = AZ_calc;
   msr_array->options[AZ_keep_info] = 1;
}
/*------------------------------------------------------ not first call */
else
{
   if (msr_array->is_factored==0) 
   {
      msr_array->options[AZ_pre_calc] = AZ_recalc;
   }
   else                           
   {
      msr_array->options[AZ_pre_calc] = AZ_reuse;
      msr_array->options[AZ_keep_info] = 1;
   }
}
/*--------------------------------------------------------- call solver */
AZ_iterate(
           tmpsol,
           tmprhs,
           msr_array->options,
           msr_array->params,
           msr_array->status,
           msr_array->proc_config,
           msr_array->Amat,
           NULL,
           NULL
          );
/*------------------------------------------------ delete temporary rhs */          
amdel(&tmprhs_a);
/*-------------------------------------------- recover unpermuted bindx */
amcopy(&(msr_array->bindx_backup),&(msr_array->bindx));
/*------------------------------------------------ invorder solv vector */
AZ_invorder_vec(
                tmpsol,
                msr_array->data_org,
                msr_array->update_index,
                NULL,
                sol->vec.a.dv
               ); 
/*------------------------------------ delete temporary solution vector */
amdel(&tmpsol_a);
/*---------------------------------------- destroy the Aztec structures */
AZ_matrix_destroy(&(msr_array->Amat)); msr_array->Amat          =NULL;
free(msr_array->external);             msr_array->external      =NULL;
free(msr_array->update_index);         msr_array->update_index  =NULL;
free(msr_array->extern_index);         msr_array->extern_index  =NULL;
free(msr_array->data_org);             msr_array->data_org      =NULL;
/*----------------------------------------- check for success of solver */
if ( (double)(msr_array->status[AZ_why]) != AZ_normal )
{
   if (actintra->intra_rank==0)
   {
      /*----------------------------------------------------- breakdown */
      if ( (double)(msr_array->status[AZ_why]) == AZ_breakdown )
      dserror("Numerical breakdown occured in solver Aztec -> Abort");
      /*------------------------------------numerical loss of precision */
      if ( (double)(msr_array->status[AZ_why]) == AZ_loss )
      {
          printf("RANK 0: AZTEC: Numerical loss of precision occured! continue...\n");
          fprintf(allfiles.out_err,"RANK 0: AZTEC: Numerical loss of precision occured, continue...\n");
      }
      /*------------------------------------------------------ ill cond */
      if ( (double)(msr_array->status[AZ_why]) == AZ_ill_cond )
      {
         printf("RANK 0: AZTEC: Preconditioning ill-conditioned or singular,\n");
         printf("               solution is least square ! continue...\n");
         fprintf(allfiles.out_err,"RANK 0: AZTEC: Preconditioning ill-conditioned or singular,\n");
         fprintf(allfiles.out_err,"               solution is least square ! continue...\n");
      }
      /*-------------------------- maximum number of iterations reached */
      if ( (double)(msr_array->status[AZ_why]) == AZ_maxits )
      {
         printf("RANK 0: AZTEC: Maximum number of iterations %d reached \n",msr_array->options[AZ_max_iter]);
         fprintf(allfiles.out_err,"RANK 0: AZTEC: Maximum number of iterations %d reached \n",msr_array->options[AZ_max_iter]);
         fflush(allfiles.out_err);
         fflush(stdout);
      }
   }
}
/*------------------------------------ print solver iterations and time */             
if (actintra->intra_rank==0)
{
   if (actsolv->fieldtyp==structure) fprintf(allfiles.out_err,"Structure:\n");
   if (actsolv->fieldtyp==fluid)     fprintf(allfiles.out_err,"Fluid:\n");
   if (actsolv->fieldtyp==ale)       fprintf(allfiles.out_err,"Ale:\n");
   fprintf(allfiles.out_err,"AZTEC: %d unknowns %d iterations %f solving time\n",
   sol->numeq_total,
   (int)(msr_array->status[AZ_its]),
   msr_array->status[AZ_solve_time]);
}
/*----------------------------------------------------------- set flags */
msr_array->ncall++;
msr_array->is_factored=1;
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to Aztec");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef AZTEC_PACKAGE */
return;
} /* end of solver_az_msr */




