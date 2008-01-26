/*----------------------------------------------------------------------*/
/*!
\file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "fsi_fluid.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::Fluid::Fluid(Teuchos::RCP<DRT::Discretization> dis,
                  Teuchos::RCP<LINALG::Solver> solver,
                  Teuchos::RCP<ParameterList> params,
                  Teuchos::RCP<IO::DiscretizationWriter> output)
  : FluidImplicitTimeInt(dis,*solver,*params,*output,true),
    interface_(dis),
    meshmap_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  stepmax_    = params_->get<int>   ("max number timesteps");
  maxtime_    = params_->get<double>("total time");
  theta_      = params_->get<double>("theta");
  timealgo_   = params_->get<FLUID_TIMEINTTYPE>("time int algo");
  dtp_ = dta_ = params_->get<double>("time step size");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  relax_ = LINALG::CreateVector(*dofrowmap,true);
  griddisp_ = LINALG::CreateVector(*dofrowmap,true);

  interface_.Setup(DRT::UTILS::CondAnd(DRT::UTILS::ExtractorCondMaxPos(genprob.ndim),
                                       DRT::UTILS::ExtractorCondInCondition(dis,"FSICoupling")));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::Fluid::Itemax() const
{
  return params_->get<int>("max nonlin iter steps");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::SetItemax(int itemax)
{
  params_->set<int>("max nonlin iter steps", itemax);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::SetInterfaceMap(Teuchos::RCP<Epetra_Map> im)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Fluid::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(trueresidual_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,velnp_);

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,dirichtoggle_);

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::SetMeshMap(Teuchos::RCP<Epetra_Map> mm)
{
  meshmap_.SetupMaps(Teuchos::rcp(discret_->DofRowMap(),false),mm);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
{
  meshmap_.InsertCondVector(fluiddisp,dispnp_);

  // new grid velocity
  UpdateGridv();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Fluid::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
{
  meshmap_.InsertCondVector(gridvel,gridv_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Fluid::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  //
  // This method is really stupid, but simple. We calculate the fluid
  // elements twice here. First because we need the global matrix to
  // do a linear solve. We do not have any RHS other than the one from
  // the Dirichlet condition at the FSI interface.
  //
  // After that we recalculate the matrix so that we can get the
  // reaction forces at the FSI interface easily.
  //
  // We do more work that required, but we do not need any special
  // element code to perform the steepest descent calculation. This is
  // quite a benefit as the special code in the old discretization was
  // a real nightmare.
  //

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  relax_->PutScalar(0.0);
  interface_.InsertCondVector(ivel,relax_);

  // set the grid displacement independent of the trial value at the
  // interface
  griddisp_->Update(1., *dispnp_, -1., *dispn_, 0.);

  // dirichtoggle_ has already been set up

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  ParameterList eleparams;
  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);
  eleparams.set("include reactive terms for linearisation",newton_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist"  ,zeros_);
#if 0
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);
#else
  discret_->SetState("dispnp", griddisp_);
  discret_->SetState("gridv", gridv_);
#endif

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  LINALG::Complete(*sysmat_);

  // No, we do not want to have any rhs. There cannot be any.
  residual_->PutScalar(0.0);

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual discplacements are supposed to be zero at
  //          boundary conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,relax_,dirichtoggle_);

  //-------solve for residual displacements to correct incremental displacements
  solver_->Solve(sysmat_,incvel_,residual_,true,true);

  // and now we need the reaction forces

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);

  // set vector values needed by elements
  discret_->ClearState();
  //discret_->SetState("velnp",incvel_);
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist"  ,zeros_);
#if 0
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);
#else
  discret_->SetState("dispnp", griddisp_);
  discret_->SetState("gridv", gridv_);
#endif

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  double density = eleparams.get("density", 0.0);

//   double norm;
//   trueresidual_->Norm2(&norm);
//   if (trueresidual_->Map().Comm().MyPID()==0)
//     cout << "==> residual norm = " << norm << " <==\n";

  // finalize the system matrix
  LINALG::Complete(*sysmat_);

  if (sysmat_->Apply(*incvel_, *trueresidual_)!=0)
    dserror("sysmat_->Apply() failed");
  trueresidual_->Scale(-density/dta_/theta_);

  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Fluid::IntegrateInterfaceShape()
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","integrate_Shapefunction");

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

  // call loop over elements
  discret_->ClearState();
  discret_->SetState("dispnp", dispnp_);
  discret_->EvaluateCondition(eleparams,integratedshapefunc,"FSICoupling");
  discret_->ClearState();

  return interface_.ExtractCondVector(integratedshapefunc);
}

#endif
