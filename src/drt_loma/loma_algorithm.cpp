/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.cpp

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "loma_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::Algorithm(
    Epetra_Comm&                  comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false)
{
  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_    = prbdyn.get<double>("CONVTOL");
  itmaxpre_ = prbdyn.get<int>("ITEMAX");

  // flag for constant thermodynamic pressure
  consthermpress_ = prbdyn.get<string>("CONSTHERMPRESS");

  // thermodynamic pressure and specific gas constant (default: 98100.0/287.0)
  // in case of non-constant thermodynamic pressure: initial value
  thermpress_  = prbdyn.get<double>("THERMOPRESS");
  gasconstant_ = prbdyn.get<double>("GASCONSTANT");

  // factor for equation of state (i.e., therm. press. / gas constant)
  // constant if thermodynamic pressure is constant
  eosfac_ = thermpress_/gasconstant_;

  // initialization of variable for initial total mass
  // (required if thermodynamic pressure is calculated from mass conservation)
  initialmass_ = 0.0;

  // flag for special flow and start of sampling period from fluid parameter list
  special_flow_ = prbdyn.get<string>("CANONICAL_FLOW");
  samstart_     = prbdyn.get<int>("SAMPLING_START");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeLoop()
{
  // decide whether stationary loop or particular time loop for solving the
  // low-Mach-number flow problem based on FLUID time-integration scheme
  if (FluidField().TimIntScheme() == timeint_stationary)
    SolveStationaryProblem();
  else if (FluidField().TimIntScheme() == timeint_afgenalpha or
           FluidField().TimIntScheme() == timeint_gen_alpha)
    GenAlphaTimeLoop();
  else OSTBDF2TimeLoop();

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
    cout<<"\n**********************\n STATIONARY LOW-MACH-NUMBER FLOW SOLVER \n**********************\n";

  // prepare time step (using one-step-theta/BDF2 procedure)
  OSTBDF2PrepareTimeStep();

  // do outer iteration loop (using one-step-theta/BDF2 procedure)
  OSTBDF2OuterLoop();

  // write output to screen and files
  Output();

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::GenAlphaTimeLoop()
{
  // do initial calculations
  InitialCalculations();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    GenAlphaPrepareTimeStep();

    // do outer iteration loop
    GenAlphaOuterLoop();

    // update all single field solvers
    GenAlphaUpdate();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OSTBDF2TimeLoop()
{
  // do initial calculations
  InitialCalculations();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    OSTBDF2PrepareTimeStep();

    // do outer iteration loop
    OSTBDF2OuterLoop();

    // update all single field solvers
    OSTBDF2Update();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::InitialCalculations()
{
  // compute initial density field using initial temperature + therm. pressure
  // temperature stored at phinp in SCATRA is used -> densnp is set
  ScaTraField().ComputeDensity(thermpress_,gasconstant_);

  // initially set density at 0 (i.e., densn)
  // furthermore, set density at -1 (i.e., densnm) for BDF2 (zero vector)
  ScaTraField().UpdateDensity();

  // store temperature and velocity of previous iteration for convergence check
  ScaTraField().TempIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);
  //ScaTraField().VelIncNp()->Update(1.0,*ConVel(),0.0);

  // compute initial convective velocity field for scalar transport solver
  // using initial fluid velocity (and pressure) field
  // (For generalized-alpha time-integration scheme, velocity at n+1
  // is weighted by density at n+alpha_F, which is identical to density
  // at n+1, since density at n was set equal to n+1 above.)
  // store temperature and velocity of previous iteration for convergence check
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 FluidField().SgVelVisc(),
                                 FluidField().Discretization());

  // set initial value of thermodynamic pressure in SCATRA
  ScaTraField().SetInitialThermPressure(thermpress_);

  if (consthermpress_=="No_energy")
  {
    // compute initial time derivative of thermodynamic pressure
    ScaTraField().ComputeInitialThermPressureDeriv();
    eosfac_ = thermpress_/gasconstant_;
  }
  else if (consthermpress_=="No_mass")
  {
    // compute initial mass in flow domain which is to be conserved
    initialmass_ = ScaTraField().ComputeInitialMass(thermpress_);
    // compute intial thermodynamic pressure and time derivative
    thermpress_  = ScaTraField().ComputeThermPressureFromMassCons(initialmass_,gasconstant_);
    eosfac_ = thermpress_/gasconstant_;
  }

  // write initial fields
  Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::GenAlphaPrepareTimeStep()
{
  // prepare temperature time step (+ initialize one-step-theta and
  // generalized-alpha schemes correctly by computing initial time
  // derivatives of temperature dphi/dt_(0) in first time step)
  ScaTraField().PrepareTimeStep();

  // compute initial time derivative of density in first time step for
  // one-step-theta and generalized-alpha schemes, for which initial
  // time derivatives of temperature dphi/dt_(0) are required
  if (Step() == 1) ScaTraField().ComputeInitialDensityDerivative();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy") ScaTraField().PredictThermPressure();

  // predict density field and time derivative
  ScaTraField().PredictDensity();

  // set density at n+1 and n, density time derivative at n, SCATRA trueresidual
  // and eos factor
  FluidField().SetTimeLomaFields(ScaTraField().DensNp(),
                                 ScaTraField().DensN(),
                                 ScaTraField().DensDtN(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().NumScal(),
                                 eosfac_);

  // prepare fluid time step, particularly predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OSTBDF2PrepareTimeStep()
{
  // prepare temperature time step (+ initialize one-step-theta and
  // generalized-alpha schemes correctly by computing initial time
  // derivatives of temperature dphi/dt_(0) in first time step)
  ScaTraField().PrepareTimeStep();

  // compute initial time derivative of density in first time step for
  // one-step-theta and generalized-alpha schemes, for which initial
  // time derivatives of temperature dphi/dt_(0) are required
  if (Step() == 1) ScaTraField().ComputeInitialDensityDerivative();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy") ScaTraField().PredictThermPressure();

  // predict density field and time derivative
  ScaTraField().PredictDensity();

  // set density at n+1, n and n-1, SCATRA trueresidual and eos factor
  FluidField().SetTimeLomaFields(ScaTraField().DensNp(),
                                 ScaTraField().DensN(),
                                 ScaTraField().DensNm(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().NumScal(),
                                 eosfac_);

  // prepare fluid time step, particularly predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::GenAlphaOuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
    cout<<"\n******************************************\n  OUTER GENERALIZED-ALPHA ITERATION LOOP\n******************************************\n";

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  while (stopnonliniter==false)
  {
    itnum++;

    // compute values at intermediate time steps
    ScaTraField().ComputeIntermediateValues();

    // store temperature and velocity of previous iteration for convergence check
    ScaTraField().TempIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);
    //ScaTraField().VelIncNp()->Update(1.0,*ConVel(),0.0);

    // set field vectors: velocity, subgrid viscosity, (negative) fluid trueresidual
    ScaTraField().SetVelocityField(FluidField().Velaf(),
                                   FluidField().SgVelVisc(),
                                   FluidField().Discretization());

    // solve transport equation for temperature
    if (Comm().MyPID()==0) cout<<"\n******************************************\n   GENERALIZED-ALPHA TEMPERATURE SOLVER\n******************************************\n";
    ScaTraField().Solve();

    // in case of non-constant thermodynamic pressure: compute
    if (consthermpress_=="No_energy")
    {
      thermpress_ = ScaTraField().ComputeThermPressure();
      eosfac_ = thermpress_/gasconstant_;
    }
    else if (consthermpress_=="No_mass")
    {
      thermpress_ = ScaTraField().ComputeThermPressureFromMassCons(initialmass_,gasconstant_);
      eosfac_ = thermpress_/gasconstant_;
    }

    // compute density using current temperature + thermodynamic pressure
    ScaTraField().ComputeDensity(thermpress_,gasconstant_);

    // compute time derivative of density
    ScaTraField().ComputeDensityDerivative();

    // set density (and density time derivative) at n+1 as well as eos factor
    FluidField().SetIterLomaFields(ScaTraField().DensNp(),
                                   ScaTraField().DensDtNp(),
                                   ScaTraField().NumScal(),
                                   eosfac_);

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n******************************************\n      GENERALIZED-ALPHA FLOW SOLVER\n******************************************\n";
    FluidField().MultiCorrector();

    // check convergence of temperature field
    stopnonliniter = ScaTraField().LomaConvergenceCheck(itnum,itmax_,ittol_);
  }

  // compute values at intermediate time steps
  ScaTraField().ComputeIntermediateValues();

  // set field vectors: velocity, subgrid viscosity, (negative) fluid trueresidual
  ScaTraField().SetVelocityField(FluidField().Velaf(),
                                 FluidField().SgVelVisc(),
                                 FluidField().Discretization());

  // solve transport equation for temperature
  if (Comm().MyPID()==0) cout<<"\n******************************************\n   GENERALIZED-ALPHA TEMPERATURE SOLVER\n******************************************\n";
  ScaTraField().Solve();

  // in case of non-constant thermodynamic pressure: compute
  if (consthermpress_=="No_energy")
  {
    thermpress_ = ScaTraField().ComputeThermPressure();
    eosfac_ = thermpress_/gasconstant_;
  }
  else if (consthermpress_=="No_mass")
  {
    thermpress_ = ScaTraField().ComputeThermPressureFromMassCons(initialmass_,gasconstant_);
    eosfac_ = thermpress_/gasconstant_;
  }

  // compute density using current temperature + thermodynamic pressure
  ScaTraField().ComputeDensity(thermpress_,gasconstant_);

  // compute time derivative of density
  ScaTraField().ComputeDensityDerivative();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OSTBDF2OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
    cout<<"\n******************************************\n OUTER ONE-STEP-THETA/BDF2 ITERATION LOOP\n******************************************\n";

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step
  //FluidField().Predictor();

  while (stopnonliniter==false)
  {
    itnum++;

    // store temperature and velocity of previous iteration for convergence check
    ScaTraField().TempIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);
    //ScaTraField().VelIncNp()->Update(1.0,*ConVel(),0.0);

    // set field vectors: velocity, subgrid viscosity, (negative) fluid trueresidual
    ScaTraField().SetVelocityField(FluidField().Velnp(),
                                   FluidField().SgVelVisc(),
                                   FluidField().Discretization());

    // solve transport equation for temperature
    if (Comm().MyPID()==0) cout<<"\n******************************************\n  ONE-STEP-THETA/BDF2 TEMPERATURE SOLVER\n******************************************\n";
    ScaTraField().Solve();

    // in case of non-constant thermodynamic pressure: compute
    if (consthermpress_=="No_energy")
    {
      thermpress_ = ScaTraField().ComputeThermPressure();
      eosfac_ = thermpress_/gasconstant_;
    }
    else if (consthermpress_=="No_mass")
    {
      thermpress_ = ScaTraField().ComputeThermPressureFromMassCons(initialmass_,gasconstant_);
      eosfac_ = thermpress_/gasconstant_;
    }

    // compute density using current temperature + thermodynamic pressure
    ScaTraField().ComputeDensity(thermpress_,gasconstant_);

    // set density (and density time derivative) at n+1 as well as eos factor
    FluidField().SetIterLomaFields(ScaTraField().DensNp(),
                                   ScaTraField().DensDtNp(),
                                   ScaTraField().NumScal(),
                                   eosfac_);

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n******************************************\n     ONE-STEP-THETA/BDF2 FLOW SOLVER\n******************************************\n";
    FluidField().MultiCorrector();

    // check convergence of temperature field
    stopnonliniter = ScaTraField().LomaConvergenceCheck(itnum,itmax_,ittol_);
  }

  // set field vectors: velocity, subgrid viscosity, (negative) fluid trueresidual
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 FluidField().SgVelVisc(),
                                 FluidField().Discretization());

  // solve transport equation for temperature
  if (Comm().MyPID()==0) cout<<"\n******************************************\n  ONE-STEP-THETA/BDF2 TEMPERATURE SOLVER\n******************************************\n";
  ScaTraField().Solve();

  // in case of non-constant thermodynamic pressure: compute
  if (consthermpress_=="No_energy")
  {
    thermpress_ = ScaTraField().ComputeThermPressure();
    eosfac_ = thermpress_/gasconstant_;
  }
  else if (consthermpress_=="No_mass")
  {
    thermpress_ = ScaTraField().ComputeThermPressureFromMassCons(initialmass_,gasconstant_);
    eosfac_ = thermpress_/gasconstant_;
  }

  // compute density using current temperature + thermodynamic pressure
  ScaTraField().ComputeDensity(thermpress_,gasconstant_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::GenAlphaUpdate()
{
  // update temperature
  ScaTraField().Update();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    ScaTraField().UpdateThermPressure();

  // update density
  ScaTraField().UpdateDensity();

  // update fluid
  FluidField().Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OSTBDF2Update()
{
  // set density at n+1, n and n-1, SCATRA trueresidual and eos factor
  // -> required for computing acceleration of present time step in FLUID
  FluidField().SetTimeLomaFields(ScaTraField().DensNp(),
                                 ScaTraField().DensN(),
                                 ScaTraField().DensNm(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().NumScal(),
                                 eosfac_);

  // update temperature
  ScaTraField().Update();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    ScaTraField().UpdateThermPressure();

  // update density
  ScaTraField().UpdateDensity();

  // update fluid
  FluidField().Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Output()
{
  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();

  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
