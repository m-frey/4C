/*----------------------------------------------------------------------*/
/*!
\file strtimint.cpp
\brief Time integration for spatially discretised structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "strtimint_mstep.H"
#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* Map damping input string to enum term */
enum STR::TimInt::DampEnum STR::TimInt::MapDampStringToEnum
(
  const std::string name  //!< identifier
)
{
  if ( (name == "no") or (name == "No") or (name == "NO") )
  {
    return damp_none;
  }
  else if ( (name == "yes") or (name == "Yes") or (name ==  "YES") 
            or (name == "Rayleigh") )
  {
    return damp_rayleigh;
  }
  else if (name == "Material")
  {
    return damp_material;
  }
  else
  {
    dserror("Cannot cope with damping type %s", name.c_str());
    return damp_none;
  }
}

/*----------------------------------------------------------------------*/
/* Map stress input string to enum */
enum STR::TimInt::StressEnum STR::TimInt::MapStressStringToEnum
(
  const std::string name  //!< identifier
)
{
  if ( (name == "cauchy") or (name == "Cauchy") )
  {
    return stress_cauchy;
  }
  else if ( (name == "2pk") or (name == "2PK")
            or (name == "Yes") or (name == "yes") or (name == "YES") )
  {
    return stress_pk2;
  }
  else if ( (name == "No") or (name == "NO") or (name == "No") )
  {
    return stress_none;
  }
  else
  {
    dserror("Cannot handle (output) stress type %s", name.c_str());
    return stress_none;
  }
}

/*----------------------------------------------------------------------*/
/* Map strain input string to enum */
enum STR::TimInt::StrainEnum STR::TimInt::MapStrainStringToEnum
(
  const std::string name  //!< identifier
)
{
  if ( (name == "ea") or (name == "EA") )
  {
    return strain_ea;
  }
  else if ( (name == "gl") or (name == "GL")
            or (name == "Yes") or (name == "yes") or (name == "YES") )
  {
    return strain_gl;
  }
  else if ( (name == "No") or (name == "NO") or (name == "No") )
  {
    return strain_none;
  }
  else
  {
    dserror("Cannot handle (output) strain type %s", name.c_str());
    return strain_none;
  }
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimInt::TimInt
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: discret_(actdis),
  myrank_(actdis->Comm().MyPID()),
  dofrowmap_(actdis->Filled() ? actdis->DofRowMap() : NULL),
  solver_(solver),
  solveradapttol_(Teuchos::getIntegralValue<int>(sdynparams,"ADAPTCONV")==1),
  solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
  output_(output),
  printscreen_(true),  // ADD INPUT PARAMETER
  errfile_(xparams.get<FILE*>("err file")), 
  printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
  printiter_(true),  // ADD INPUT PARAMETER
  writerestartevery_(sdynparams.get<int>("RESTARTEVRY")),
  writestate_((bool) Teuchos::getIntegralValue<int>(ioparams,"STRUCT_DISP")),
  writestateevery_(sdynparams.get<int>("RESEVRYDISP")),
  writestrevery_(sdynparams.get<int>("RESEVRYSTRS")),
  writestress_(MapStressStringToEnum(ioparams.get<std::string>("STRUCT_STRESS"))),
  writestrain_(MapStrainStringToEnum(ioparams.get<std::string>("STRUCT_STRAIN"))),
  writeenergyevery_(sdynparams.get<int>("RESEVRYERGY")),
  energyfile_(),
  damping_(MapDampStringToEnum(sdynparams.get<std::string>("DAMPING"))),
  dampk_(sdynparams.get<double>("K_DAMP")),
  dampm_(sdynparams.get<double>("M_DAMP")),
  conman_(Teuchos::null),
  uzawasolv_(Teuchos::null),
  surfstressman_(Teuchos::null),
  potman_(Teuchos::null),
  time_(Teuchos::null),
  timen_(0.0),
  dt_(Teuchos::null),
  timemax_(sdynparams.get<double>("MAXTIME")),
  stepmax_(sdynparams.get<int>("NUMSTEP")),
  step_(Teuchos::null),
  stepn_(0),
  dirichtoggle_(Teuchos::null),
  invtoggle_(Teuchos::null),
  zeros_(Teuchos::null),
  dis_(Teuchos::null),
  vel_(Teuchos::null),
  acc_(Teuchos::null),
  disn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  stiff_(Teuchos::null),
  mass_(Teuchos::null),
  damp_(Teuchos::null)
{
  // welcome user
  if (myrank_ == 0)
  {
    std::cout << "Welcome to Structural Time Integration " << std::endl;
    std::cout << "     __o__                          __o__       " << std::endl;
    std::cout << "__  /-----\\__                  __  /-----\\__           " << std::endl;
    std::cout << "\\ \\/       \\ \\    |       \\    \\ \\/       \\ \\          " << std::endl;
    std::cout << " \\ |  tea  | |    |-------->    \\ |  tea  | |          " << std::endl;
    std::cout << "  \\|       |_/    |       /      \\|       |_/          " << std::endl;
    std::cout << "    \\_____/   ._                   \\_____/   ._ _|_ /| " << std::endl;
    std::cout << "              | |                            | | |   | " << std::endl;
    std::cout << std::endl;
  }

  // check wether discretisation has been completed
  if (not discret_->Filled())
  {
    dserror("Discretisation is not complete!");
  }


  // time state 
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, 0.0));  // HERE SHOULD BE SOMETHING LIKE (sdynparams.get<double>("TIMEINIT"))
  timen_ = (*time_)[0];  // set target time to initial time
  dt_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, sdynparams.get<double>("TIMESTEP")));
  step_ = 0;

  // output file for energy
  if ( (writeenergyevery_ != 0) and (myrank_ == 0) )
  {
    std::string energyname 
      = DRT::Problem::Instance()->OutputControlFile()->FileName()
      + ".energy";
    energyfile_ = new std::ofstream(energyname.c_str());
    *energyfile_ << "# timestep       time total_energy kinetic_energy internal_energy external_energy" << std::endl;
  }

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //if (not discret_->Filled())
  //{
  //  discret_->FillComplete();
  //}
  //dofrowmap_ = discret_->DofRowMap();

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap_, true);

  // Dirichlet vector
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap_, true);
  // set Dirichlet toggle vector
  {
    Teuchos::ParameterList p;
    p.set("total time", timen_);
    discret_->EvaluateDirichlet(p, zeros_, null, null, dirichtoggle_);
    zeros_->PutScalar(0.0); // just in case of change
  }
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  invtoggle_ = LINALG::CreateVector(*dofrowmap_, false);
  // compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0, *dirichtoggle_, 1.0);

  // displacements D_{n}
  // cout << "we are here" << endl;
  dis_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));

  // displacements D_{n+1} at t_{n+1}
  disn_ = LINALG::CreateVector(*dofrowmap_, true);
  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*dofrowmap_, true);
  // accelerations A_{n+1} at t_{n+1}
  accn_ = LINALG::CreateVector(*dofrowmap_, true);

  // create empty matrices
  stiff_ = Teuchos::rcp(
    new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
  );
  mass_ = Teuchos::rcp(
    new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
  );
  if (damping_ == damp_rayleigh)
  {
    damp_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
    );
  }

  // initialize constraint manager
  conman_ = Teuchos::rcp(new UTILS::ConstrManager(discret_, 
                                                  (*dis_)(0), 
                                                  sdynparams));
  // initialize Uzawa solver
  uzawasolv_ = Teuchos::rcp(new UTILS::UzawaSolver(discret_, 
                                                   *solver_, 
                                                   dirichtoggle_,
                                                   invtoggle_, 
                                                   sdynparams));
  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_->DofRowMap();

  // Check for surface stress conditions due to interfacial phenomena
  {
    vector<DRT::Condition*> surfstresscond(0);
    discret_->GetCondition("SurfaceStress",surfstresscond);
    if (surfstresscond.size())
    {
      surfstressman_ = rcp(new DRT::SurfStressManager(*discret_));
    }
  }
  
  // Check for potential conditions 
  {
    vector<DRT::Condition*> potentialcond(0);
    discret_->GetCondition("Potential",potentialcond);
    if (potentialcond.size())
    {
      potman_ = rcp(new DRT::PotentialManager(Discretization(), *discret_));
    }
  }

  // determine mass, damping and initial accelerations
  DetermineMassDampConsistAccel();

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state
 * and identify consistent accelerations */
void STR::TimInt::DetermineMassDampConsistAccel()
{
  // temporary force vectors in this routine
  Teuchos::RCP<Epetra_Vector> fext 
    = LINALG::CreateVector(*dofrowmap_, true); // external force
  Teuchos::RCP<Epetra_Vector> fint 
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  // overwrite initial state vectors with DirichletBCs
  ApplyDirichletBC((*time_)[0], (*dis_)(0), (*vel_)(0), (*acc_)(0));

  // get external force
  ApplyForceExternal((*time_)[0], (*dis_)(0), (*vel_)(0), fext);
  
  // initialise matrices
  stiff_->Zero();
  mass_->Zero();

  // get initial internal force and stiffness and mass
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_nlnstiffmass");
    // other parameters that might be needed by the elements
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement", zeros_);
    discret_->SetState("displacement", (*dis_)(0));
    if (damping_ == damp_material) discret_->SetState("velocity", (*vel_)(0));
    discret_->Evaluate(p, stiff_, mass_, fint, null, null);
    discret_->ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_ == damp_rayleigh)
  {
    damp_->Add(*stiff_, false, dampk_, 0.0);
    damp_->Add(*mass_, false, dampm_, 1.0);
    damp_->Complete();
  }

  // calculate consistent initial accelerations
  // WE MISS:
  //   - surface stress forces
  //   - potential forces
  {
    Teuchos::RCP<Epetra_Vector> rhs 
      = LINALG::CreateVector(*dofrowmap_, true);
    if (damping_ == damp_rayleigh)
    {
      damp_->Multiply(false, (*vel_)[0], *rhs);
    }
    rhs->Update(-1.0, *fint, 1.0, *fext, -1.0);
    Epetra_Vector rhscopy = Epetra_Vector(*rhs);
    rhs->Multiply(1.0, *invtoggle_, rhscopy, 0.0);
    solver_->Solve(mass_->EpetraMatrix(), (*acc_)(0), rhs, true, true);
  }

  // leave this
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate Dirichlet BC at t_{n+1} */
void STR::TimInt::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc
)
{
  // apply DBCs
  // needed parameters
  ParameterList p;
  p.set("total time", time);  // target time
  
  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_->ClearState();
  discret_->EvaluateDirichlet(p, dis, vel, acc, dirichtoggle_);
  discret_->ClearState();

  // compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0, *dirichtoggle_, 1.0);

  // ciao
  return;
}

/*----------------------------------------------------------------------*/
/* Reset configuration after time step */
void STR::TimInt::ResetStep()
{
  // reset state vectors
  disn_->Update(1.0, (*dis_)[0], 0.0);
  veln_->Update(1.0, (*vel_)[0], 0.0);
  accn_->Update(1.0, (*acc_)[0], 0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_struct_reset_istep");    
    // go to elements
    discret_->Evaluate(p, null, null, null, null, null);
    discret_->ClearState();
  }

  // I am gone
  return;
}

/*----------------------------------------------------------------------*/
/* output to file
 * originally by mwgee 03/07 */
void STR::TimInt::OutputStep()
{
  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (writerestartevery_ and (step_%writerestartevery_ == 0) )
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if ( writestate_ 
       and writestateevery_ and (step_%writestateevery_ == 0)
       and (not datawritten) )
  {
    OutputState(datawritten);
  }

  // output stress & strain
  if ( writestrevery_
       and ( (writestress_ != stress_none)
             or (writestrain_ != strain_none) )
       and (step_%writestrevery_ == 0) )
  {
    OutputStressStrain(datawritten);
  }

  // output energy
  if ( writeenergyevery_ and (step_%writeenergyevery_ == 0) )
  {
    OutputEnergy();
  }

  // what's next?
  return;
}

/*----------------------------------------------------------------------*/
/* write restart
 * originally by mwgee 03/07 */
void STR::TimInt::OutputRestart
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // write restart output, please
  output_->WriteMesh(step_, (*time_)[0]);
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  output_->WriteVector("acceleration", (*acc_)(0));
  //output_->WriteVector("fexternal", fext_);  // CURRENTLY NOT AVAILABLE THINK OF SCENARIO

  // surface stress
  if (surfstressman_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Map> surfrowmap 
      = surfstressman_->GetSurfRowmap();
    Teuchos::RCP<Epetra_Vector> A
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    Teuchos::RCP<Epetra_Vector> con 
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    surfstressman_->GetHistory(A,con);
    output_->WriteVector("Aold", A);
    output_->WriteVector("conquot", con);
  }
  
  // potential forces
  if (potman_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Map> surfrowmap = potman_->GetSurfRowmap();
    Teuchos::RCP<Epetra_Vector> A 
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    potman_->GetHistory(A);
    output_->WriteVector("Aold", A);
  }

  // constraints
  if (conman_->HaveConstraint())
  {
    output_->WriteDouble("uzawaparameter",
                        uzawasolv_->GetUzawaParameter());
  }

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_)
  { 
    printf("====== Restart written in step %d\n", step_);
    fflush(stdout);
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart written in step %d\n", step_);
    fflush(errfile_);
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations
 * originally by mwgee 03/07 */
void STR::TimInt::OutputState
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // write now
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  output_->WriteVector("acceleration", (*acc_)(0));
  //output_->WriteVector("fexternal",fext_);  // CURRENTLY NOT AVAILABLE
  output_->WriteElementData(); 

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* stress calculation and output
 * originally by lw */
void STR::TimInt::OutputStressStrain
(
  bool& datawritten
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  p.set("action", "calc_struct_stress");
  // other parameters that might be needed by the elements
  p.set("total time", (*time_)[0]);
  p.set("delta time", (*dt_)[0]);
  
  // stress
  if (writestress_ == stress_cauchy)
  {
    // output of Cauchy stresses instead of 2PK stresses
    p.set("cauchy", true);
  }
  else
  {
    // this will produce 2nd PK stress ????
    p.set("cauchy", false);
  }
  Teuchos::RCP<std::vector<char> > stressdata
    = Teuchos::rcp(new std::vector<char>());
  p.set("stress", stressdata);

  // strain
  if (writestrain_ == strain_ea)
  {
    p.set("iostrain", "euler_almansi");
  }
  else if (writestrain_ == strain_gl)
  {
    // WILL THIS CAUSE TROUBLE ????
    // THIS STRING DOES NOT EXIST IN SO3
    p.set("iostrain", "green_lagrange");
  }
  else
  {
    p.set("iostrain", "none");
  }
  Teuchos::RCP<std::vector<char> > straindata
    = Teuchos::rcp(new std::vector<char>());
  p.set("strain", straindata);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement", zeros_);
  discret_->SetState("displacement", (*dis_)(0));
  discret_->Evaluate(p, null, null, null, null, null);
  discret_->ClearState();

  // Make new step
  if (not datawritten)
  {
    output_->NewStep(step_, (*time_)[0]);
  }
  datawritten = true;

  // write stress
  if (writestress_ != stress_none)
  {
    std::string stresstext = "";
    if (writestress_ == stress_cauchy)
    {
      stresstext = "gauss_cauchy_stresses_xyz";
    }
    else if (writestress_ == stress_pk2)
    {
      stresstext = "gauss_2PK_stresses_xyz";
    }
    output_->WriteVector(stresstext, *stressdata, 
                         *(discret_->ElementColMap()));
  }

  // write strain
  if (writestrain_ != strain_none)
  {
    std::string straintext = "";
    if (writestrain_ == strain_ea)
    {
      straintext = "gauss_EA_strains_xyz";
    }
    else
    {
      straintext = "gauss_GL_strains_xyz";
    }
    output_->WriteVector(straintext, *straindata,
                         *(discret_->ElementColMap()));
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
/* output system energies */
void STR::TimInt::OutputEnergy()
{
  // internal/strain energy
  double intergy = 0.0;  // total internal energy
  {
    ParameterList p;
    // other parameters needed by the elements
    p.set("action", "calc_struct_energy");
    
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("displacement", (*dis_)(0));
    // get energies
    Teuchos::RCP<Epetra_SerialDenseVector> energies 
      = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, energies);
    discret_->ClearState();
    intergy = (*energies)(0);
  }

  // global calculation of kinetic energy
  double kinergy = 0.0;  // total kinetic energy
  {
    Teuchos::RCP<Epetra_Vector> linmom 
      = LINALG::CreateVector(*dofrowmap_, true);
    mass_->Multiply(false, (*vel_)[0], *linmom);
    linmom->Dot((*vel_)[0], &kinergy);
    kinergy *= 0.5;
  }

  // external energy
  double extergy = 0.0;  // total external energy
  {
    // WARNING: This will only work with dead loads!!!
    Teuchos::RCP<Epetra_Vector> fext = Fext();
    fext->Dot((*dis_)[0], &extergy);
  }

  // total energy
  double totergy = kinergy + intergy - extergy;

  // the output
  if (myrank_ == 0)
  {
    *energyfile_ << " " << std::setw(9) << step_
                 << std::scientific  << std::setprecision(16)
                 << " " << (*time_)[0]
                 << " " << totergy
                 << " " << kinergy
                 << " " << intergy
                 << " " << extergy
                 << std::endl;
  }

  // in god we trust
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void STR::TimInt::ApplyForceExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
  const Teuchos::RCP<Epetra_Vector> vel,  //!< velocity state
  Teuchos::RCP<Epetra_Vector>& fext  //!< external force
)
{
  ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("displacement", dis);
  if (damping_ == damp_material) discret_->SetState("velocity", vel);
  // get load vector
  discret_->EvaluateNeumann(p, *fext);
  discret_->ClearState();

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at state */
void STR::TimInt::ApplyForceStiffInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // residual displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  Teuchos::RCP<Epetra_Vector> fint,  // internal force
  Teuchos::RCP<LINALG::SparseMatrix> stiff  // stiffness matrix
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  const std::string action = "calc_struct_nlnstiff";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement", disi);
  discret_->SetState("displacement", dis);
  if (damping_ == damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, stiff, null, fint, null, null);
  discret_->ClearState();
  
  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force */
void STR::TimInt::ApplyForceInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // incremental displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  Teuchos::RCP<Epetra_Vector> fint  // internal force
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements 
  const std::string action = "calc_struct_internalforce";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement", disi);  // these are incremental
  discret_->SetState("displacement", dis);
  if (damping_ == damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, null, null, fint, null, null);
  discret_->ClearState();
  
  // where the fun starts
  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
void STR::TimInt::Integrate()
{
  // set target time and step
  timen_ = (*time_)[0] + (*dt_)[0];
  stepn_ = step_ + 1;

  // time loop
  while ( (timen_ <= timemax_) and (stepn_ <= stepmax_) )
  {
    // integrate time step
    // after this step we hold disn_, etc
    IntegrateStep();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    UpdateStep();

    // update time and step
    time_->UpdateSteps(timen_);
    step_ = stepn_;
    // 
    timen_ += (*dt_)[0];
    stepn_ += 1;

    // print info about finished time step
    PrintStep();

    // write output
    OutputStep();
  }

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
