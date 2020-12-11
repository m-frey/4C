/*----------------------------------------------------------------------*/
/*! \file
\brief strategies for Newton-Raphson convergence check for monolithic scalar-structure interaction
problems

To keep the time integrator class for monolithic scalar-structure interaction problems as plain as
possible, the convergence check for the Newton-Raphson iteration has been encapsulated within
separate strategy classes. Every specific convergence check strategy (e.g., for monolithic
scalar-structure interaction problems involving standard scalar transport or electrochemistry)
computes, checks, and outputs different relevant vector norms and is implemented in a subclass
derived from an abstract, purely virtual interface class.

\level 2


 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_convcheck_strategies.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 11/17 |
 *----------------------------------------------------------------------*/
SSI::SSIMono::ConvCheckStrategyBase::ConvCheckStrategyBase(
    const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
    )
    : itermax_(parameters.get<int>("ITEMAX")),
      itertol_(parameters.sublist("MONOLITHIC").get<double>("CONVTOL")),
      restol_(parameters.sublist("MONOLITHIC").get<double>("ABSTOLRES"))
{
}


/*-----------------------------------------------------------------------*
 | check termination criterion for Newton-Raphson iteration   fang 11/17 |
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyStd::ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  //! compute L2 norm of scalar transport state vector
  double scatradofnorm(0.);
  ssi_mono.ScaTraField()->Phinp()->Norm2(&scatradofnorm);

  //! compute L2 norm of scalar transport increment vector
  double scatraincnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(
          ssi_mono.increment_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&scatraincnorm);

  //! compute L2 norm of scalar transport residual vector
  double scatraresnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.residual_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&scatraresnorm);

  //! compute L2 norm of structural state vector
  double structuredofnorm(0.);
  ssi_mono.StructureField()->Dispnp()->Norm2(&structuredofnorm);

  //! compute L2 norm of structural residual vector
  double structureresnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.residual_, ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureresnorm);

  //! compute L2 norm of structural increment vector
  double structureincnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.increment_, ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureincnorm);

  // safety checks
  if (std::isnan(scatradofnorm) or std::isnan(scatraresnorm) or std::isnan(scatraincnorm) or
      std::isnan(structuredofnorm) or std::isnan(structureresnorm) or std::isnan(structureincnorm))
    dserror("Vector norm is not a number!");
  if (std::isinf(scatradofnorm) or std::isinf(scatraresnorm) or std::isinf(scatraincnorm) or
      std::isinf(structuredofnorm) or std::isinf(structureresnorm) or std::isinf(structureincnorm))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if (scatradofnorm < 1.e-10) scatradofnorm = 1.e-10;
  if (structuredofnorm < 1.e-10) structuredofnorm = 1.e-10;

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|- scatra-res -|- scatra-inc -|- struct-res "
                   "-|- struct-inc -|"
                << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatraresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | "
                << "(       --      , te = " << std::setw(10) << std::setprecision(3)
                << ssi_mono.dtele_ << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatraresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << scatraincnorm / scatradofnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << structureresnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureincnorm / structuredofnorm
                << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
                << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
                << std::endl;
    }

    // convergence check
    if (scatraresnorm <= itertol_ and structureresnorm <= itertol_ and
        scatraincnorm / scatradofnorm <= itertol_ and
        structureincnorm / structuredofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (scatraresnorm < restol_ and structureresnorm < restol_) exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence
  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+"
                << std::endl;
      std::cout << "|      Newton-Raphson method has not converged after a maximum number of "
                << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+"
              << std::endl;
  }

  return exit;
}


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 11/17 |
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElch::ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  //! compute L2 norm of concentration state vector
  double concdofnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.ScaTraField()->Phinp())
      ->Norm2(&concdofnorm);

  //! compute L2 norm of concentration increment vector
  double concincnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.increment_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&concincnorm);

  //! compute L2 norm of concentration residual vector
  double concresnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.residual_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&concresnorm);

  //! compute L2 norm of potential state vector
  double potdofnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.ScaTraField()->Phinp())
      ->Norm2(&potdofnorm);

  //! compute L2 norm of potential increment vector
  double potincnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.increment_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potincnorm);

  //! compute L2 norm of potential residual vector
  double potresnorm(0.);
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.residual_, ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potresnorm);

  //! compute L2 norm of structural state vector
  double structuredofnorm(0.);
  ssi_mono.StructureField()->Dispnp()->Norm2(&structuredofnorm);

  //! compute L2 norm of structural residual vector
  double structureresnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.residual_, ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureresnorm);

  //! compute L2 norm of structural increment vector
  double structureincnorm(0.);
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.increment_, ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureincnorm);

  // safety checks
  if (std::isnan(concdofnorm) or std::isnan(concresnorm) or std::isnan(concincnorm) or
      std::isnan(potdofnorm) or std::isnan(potresnorm) or std::isnan(potincnorm) or
      std::isnan(structuredofnorm) or std::isnan(structureresnorm) or std::isnan(structureincnorm))
    dserror("Vector norm is not a number!");
  if (std::isinf(concdofnorm) or std::isinf(concresnorm) or std::isinf(concincnorm) or
      std::isinf(potdofnorm) or std::isinf(potresnorm) or std::isinf(potincnorm) or
      std::isinf(structuredofnorm) or std::isinf(structureresnorm) or std::isinf(structureincnorm))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if (concdofnorm < 1.e-10) concdofnorm = 1.e-10;
  if (potdofnorm < 1.e-10) potdofnorm = 1.e-10;
  if (structuredofnorm < 1.e-10) structuredofnorm = 1.e-10;

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res "
                   "---|-- pot-inc ---|- struct-res -|- struct-inc -|"
                << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << potresnorm << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | "
                << "(       --      , te = " << std::setw(10) << std::setprecision(3)
                << ssi_mono.dtele_ << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << concincnorm / concdofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << potincnorm / potdofnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureincnorm / structuredofnorm
                << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
                << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
                << std::endl;
    }

    // convergence check
    if (concresnorm <= itertol_ and potresnorm <= itertol_ and structureresnorm <= itertol_ and
        concincnorm / concdofnorm <= itertol_ and potincnorm / potdofnorm <= itertol_ and
        structureincnorm / structuredofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (concresnorm < restol_ and potresnorm < restol_ and structureresnorm < restol_) exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence

  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|                     Newton-Raphson method has not converged after a maximum "
                   "number of "
                << std::setw(2) << itermax_ << " iterations!                     |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+"
              << std::endl;
  }

  return exit;
}
