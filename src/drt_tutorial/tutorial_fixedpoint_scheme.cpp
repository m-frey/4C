/*---------------------------------------------------------------------*/
/*! \file

\brief student's c++/baci tutorial simple fixed-point problem


\level 2

*/
/*---------------------------------------------------------------------*/

#include "tutorial_fixedpoint_scheme.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"


/// ctor
TUTORIAL::FixedPointScheme::FixedPointScheme()
{
  // print that we are now in the FixedPointScheme tutorial
  PrintTutorialType();

  // in this function the parameters of the toy problem are defined
  ProblemDefinition();

  // initialize member variables
  time_ = DRT::Problem::Instance()->TutorialParams().get<double>("TIMESTEP");
  dt_ = DRT::Problem::Instance()->TutorialParams().get<double>("TIMESTEP");
  Tend_ = DRT::Problem::Instance()->TutorialParams().get<double>("MAXTIME");
  convtol_ = DRT::Problem::Instance()
                 ->TutorialParams()
                 .sublist("FIXED POINT SCHEME")
                 .get<double>("CONVTOL");
  omega_ = DRT::Problem::Instance()
               ->TutorialParams()
               .sublist("FIXED POINT SCHEME")
               .get<double>("RELAX_PARAMETER");

  // initialize displacement state
  disp_ = 0.0;
  // initialize active force
  Fact_ = 0.0;
  // initialize state
  x_ = 0.0;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::FixedPointScheme::PrintTutorialType()
{
  std::cout << "\n YOU CHOSE THE PARTITIONED FIXED POINT TUTORIAL ! \n\n" << std::endl;

  return;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::FixedPointScheme::ProblemDefinition()
{
  std::cout << "  Body 1 (Spring 1)  Body 1 (Spring 2)   " << std::endl;
  std::cout << "          E1                E2           " << std::endl;
  std::cout << "    \\               disp             /   " << std::endl;
  std::cout << "    \\|  /\\    /\\    O--> /\\    /\\   |/   " << std::endl;
  std::cout << "    \\|_/  \\  /  \\___|___/  \\  /  \\__|/   " << std::endl;
  std::cout << "    \\|     \\/       |       \\/      |/   " << std::endl;
  std::cout << "                    |               |/   " << std::endl;
  std::cout << "                    |               |/   " << std::endl;
  std::cout << "                    |      Body 2   |/   " << std::endl;
  std::cout << "                    |        K      |/   " << std::endl;
  std::cout << "                    |    /\\    /\\   |/   " << std::endl;
  std::cout << "                    |___/  \\  /  \\__|/   " << std::endl;
  std::cout << "                            \\/      |/   " << std::endl;
  std::cout << "                   Body 2 generates  " << std::endl;
  std::cout << "                   active force 'Fact'  " << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "      We search the equilibrium displacement 'disp'  " << std::endl;
  std::cout << "      under action of the active force generated by  " << std::endl;
  std::cout << "      Body 2.  " << std::endl;
  std::cout << "\n" << std::endl;

  // The active force is slowly driven up by an "1-cosine" function
  // from 0 to 'Fmax'. This is done in the TimeLoop().
  //
  // Maximum active force
  Fmax_ = 7500.0;
  // Time after which maximum active force 'Fmax' is reached
  Tm_ = 20.0;

  // Stiffness of Body 2
  K_cell_ = 57.0;

  // Stiffness of Body 1
  E_1_ = 1500.0;
  E_2_ = 1400.0;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::FixedPointScheme::TimeLoop()
{
  int step = 1;

  // loop over
  while (Time() <= MaxTime())
  {
    // ramp active force from zero to maximum value
    if (Time() <= Tm_)
      Fact_ = Fmax_ * 0.5 * (1.0 - cos((M_PI / Tm_) * Time()));
    else
      Fact_ = Fmax_;

    // Predict new coupling force
    x_ = InitialGuess();

    std::cout << "\nTIMESTEP " << step << "/" << MaxTime() / Dt() << " time=" << Time() << "/"
              << MaxTime() << " Fact=" << Fact_ << std::endl;

    // update step
    step++;

    // equilibrium iterations for this time step
    IterateFixedPoint();

    // update time n->n+1
    IncrementTime(Dt());
  }

}  // TimeLoop


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
void TUTORIAL::FixedPointScheme::IterateFixedPoint()
{
  // initialize fixed-point scheme to be unconverged
  bool converged = false;

  // initialize iteration counter
  int iter = 0;

  // equilibrium loop
  while (not converged)
  {
    iter++;

    /// get current state
    double Fadh = x_;

    // call Operator1 (solve Body 1) (Force/Neumann-Partitions)
    disp_ = Operator1(Fadh);

    // call Operator2 (solve Body 2) (displacement/Dirichlet-Partition)
    double Fadhnew = Operator2(disp_);

    // calc increment between previous solution and new solution
    double inc = Fadhnew - Fadh;
    // update state state
    x_ += Omega() * inc;

    // check if converged
    converged = ConvergenceCheck(inc, convtol_);


    std::cout << "iter: " << iter << " inc=" << inc << std::endl;

    if (iter == 100) break;
  }

  if (converged)
    std::cout << "Time=" << Time() << " Converged in " << iter << " steps. disp=" << disp_
              << " Fact=" << Fact_ << std::endl;


  return;
}  // iterate fixed point


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
double TUTORIAL::FixedPointScheme::Operator1(double Fadh) { return Fadh / (E_1_ + E_2_); }


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
double TUTORIAL::FixedPointScheme::Operator2(double disp) { return Fact_ - K_cell_ * disp; }


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
bool TUTORIAL::FixedPointScheme::ConvergenceCheck(double inc, double tol)
{
  if (abs(inc) <= tol)
    return true;
  else
    return false;
}


/*-----------------------------------------------------------------------/
/-----------------------------------------------------------------------*/
double TUTORIAL::FixedPointScheme::InitialGuess() { return Fact_ - K_cell_ * disp_; }
