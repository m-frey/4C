/*----------------------------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_filter.cpp

\brief Inner status test class for constraint problems. Filter
       techniques are based on ideas from multi-objective optimization:

       - Control of the two distinct goals of minimization of the objective
         function and satisfaction of the constraints.

       - Unlike merit functions, filter methods keep these two goals separate

\maintainer Michael Hiermeier

\date Mar 6, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "nox_nln_inner_statustest_filter.H"
#include "nox_nln_inner_statustest_interface_required.H"
#include "nox_nln_meritfunction_lagrangian.H"
#include "nox_nln_linesearch_generic.H"
#include "nox_nln_solver_linesearchbased.H"
#include "nox_nln_statustest_activeset.H"

#include "../drt_lib/drt_dserror.H"

#include <NOX_MeritFunction_Generic.H>
#include <NOX_Solver_Generic.H>

/*----------------------------------------------------------------------------*/
// definition and initialization of static members
unsigned NOX::NLN::INNER::StatusTest::Filter::Point::num_coords_( 0 );
unsigned NOX::NLN::INNER::StatusTest::Filter::Point::num_obj_coords_( 0 );
double NOX::NLN::INNER::StatusTest::Filter::Point::gamma_obj_( 0.0 );
double NOX::NLN::INNER::StatusTest::Filter::Point::gamma_theta_( 0.0 );
std::vector<bool> NOX::NLN::INNER::StatusTest::Filter::Point::isvalid_scaling_;
LINALG::SerialDenseVector NOX::NLN::INNER::StatusTest::Filter::Point::scale_;
LINALG::SerialDenseVector NOX::NLN::INNER::StatusTest::Filter::Point::weights_;

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Filter::Filter(
    const Teuchos::RCP<Generic>& armijo,
    const plain_merit_func_set& infeasibility_vec,
    const double weight_objective_func,
    const double weight_infeasibility_func,
    const double sf,
    const double st,
    const double gamma_alpha,
    const NOX::Utils& utils )
    : status_( status_unevaluated ),
      theta_( infeasibility_vec ),
      curr_points_( plain_const_point_pair( Teuchos::null, Teuchos::null ) ),
      curr_fpoints_( plain_point_pair( Teuchos::null, Teuchos::null ) ),
      filter_(),
      non_dominated_filter_points_(),
      gamma_alpha_( gamma_alpha ),
      amin_obj_( 1.0 ),
      amin_theta_( 1.0 ),
      amin_ftype_( 1.0 ),
      amin_( 1.0 ),
      sf_( sf ),
      st_( st ),
      model_lin_terms_(),
      model_mixed_terms_(),
      armijo_test_( armijo ),
      is_ftype_step_( false ),
      utils_( utils )
{
  Point::resetStaticMembers( 1, theta_.number_, weight_objective_func,
      weight_infeasibility_func );

  Point::setMarginSafetyFactors();

  model_lin_terms_.Resize( Point::num_coords_ );
  model_mixed_terms_.Resize( Point::num_coords_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setMarginSafetyFactors()
{
  gamma_obj_ = std::min<double>( 0.001,
      1.0 / ( 2.0*std::sqrt<double>( static_cast<double>( num_coords_ ) ).real() ) );

  gamma_theta_ = gamma_obj_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::resetStaticMembers(
    const unsigned num_obj_coords,
    const unsigned num_theta_coords,
    const double weight_objective_func,
    const double weight_infeasibility_func )
{
  num_obj_coords_ = num_obj_coords;
  num_coords_ = num_theta_coords + num_obj_coords;

  isvalid_scaling_.resize( num_coords_, false );
  scale_.LightResize( num_coords_ );
  weights_.LightResize( num_coords_ );

  // initialize the scale_ vector
  std::fill( scale_.A(), scale_.A()+num_coords_, 1.0 );

  // initialize the weights vector
  std::fill( weights_.A(), weights_.A()+num_obj_coords_, weight_objective_func );
  std::fill( weights_.A()+num_obj_coords_, weights_.A()+num_coords_, weight_infeasibility_func );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter::Point>
NOX::NLN::INNER::StatusTest::Filter::Point::create(
    const NOX::MeritFunction::Generic& merit_func,
    const Infeasibility& infeasibility_func,
    const NOX::Abstract::Group& grp )
{
  Teuchos::RCP<Point> point_ptr = Teuchos::rcp( new Point );
  Point& point = *point_ptr;

  point(0) = merit_func.computef( grp );
  infeasibility_func.computef( point.A()+Point::num_obj_coords_, grp );

  point.max_theta_id_ = infeasibility_func.findMaxThetaId( point.A()+Point::num_obj_coords_ );

  return point_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter::Point>
NOX::NLN::INNER::StatusTest::Filter::Point::makeFilterPoint( const Point& p )
{
  Teuchos::RCP<Point> fp_ptr = Teuchos::rcp( new Point( p ) );
  Point& fp = *fp_ptr;

  if ( not fp.is_filter_point_ )
  {
    fp.scale();
    fp.setNorm();
    fp.setMargin();
    fp.is_filter_point_ = true;
  }

  return fp_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::scale()
{

  for ( int i=0; i < coords_.Length(); ++i)
  {
    // ToDo Don't use the objective function value after the
    // tangential displacement prediction step as scaling
    // factor for the filter coordinate. This would lead
    // to very high coordinate entries.
//    if (i==0 and (init and pred_ == INPAR::STR::pred_tangdis))
//      continue;

    if (coords_(i) != 0.0)
    {
      if ( not isvalid_scaling_[i] )
      {
        scale_(i) = weights_(i) / std::abs( coords_(i) );
        isvalid_scaling_[i] = true;
      }

      coords_(i) *= scale_(i);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setMargin()
{
  const double max_theta = maxTheta();

  for ( unsigned i=0; i<num_obj_coords_; ++i )
  {
    margin_(i) = gamma_obj_ * max_theta;
  }

  for ( unsigned i=num_obj_coords_; i<num_coords_; ++i )
  {
    margin_(i) = gamma_theta_ * max_theta;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Point::setNorm()
{
  norm_ = coords_.Norm2();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::Filter::Point::print(
    std::ostream& stream, int par_indent_length, const NOX::Utils& u ) const
{
  std::string par_indent;
  par_indent.assign(par_indent_length,' ');

  if ( is_filter_point_ )
  {
    stream << "Filter-Point -- { ";
  }
  else
  {
    stream << "Point -- { ";
  }

  for ( unsigned i=0; i<num_coords_; ++i )
  {
    if ( i!= 0)
      stream << ", ";
    stream << NOX::Utils::sciformat( coords_(i), OUTPUT_PRECISION );
  }
  stream << " } with norm = " << NOX::Utils::sciformat( norm_, OUTPUT_PRECISION ) << ";\n";

  if ( u.isPrintType( NOX::Utils::Details ) )
  {
    stream << par_indent << "margin = { ";
    for ( unsigned i=0; i<num_coords_; ++i )
      {
        if ( i!= 0)
          stream << ", ";
        stream << NOX::Utils::sciformat( margin_(i), OUTPUT_PRECISION );
      }
      stream << " };\n";

      stream << par_indent << "MaxTheta = { "
          << "id = " << max_theta_id_ << ", value = " << maxTheta()
          << ", scale = " << NOX::Utils::sciformat( scaleOfMaxTheta(),
              OUTPUT_PRECISION ) << " };\n";
  }


  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computef(
    double* theta_values,
    const NOX::Abstract::Group& grp ) const
{
  plain_merit_func_set::const_iterator cit = vector_.begin();

  for ( unsigned i=0; cit != vector_.end(); ++cit, ++i )
  {
    const NOX::MeritFunction::Generic& func = **cit;
    theta_values[i] = func.computef( grp );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned NOX::NLN::INNER::StatusTest::Filter::Infeasibility::findMaxThetaId(
    double* theta_values ) const
{
  unsigned max_theta_id = 0;
  unsigned i = 1;

  while ( i < number_ )
  {
   if ( theta_values[max_theta_id] < theta_values[i] )
     max_theta_id = i;

   ++i;
  }

  return max_theta_id;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::CheckStatus(
    const Interface::Required &  interface,
    const NOX::Solver::Generic & solver,
    const NOX::Abstract::Group & grp,
    NOX::StatusTest::CheckType   checkType)
{
  const NOX::NLN::LineSearch::Generic* linesearch_ptr =
      dynamic_cast<const NOX::NLN::LineSearch::Generic*>( &interface );
  if ( not linesearch_ptr )
    dserror( "Dynamic cast failed!" );
  const NOX::NLN::LineSearch::Generic& linesearch = *linesearch_ptr;

  const NOX::Abstract::Vector& dir = linesearch.GetSearchDirection();

  const int iter_ls = interface.GetNumIterations();
  const int iter_newton = solver.getNumIterations();


  const NOX::MeritFunction::Generic& merit_func = interface.GetMeritFunction();

  // do stuff at the beginning of a line search call
  if ( iter_ls == 0 )
  {
    switch ( iter_newton )
    {
      // --------------------------------------------------------------------
      // compute the point coordinates of the first reference state
      // which is accepted by default
      // --------------------------------------------------------------------
      case 0:
      {
        curr_points_.first = Point::create( merit_func, theta_, grp );
        curr_fpoints_.first = Point::makeFilterPoint( *curr_points_.first );
        break;
      }
      // --------------------------------------------------------------------
      // Move the accepted point to the first position at the very beginning
      // of each Newton step (except for the first Newton step)
      // --------------------------------------------------------------------
      default:
      {
        // set accepted  trial point at the first position
        curr_points_.first = curr_points_.second;
        curr_fpoints_.first = curr_fpoints_.second;

        break;
      }
    }

    // setup the linear and quadratic model terms in the beginning
    SetupModelTerms( dir, grp, merit_func );

    // compute the minimal step length estimates
    ComputeMinimalStepLengthEstimates();

    // setup armijo test
    armijo_test_->CheckStatus( interface, solver, grp, checkType );

    status_ = status_unevaluated;

    return status_;
  }

  // reset the f-type flag
  is_ftype_step_ = false;

  // compute the new point coordinates of the trial point
  curr_points_.second = Point::create( merit_func, theta_, grp );
  curr_fpoints_.second = Point::makeFilterPoint( *curr_points_.second );

  // trial filter point
  const Point& trial_fp = *curr_fpoints_.second;
  status_ = AcceptabilityCheck( trial_fp );

  // get current step length
  const double step = linesearch.GetStepLength();

  switch ( status_ )
  {
    // if the current trial point is not in the taboo region, we will check a
    // 2-nd criterion
    case status_converged:
    {
      // ------------------------------------------
      // Final F-Type check
      // ------------------------------------------
      if ( CheckFTypeSwitchingCondition( step ) )
      {
        is_ftype_step_ = true;
        status_ = armijo_test_->CheckStatus( interface, solver, grp, checkType );
      }
      // ------------------------------------------
      // Final filter check
      // ------------------------------------------
      else
      {
        status_ = SufficientReductionCheck( trial_fp );

        if ( status_ == status_converged and trial_fp.norm_ >= 0.0 )
          AugmentFilter();
      }

      break;
    }
    default:
    {
      const enum NOX::StatusTest::StatusType active_set_status =
          GetActiveSetStatus( solver );

      if ( step < gamma_alpha_ * amin_ and active_set_status != NOX::StatusTest::Unconverged )
        dserror("The step-length is too short! We can't find a feasible solution "
            "in the current search direction! (active-set status = %s)",
            NOX::NLN::StatusTest::StatusType2String( active_set_status ).c_str() );

      break;
    }
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::GetActiveSetStatus(
    const NOX::Solver::Generic & solver ) const
{
  const NOX::NLN::Solver::LineSearchBased* ls_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased* >( &solver );
  if ( not ls_solver )
    dserror( "The given non-linear solver is not line search based!" );

  NOX::StatusTest::Generic* active_set_test =
      ls_solver->GetOuterStatusTest<NOX::NLN::StatusTest::ActiveSet>();

  if ( not active_set_test )
    return NOX::StatusTest::Unevaluated;

  return active_set_test->checkStatus( solver, NOX::StatusTest::Complete );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::SufficientReductionCheck(
    const Point& trial_fp ) const
{
  const Point& previous_fp = *curr_fpoints_.first;

  for ( unsigned i = 0; i<Point::num_coords_; ++i )
  {
    if ( trial_fp(i) <= previous_fp(i) - previous_fp.margin_(i) )
    {
      return status_converged;
    }
  }

  return status_step_too_long;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::AugmentFilter()
{
  const Teuchos::RCP<Point> new_fp_ptr = curr_fpoints_.second;

  std::fill( filter_.begin(), filter_.end(), Teuchos::null );
  filter_.resize( non_dominated_filter_points_.size() + 1, Teuchos::null );
  plain_point_set::iterator it = filter_.begin();

  for ( plain_point_set::const_iterator cit = non_dominated_filter_points_.begin();
      cit != non_dominated_filter_points_.end(); ++cit, ++it )
  {
    const Teuchos::RCP<Point>& fp_ptr = *cit;

    if ( new_fp_ptr->norm_ < fp_ptr->norm_ )
    {
      *it = new_fp_ptr;
      ++it;
    }
    *it = *cit;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::AcceptabilityCheck( const Point& trial_fp )
{
  if ( Prefiltering( trial_fp ) )
    return status_converged;

  /* If prefiltering did not succeed, we perform the acceptability check point
   * by point */
  bool passed_check = false;
  for ( plain_point_set::const_iterator cit = filter_.begin();
      cit != filter_.end(); ++cit )
  {
    const Point& fp = **cit;

    passed_check = false;
    for ( unsigned i=0; i<fp.num_coords_; ++i )
    {
      if ( fp.norm_ < 0.0 or
          trial_fp(i) < fp(i) - fp.margin_(i) )
      {
        passed_check = true;
        break;
      }
    }

    // If the check failed for one filter point, the whole test failed.
    if ( not passed_check )
      break;
  }

  return ( passed_check ? status_converged : status_step_too_long );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Filter::Prefiltering(
    const Point& trial_fp )
{
  bool success = false;

  const double sqrt_theta_num = std::sqrt<double>( static_cast<double>( theta_.number_ ) ).real();

  unsigned prefiltering_index = 0;
  for ( plain_point_set::const_iterator cit = filter_.begin();
        cit < filter_.end(); ++cit, ++prefiltering_index )
  {
    const Point& fp = **cit;

    if ( fp.norm_ < 0.0 or
         trial_fp.norm_ < fp.norm_ - sqrt_theta_num * Point::gamma_theta_ * fp.maxTheta() )
    {
      success = true;
      break;
    }
  }

  non_dominated_filter_points_.clear();

  if ( success )
    IdentifyNonDominatedFilterPoints( trial_fp, prefiltering_index );

  return success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::IdentifyNonDominatedFilterPoints(
    const Point& trial_fp,
    const unsigned prefiltering_index )
{
  non_dominated_filter_points_.reserve( filter_.size() );

  std::copy( filter_.begin(), filter_.begin() +  prefiltering_index,
      std::back_inserter( non_dominated_filter_points_ ) );

  for ( plain_point_set::const_iterator cit = filter_.begin() + prefiltering_index;
        cit < filter_.end(); ++cit )
  {
    const Point& fp = **cit;

    unsigned num_dominated_coords = 0;

    if ( fp.norm_ < 0.0 )
      num_dominated_coords = fp.num_coords_;
    else
    {
      for ( unsigned i=0; i<fp.num_coords_; ++i )
      {
        if ( fp(i) - fp.margin_(i) >= trial_fp(i) )
          ++num_dominated_coords;
      }
    }

    // If the trial filter point does not dominate the current filter point, we
    // will keep the current filter point in our filter point set.
    if ( num_dominated_coords < fp.num_coords_ )
      non_dominated_filter_points_.push_back( *cit );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::SetupModelTerms(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp,
    const NOX::MeritFunction::Generic& merit_func )
{
  if ( merit_func.name() == "Lagrangian" )
  {
    const NOX::NLN::MeritFunction::Lagrangian& lagrangian =
        dynamic_cast<const NOX::NLN::MeritFunction::Lagrangian&>( merit_func );

    model_lin_terms_(0) = lagrangian.computeSlope( dir, grp );
    model_mixed_terms_(0) = lagrangian.computeMixed2ndOrderTerms( dir, grp );
  }
  else
    dserror("Currently unsupported merit function type: \"%s\"",
        merit_func.name().c_str() );

  theta_.computeSlope( dir, grp, model_lin_terms_.A()+Point::num_obj_coords_ );
  theta_.computeMixed2ndOrderTerms( dir, grp, model_mixed_terms_.A()+Point::num_obj_coords_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::ComputeMinimalStepLengthEstimates()
{
  /* compute minimal step length estimate based on the 2nd objective function
   * filter acceptability check */
  amin_obj_ = MinimalStepLengthEstimateOfObjFuncFilterCheck();

  /* compute minimal step length estimate based on the 2nd constraint violation
   * filter acceptability check */
  amin_theta_ = theta_.minimalStepLengthEstimate(
      curr_points_.first->A() + Point::num_obj_coords_,
      model_lin_terms_.A() + Point::num_obj_coords_ );

  /* compute minimal step length estimate based on the ftype switching condition */
  amin_ftype_ = 1.0;
  if ( CheckFTypeSwitchingCondition( 1.0 ) )
  {
    amin_ftype_ = MinimalStepLengthEstimateOfFTypeCondition();
  }

  amin_ = std::min( std::min( amin_obj_, amin_theta_ ), amin_ftype_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double
NOX::NLN::INNER::StatusTest::Filter::MinimalStepLengthEstimateOfObjFuncFilterCheck() const
{
  double amin_obj = 1.0;

  // Is the current search direction a descent direction for the objective model?
  if ( model_lin_terms_(0) < 0.0 )
  {
    // get the maximal value of the accepted infeasibility measurements
    const double max_theta = curr_points_.first->maxTheta();
    // get the accepted objective function value
    const double obj_slope = model_lin_terms_(0);

    // check the 2nd order mixed derivative term
    const double obj_mixed_term = model_mixed_terms_(0);
    const bool is_linear_obj_model = ( std::abs( obj_mixed_term ) < 1.0e-12 );

    if ( is_linear_obj_model )
    {
      /*-----------------------------*
       | Filter Check (linear model) |
       *-----------------------------*-------------------------------------*
       |  Linear model:                                                    |
       |  s_f* (L_k + a * LIN(L_k)) < s_f * L_k - s_t * gamma_f * theta_k, |
       |                                                                   |
       |  where s_f and s_t are the scaling factors for the objective      |
       |  and constraint values, respectively.                             |
       |                                                                   |
       | => a > - (s_t/s_f) * (gamma_f * theta_k)/LIN(L_k).                |
       *-------------------------------------------------------------------*/

      amin_obj = - ( Point::gamma_obj_ * max_theta * Point::scale_(1) ) /
          ( obj_slope * Point::scale_(0) );
    }
    else
    {
      /*-------------------------------*
       | Filter Check (2nd order model |
       *-------------------------------*---------------------------------------------*
       |  Quadratic model:                                                           |
       |  s_f * (L_k + c1 * a + c2 * a^2) < s_f * L_k - s_t * gamma_f * theta_k      |
       |                                                                             |
       |  a_1/2 =   (-c1 (+-) sqrt(c1^2 + 4 * c2 * gamma_f * (s_t/s_f) * theta_k))   |
       |          / (-2 * c2).                                                       |
       |  Only the solution corresponding to the minus sign is interesting. To       |
       |  understand this, we consider two different cases. For all of them is       |
       |  c1 lower than zero (descent direction):                                    |
       |                                                                             |
       |  [1] c2 > 0. This corresponds to a parabola which opens upward:             |
       |      In this case there are normally two positive roots and we choose the   |
       |      1st/smaller one. The minimizer of the quadratic 1-D model is not       |
       |      important for us. Nevertheless, it would be possible to check the      |
       |      gradient of the 1-D model for the unity step length and extend the     |
       |      line search method by increasing the step-length if the gradient       |
       |                                                                             |
       |  r_s(x_k+d) - (z_n + dz)^T * grad[wgn(x_k+d)]^T * d - wgn(x_k+d)^T * dz     |
       |                                                                             |
       |      is lower than zero and the step is not accepted. At the moment we use  |
       |      always a backtracking strategy and need only a lower bound for the     |
       |      step length parameter.                                                 |
       |                                                                             |
       |  [2] c2 < 0. This corresponds to a parabola which opens downward. Here we   |
       |      use the same idea. We are interested in a lower bound. We need the     |
       |      2nd/right root, which corresponds to the minus sign again.             |
       *-----------------------------------------------------------------------------*/
      // If the parabola opens upward and the minimum of the quadratic model lies over
      // the specified threshold, it is not possible to find a solution, because the
      // the parabola and the constant line have no intersection point.
      if ( obj_mixed_term > 0.0 and
          ( (obj_slope * obj_slope) / (-4.0 * obj_mixed_term) >
            -Point::scale_(1)/Point::scale_(0) * Point::gamma_obj_ * max_theta ) )
        amin_obj = 1.0;
      else
      {
        double atmin_obj = std::sqrt<double>( obj_slope*obj_slope - 4.0 * obj_mixed_term
            * Point::gamma_obj_ * Point::scale_(1) / Point::scale_(0) * max_theta ).real();
        atmin_obj = ( -obj_slope - atmin_obj)/( 2.0 * obj_mixed_term );

        amin_obj = std::min( amin_obj, std::max( 0.0, atmin_obj ) );
      }
    }
  }

  return amin_obj;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double
NOX::NLN::INNER::StatusTest::Filter::MinimalStepLengthEstimateOfFTypeCondition() const
{
  // linear term of the objective model
  const double obj_slope = model_lin_terms_(0);

  // 2nd order mixed derivative term of the objective function
  const double obj_mixed = model_mixed_terms_(0);

  const double scale_of_max_theta = curr_points_.first->scaleOfMaxTheta();

  // safe-guarding strategy: lower and upper bounds for the minimal step length estimate
  double lBound = 0.0;
  double uBound = 1.0;

  const double d = std::pow( scale_of_max_theta, st_ ) / std::pow( Point::scale_(0), sf_ );

  if ( computeFTypeSwitchingCondition( lBound, d ) > 0.0 )
    dserror("The function value for the lower bound is greater than zero!");

  if ( computeFTypeSwitchingCondition( uBound, d ) < 0.0 )
    dserror("The function value for the upper bound is lower than zero!");

  // set initial value
  double amin = lBound;

  // newton control parameters
  const unsigned itermax = 10;
  unsigned iter = 0;
  bool isconverged = false;

  const double TOL_LOCAL_NEWTON = 1.0e-8;

  while ( not isconverged and iter < itermax )
  {
    double f = computeFTypeSwitchingCondition( amin, d );

    // update lower bound
    if ( f < 0.0 and amin > lBound )
      lBound = amin;
    // update upper bound
    else if ( f > 0.0 and amin < uBound )
      uBound = amin;

    double da = - (obj_slope + obj_mixed * amin );
    da = -f / ( std::pow( da, sf_ ) *
        ( 1.0 - sf_ * obj_mixed * amin / ( da ) ) );

    amin += da;

    // safe-guarding strategy
    if ( amin < lBound or amin > uBound )
      amin = 0.5 * ( lBound + uBound );

    // relative convergence check
    isconverged = std::abs( da ) < TOL_LOCAL_NEWTON * std::max( amin, 1.0e-12 );

    ++iter;
  }
  if ( not isconverged )
    dserror( "The local Newton did not converge! " );

  return amin;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Filter::CheckFTypeSwitchingCondition(
    const double step ) const
{
  const double obj_model = GetObjModel( step );

  // descent direction?
  if ( obj_model < 0.0 )
  {
    const double scale_of_max_theta = curr_points_.first->scaleOfMaxTheta();

    const double d  = std::pow(scale_of_max_theta,st_) / std::pow(Point::scale_(0),sf_);

    return ( computeFTypeSwitchingCondition( step, d ) > 0.0 );
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::computeFTypeSwitchingCondition(
    const double step, const double d ) const
{
  dsassert( d>0.0, "The scaling factor d is smaller than / equal to zero!" );

  const double max_theta = curr_points_.first->maxTheta();

  // linear term of the objective model
  const double obj_slope = model_lin_terms_(0);

  // 2nd order mixed derivative term of the objective function
  const double obj_mixed = model_mixed_terms_(0);

  return std::pow( -( obj_slope + step * obj_mixed ), sf_ ) * step -
      d * std::pow( max_theta, st_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::GetObjModel( const double step ) const
{
  return step * model_lin_terms_(0) + step * step * model_mixed_terms_(0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::Filter::Infeasibility::minimalStepLengthEstimate(
    const double* accepted_theta,
    const double* theta_slope ) const
{
  double amin = 1.0;
  /*------------------*
   | Filter Check     |
   *------------------*------------------------------------*
   |    theta_k + a * LIN(theta_k) < (1-gamma_t) * theta_k |
   |                                                       |
   | => a > -gamma_t*theta_k / LIN(theta_k)                |
   *-------------------------------------------------------*/
  for ( unsigned i=0; i < number_; ++i )
  {
    double amin_theta = 1.0;

    // Is the current search direction a descent direction for the infeasibility measure?
    if ( theta_slope[i] < 0.0 )
    {
      amin_theta = -Point::gamma_theta_ * accepted_theta[i] / theta_slope[i];
    }

    amin = std::min( amin_theta, amin );
  }

  return amin;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computeSlope(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp,
    double * theta_slope_values ) const
{
  plain_merit_func_set::const_iterator cit = vector_.begin();

  for ( unsigned i=0; cit != vector_.end(); ++cit, ++i )
  {
    const NOX::MeritFunction::Generic& func = **cit;

    theta_slope_values[i] = func.computeSlope( dir, grp );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Filter::Infeasibility::computeMixed2ndOrderTerms(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp,
    double * theta_mixed_values ) const
{
  // no mixed 2nd order terms for the infeasibility measures
  std::fill( theta_mixed_values, theta_mixed_values+number_, 0.0 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream & NOX::NLN::INNER::StatusTest::Filter::Print(
    std::ostream& stream, int indent ) const
{
  std::string indent_str;
  indent_str.assign(indent,' ');

  stream << indent_str;
  stream << status_;

  std::string par_indent( "    " );
  par_indent += indent_str;
  const int par_length = par_indent.size();

  stream << indent_str << "CURRENT POINT PAIR \n{\n";
  stream << par_indent << "Accepted previous point:\n";
  curr_points_.first->print( stream, par_length, utils_ );

  stream << par_indent << "Current trial point:\n";
  curr_points_.second->print( stream, par_length, utils_ );
  stream << "}\n";

  stream << indent_str << "FILTER\n{\n";
  unsigned id=0;
  for ( auto cit= filter_.begin(); cit != filter_.end(); ++cit, ++id )
  {
    const Point& fp = **cit;
    stream << "(" << id << ") ";
    fp.print( stream, par_length, utils_ );
  }
  stream << "}\n";

  stream << indent_str << "F-Type condition is "
      << ( is_ftype_step_ ? "" : "not " ) << "fulfilled.\n";
  if ( is_ftype_step_ )
    armijo_test_->Print( stream, par_indent.size() );

  if ( not utils_.isPrintType( NOX::Utils::Details ) )
    return stream;

  // --- Detailed filter output -----------------------------------------------

  stream << indent_str << "MINIMAL STEP LENGTH ESTIMATES\n{\n";
  stream << par_indent << "Objective estimate = "
      << NOX::Utils::sciformat( amin_obj_, OUTPUT_PRECISION ) << "\n";
  stream << par_indent << "Theta estimate     = "
      << NOX::Utils::sciformat( amin_theta_, OUTPUT_PRECISION ) << "\n";
  stream << par_indent << "F-type estimate    = "
      << NOX::Utils::sciformat( amin_ftype_, OUTPUT_PRECISION ) << "\n";
  stream << par_indent << "-------------------- " << "\n";
  stream << par_indent << "Over-all estimate  = "
      << NOX::Utils::sciformat( amin_, OUTPUT_PRECISION ) << "\n";
  stream << "}\n";

  stream << indent_str << "INFEASIBILITY STATISTICS\n{\n";
  stream << par_indent << "Number of theta  = " << theta_.number_ << "\n";
  stream << par_indent << "Types            = {";
  for ( const auto& theta_ptr : theta_.vector_ )
  {
    stream << " \"" << theta_ptr->name() << "\"";
  }
  stream << " };\n";
  stream << "}\n";

  stream << indent_str << "GENERAL POINT STATISTICS\n{\n";
  stream << par_indent << "Number of coords = " << Point::num_coords_ << "\n";
  stream << par_indent << "Number of obj    = " << Point::num_obj_coords_ << "\n";
  stream << par_indent << "Gamma_obj        = "
      << NOX::Utils::sciformat( Point::gamma_obj_, OUTPUT_PRECISION ) << "\n";
  stream << par_indent << "Gamma_theta      = "
      << NOX::Utils::sciformat( Point::gamma_theta_, OUTPUT_PRECISION ) << "\n";
  stream << par_indent << "Scales           = { ";
  for ( unsigned i=0; i<Point::num_coords_; ++i )
  {
    if ( i!= 0)
      stream << ", ";
    stream << NOX::Utils::sciformat( Point::scale_(i), OUTPUT_PRECISION );
  }
  stream << " };\n";
  stream << par_indent << "Valid scales           = {";
  for ( const bool valid_scale : Point::isvalid_scaling_  )
  {
    stream << " ";
    stream << ( valid_scale ? "valid" : "invalid" );
  }
  stream << " };\n";
  stream << par_indent << "Weights          = { ";
  for ( unsigned i=0; i<Point::num_coords_; ++i )
  {
    if ( i!= 0)
      stream << ", ";
    stream << NOX::Utils::sciformat( Point::weights_(i), OUTPUT_PRECISION );
  }
  stream << " };\n";
  stream << "}\n";

  return stream;
}
