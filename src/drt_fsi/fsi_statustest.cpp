
#ifdef CCADISCRET

#include "fsi_statustest.H"
#include "../drt_lib/drt_dserror.H"

#include <NOX_Common.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Utils.H>

#include <NOX_Epetra_Vector.H>

#include <Thyra_DefaultProductVector.hpp>

#include "../drt_lib/linalg_utils.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GenericNormF::GenericNormF(std::string name,
                                 double tolerance,
                                 ScaleType stype)
  : status_(NOX::StatusTest::Unevaluated),
    normType_(NOX::Abstract::Vector::TwoNorm),
    scaleType_(stype),
    toleranceType_(Absolute),
    specifiedTolerance_(tolerance),
    initialTolerance_(1.0),
    trueTolerance_(tolerance),
    normF_(0.0),
    name_(name)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::computeNorm(const Epetra_Vector& v)
{
  int n = v.GlobalLength();
  double norm;
  int err;

  switch (normType_)
  {
  case NOX::Abstract::Vector::TwoNorm:
    err = v.Norm2(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= sqrt(1.0 * n);
    break;

  case NOX::Abstract::Vector::OneNorm:
    err = v.Norm1(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= n;
    break;

  case NOX::Abstract::Vector::MaxNorm:
    err = v.NormInf(&norm);
    if (err!=0)
      dserror("norm failed");
    if (scaleType_ == Scaled)
      norm /= n;
    break;

  default:
    dserror("norm type confusion");
    break;
  }

  return norm;
}


#if 0
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::GenericNormF::relativeSetup(NOX::Abstract::Group& initialGuess)
{
  NOX::Abstract::Group::ReturnType rtype;
  rtype = initialGuess.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utils.err() << "NOX::StatusTest::NormF::NormF - Unable to compute F"
		<< endl;
    throw "NOX Error";
  }

  initialTolerance_ = computeNorm(initialGuess);
  trueTolerance_ = specifiedTolerance_ / initialTolerance_;
}
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType
FSI::GenericNormF::checkStatus(const NOX::Solver::Generic& problem,
                                NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    normF_ = 0.0;
    status_ = NOX::StatusTest::Unevaluated;
  }
  else
  {
    normF_ = computeNorm( problem.getSolutionGroup() );
    if ((normF_ != -1) and (normF_ < trueTolerance_))
    {
      status_ = NOX::StatusTest::Converged;
    }
    else
    {
      status_ = NOX::StatusTest::Unconverged;
    }
  }

  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType FSI::GenericNormF::getStatus() const
{
  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& FSI::GenericNormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';

  stream << status_
         << name_ << "-Norm = " << NOX::Utils::sciformat(normF_,3)
         << " < " << NOX::Utils::sciformat(trueTolerance_, 3)
         << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';

  stream << setw(13) << " (";

  if (scaleType_ == Scaled)
    stream << "Length-Scaled";
  else
    stream << "Unscaled";

  stream << " ";

  if (normType_ == NOX::Abstract::Vector::TwoNorm)
    stream << "Two-Norm";
  else if (normType_ == NOX::Abstract::Vector::OneNorm)
    stream << "One-Norm";
  else if (normType_ == NOX::Abstract::Vector::MaxNorm)
    stream << "Max-Norm";

  stream << ", ";

  if (toleranceType_ == Absolute)
    stream << "Absolute Tolerance";
  else
    stream << "Relative Tolerance";

  stream << ")\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getNormF() const
{
  return normF_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getTrueTolerance() const
{
  return trueTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getSpecifiedTolerance() const
{
  return specifiedTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormF::getInitialTolerance() const
{
  return initialTolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PartialNormF::PartialNormF(std::string name,
                                const LINALG::MultiMapExtractor& extractor,
                                int blocknum,
                                double tolerance,
                                ScaleType stype)
  : GenericNormF(name,tolerance,stype),
    extractor_(extractor),
    blocknum_(blocknum)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::PartialNormF::computeNorm(const NOX::Abstract::Group& grp)
{
  if (!grp.isF())
    return -1.0;

  // extract the block epetra vector

  const NOX::Abstract::Vector& abstract_f = grp.getF();
  const NOX::Epetra::Vector& f = Teuchos::dyn_cast<const NOX::Epetra::Vector>(abstract_f);

  // extract the inner vector elements we are interested in

  Teuchos::RCP<Epetra_Vector> v = extractor_.ExtractVector(f.getEpetraVector(),blocknum_);

  return FSI::GenericNormF::computeNorm(*v);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GenericNormUpdate::GenericNormUpdate(double tol,
                                          NOX::Abstract::Vector::NormType ntype,
                                          ScaleType stype)
  : status_(NOX::StatusTest::Unevaluated),
    normType_(ntype),
    scaleType_(stype),
    tolerance_(tol),
    normUpdate_(0.0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GenericNormUpdate::GenericNormUpdate(double tol, ScaleType stype)
  : status_(NOX::StatusTest::Unevaluated),
    normType_(NOX::Abstract::Vector::TwoNorm),
    scaleType_(stype),
    tolerance_(tol),
    normUpdate_(0.0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::GenericNormUpdate::~GenericNormUpdate()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType FSI::GenericNormUpdate::checkStatus(const NOX::Solver::Generic& problem,
                                                                NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    status_ = NOX::StatusTest::Unevaluated;
    normUpdate_ = -1.0;
    return status_;
  }

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    status_ = NOX::StatusTest::Unconverged;
    normUpdate_ = -1.0;
    return status_;
  }

  // Check that F exists!
  if (!problem.getSolutionGroup().isF())
  {
    status_ = NOX::StatusTest::Unconverged;
    normUpdate_ = -1.0;
    return status_;
  }

  const NOX::Abstract::Vector& oldSoln = problem.getPreviousSolutionGroup().getX();
  const NOX::Abstract::Vector& curSoln = problem.getSolutionGroup().getX();

  if (Teuchos::is_null(updateVectorPtr_))
    updateVectorPtr_ = curSoln.clone();

  updateVectorPtr_->update(1.0, curSoln, -1.0, oldSoln, 0.0);

  computeNorm(Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>(updateVectorPtr_)->getEpetraVector());

  status_ = (normUpdate_ < tolerance_) ? NOX::StatusTest::Converged : NOX::StatusTest::Unconverged;
  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormUpdate::computeNorm(const Epetra_Vector& v)
{
  int n = (scaleType_ == Scaled) ? updateVectorPtr_->length() : 0;

  switch (normType_)
  {
  case NOX::Abstract::Vector::TwoNorm:
    normUpdate_ = updateVectorPtr_->norm();
    if (scaleType_ == Scaled)
      normUpdate_ /= sqrt(1.0 * n);
    break;

  default:
    normUpdate_ = updateVectorPtr_->norm(normType_);
    if (scaleType_ == Scaled)
      normUpdate_ /= n;
    break;

  }
  return normUpdate_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::StatusTest::StatusType FSI::GenericNormUpdate::getStatus() const
{
  return status_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& FSI::GenericNormUpdate::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status_
         << "Absolute Update-Norm = "
         << NOX::Utils::sciformat(normUpdate_, 3)
	 << " < "
         << NOX::Utils::sciformat(tolerance_, 3)
         << endl;
  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormUpdate::getNormUpdate() const
{
  return normUpdate_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::GenericNormUpdate::getTolerance() const
{
  return tolerance_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PartialNormUpdate::PartialNormUpdate(std::string name,
                                          const LINALG::MultiMapExtractor& extractor,
                                          int blocknum,
                                          double tolerance,
                                          ScaleType stype)
  : GenericNormUpdate(tolerance,stype),
    extractor_(extractor),
    blocknum_(blocknum)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::PartialNormUpdate::computeNorm(const Epetra_Vector& v)
{
  return FSI::GenericNormUpdate::computeNorm(*extractor_.ExtractVector(v,blocknum_));
}


#endif
