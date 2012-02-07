
#ifdef TRILINOS_PACKAGE

#include "fsi_nox_aitken.H"
#include "fsi_utils.H"

#include <NOX_Common.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>
#include <Teuchos_ParameterList.hpp>
#include <NOX_GlobalData.H>

#include "../drt_lib/drt_colors.H"

// debug output
#if 1

#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <NOX_Epetra_Vector.H>
#include "../drt_lib/standardtypes_cpp.H"
#ifdef CCADISCRET
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#else
extern struct _FILES  allfiles;
#endif
#endif

NOX::FSI::AitkenRelaxation::AitkenRelaxation(const Teuchos::RefCountPtr<NOX::Utils>& utils,
                                             Teuchos::ParameterList& params)
  : utils_(utils)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");
  nu_ = p.get("Start nu", 0.0);

  double maxstep = p.get("max step size", 0.0);
  if (maxstep > 0)
    nu_ = 1-maxstep;
}


NOX::FSI::AitkenRelaxation::~AitkenRelaxation()
{
}


bool NOX::FSI::AitkenRelaxation::reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
                                       Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Aitken");

  // do not reset the aitken factor
  //nu_ = p.get("Start nu", 0.0);

  // We might want to constrain the step size of the first relaxation
  // in a new time step.
  double maxstep = p.get("max step size", 0.0);
  if (maxstep > 0 && maxstep < 1-nu_)
    nu_ = 1-maxstep;

  if (!is_null(del_))
  {
    del_->init(1e20);
  }
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::AitkenRelaxation::compute(Abstract::Group& grp, double& step,
                                         const Abstract::Vector& dir,
                                         const Solver::Generic& s)
{
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n"
                  << "-- Aitken Line Search -- \n";
  }

  // debug output
#if 0
  {
    static int step;
    ostringstream filename;
#ifdef CCADISCRET
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
#else
    const std::string filebase = allfiles.outputfile_kenner;
#endif
    filename << filebase << "_" << step << ".aitken.QR";
    step += 1;

    ofstream out(filename.str().c_str());
    ::FSI::UTILS::MGS(s.getPreviousSolutionGroup(), dir, 10, out);
  }
#endif

  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  const NOX::Abstract::Vector& F = oldGrp.getF();

  // turn off switch
#if 1
  if (is_null(del_))
  {
    del_  = F.clone(ShapeCopy);
    del2_ = F.clone(ShapeCopy);
    del_->init(1.0e20);
    del2_->init(0.0);
  }

  del2_->update(1,*del_,1,F);
  del_ ->update(-1,F);

  const double top = del2_->innerProduct(*del_);
  const double den = del2_->innerProduct(*del2_);

  nu_ = nu_ + (nu_ - 1.)*top/den;
  step = 1. - nu_;
#else
  step = 1.;
#endif

  utils_->out() << "          RELAX = " YELLOW_LIGHT << setw(5) << step << END_COLOR "\n";

  grp.computeX(oldGrp, dir, step);

  // Calculate F anew here. This results in another FSI loop. However
  // the group will store the result, so it will be reused until the
  // group's x is changed again. We do not waste anything.
  grp.computeF();

  // is this reasonable at this point?
  double checkOrthogonality = fabs( grp.getF().innerProduct(dir) );

  if (utils_->isPrintType(Utils::InnerIteration))
  {
    utils_->out() << setw(3) << "1" << ":";
    utils_->out() << " step = " << utils_->sciformat(step);
    utils_->out() << " orth = " << utils_->sciformat(checkOrthogonality);
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }

  // write omega
#if 1
  double fnorm = grp.getF().norm();
  if (dynamic_cast<const NOX::Epetra::Vector&>(F).getEpetraVector().Comm().MyPID()==0)
  {
    static int count;
    static std::ofstream* out;
    if (out==NULL)
    {
#ifdef CCADISCRET
      std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
#else
      std::string s = allfiles.outputfile_kenner;
#endif
      s.append(".omega");
      out = new std::ofstream(s.c_str());
    }
    (*out) << count << " "
           << step << " "
           << fnorm
           << "\n";
    count += 1;
    out->flush();
  }
#endif

  // debug output
#if 0
  {
    static int step;
    ostringstream filename;
#ifdef CCADISCRET
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
#else
    const std::string filebase = allfiles.outputfile_kenner;
#endif
    filename << filebase << "_" << step << ".aitken.X";
    step += 1;

    ofstream out(filename.str().c_str());
    dir.print(out);
    //oldGrp.getX().print(out);
    //grp.getX().print(out);
  }
#endif

  return true;
}


#endif
