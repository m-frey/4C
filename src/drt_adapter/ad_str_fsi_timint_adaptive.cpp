/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint.cpp

\brief Structure field adapter for time step size adaptivity within monolithic FSI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ad_str_fsi_timint_adaptive.H"

#include "../drt_inpar/inpar_structure.H"

#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure/strtimada.H"
#include "../drt_structure/strtimint.H"
#include "../drt_structure/stru_aux.H"


/*======================================================================*/
/* constructor */
ADAPTER::StructureFSITimIntAda::StructureFSITimIntAda(
  Teuchos::RCP<STR::TimAda> sta,
  Teuchos::RCP<Structure> sti
)
: FSIStructureWrapper(sti),
  StructureTimIntAda(sta, sti)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& sada = sdyn.sublist("TIMEADAPTIVITY");

  // type of error norm
  errnorm_ = DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(sada,"LOCERRNORM");

  //----------------------------------------------------------------------------
  // Handling of Dirichlet BCs in error estimation
  //----------------------------------------------------------------------------
  // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(sti->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  numdbcdofs_ = sti->GetDBCMapExtractor()->CondMap()->NumGlobalElements();
  numdbcfsidofs_ = intersectionmap->NumGlobalElements();
  numdbcinnerdofs_ = numdbcdofs_ - numdbcfsidofs_;

}

/*----------------------------------------------------------------------------*/
/* Indicate norms of local discretization error */
void ADAPTER::StructureFSITimIntAda::IndicateErrorNorms(double& err,
                                                        double& errcond,
                                                        double& errother)
{
  // call functionality of adaptive structural time integrator
  StrAda()->EvaluateLocalErrorDis();

  // Indicate has to be done by the FSI algorithm since it depends on the interface
  IndicateErrors(err, errcond, errother);

  return;
}

/*----------------------------------------------------------------------------*/
/* Indicate local discretization error */
void ADAPTER::StructureFSITimIntAda::IndicateErrors(double& err,
                                                    double& errcond,
                                                    double& errother)
{
  // vector with local discretization error for each DOF
  Teuchos::RCP<Epetra_Vector> error = StrAda()->LocErrDis();

  // extract the condition part of the full error vector
  // (i.e. only interface displacement DOFs)
  Teuchos::RCP<Epetra_Vector> errorcond
    = Teuchos::rcp(new Epetra_Vector(*interface_->ExtractFSICondVector(error)));

  // in case of structure split: extract the other part of the full error vector
  // (i.e. only interior displacement DOFs)
  Teuchos::RCP<Epetra_Vector> errorother
    = Teuchos::rcp(new Epetra_Vector(*interface_->ExtractFSICondVector(error)));

  // calculate norms of different subsets of local discretization error vector
  err = STR::AUX::CalculateVectorNorm(errnorm_, error, numdbcdofs_);
  errcond = STR::AUX::CalculateVectorNorm(errnorm_, errorcond, numdbcfsidofs_);
  errother = STR::AUX::CalculateVectorNorm(errnorm_, errorother, numdbcinnerdofs_);

  return;
}

/*----------------------------------------------------------------------------*/
/* Do a single step with auxiliary time integration scheme */
void ADAPTER::StructureFSITimIntAda::TimeStepAuxiliar()
{
  StrAda()->IntegrateStepAuxiliar();
}

/*----------------------------------------------------------------------------*/
/* Calculate time step size suggestion */
const double ADAPTER::StructureFSITimIntAda::CalculateDt(const double norm)
{
  return StrAda()->CalculateDt(norm);
}

/*----------------------------------------------------------------------------*/
/* Get time step size of adaptive structural time integrator */
double ADAPTER::StructureFSITimIntAda::Dt() const
{
  return StrAda()->Dt();
}

/*----------------------------------------------------------------------------*/
/* Get target time \f$t_{n+1}\f$ of current time step */
double ADAPTER::StructureFSITimIntAda::GetTimeNew() const
{
  return StrAda()->GetTimeNew();
}

/*----------------------------------------------------------------------------*/
/* Set new time step size */
void ADAPTER::StructureFSITimIntAda::SetDt(const double dtnew)
{
  StrAda()->SetDt(dtnew);
}

/*----------------------------------------------------------------------------*/
/* Update step size */
void ADAPTER::StructureFSITimIntAda::UpdateStepSize(const double dtnew)
{
  StrAda()->UpdateStepSize(dtnew);
}

/*----------------------------------------------------------------------------*/
/*  Reset certain quantities to prepare repetition of current time step */
void ADAPTER::StructureFSITimIntAda::ResetStep()
{
  StrAda()->ResetStep();
}
