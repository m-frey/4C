/*----------------------------------------------------------------------*/
/*!
\file FSIStructureWrapper.cpp

\brief Structural adapter for FSI problems containing the interface
       and methods dependent on the interface

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/



#include "ad_str_fsiwrapper.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../drt_structure/stru_aux.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FSIStructureWrapper::FSIStructureWrapper(Teuchos::RCP<Structure> structure)
: StructureWrapper(structure)
{
  // set-up FSI interface
  interface_ = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface_->Setup(*Discretization(), *Discretization()->DofRowMap());

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  predictor_ = DRT::INPUT::IntegralValue<int>(fsidyn,"PREDICTOR");

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FSIStructureWrapper::UseBlockMatrix()
{
  StructureWrapper::UseBlockMatrix(interface_,interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FSIStructureWrapper::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> relax = interface_->InsertFSICondVector(iforce);
  SetForceInterface(relax);
  Teuchos::RCP<Epetra_Vector> idisi = SolveRelaxationLinear();

  // we are just interested in the incremental interface displacements
  return interface_->ExtractFSICondVector(idisi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FSIStructureWrapper::PredictInterfaceDispnp()
{
  // prestressing business
  double time = 0.0;
  double dt   = 0.0;
  double pstime = -1.0;
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = GetTime();
    dt = GetTimeStepSize();
    pstime = sdyn.get<double>("PRESTRESSTIME");
  }

  Teuchos::RCP<Epetra_Vector> idis;

  switch (predictor_)
  {
  case 1:
  {
    // d(n)
    // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
    if (pstype != INPAR::STR::prestress_none && time+dt <= pstime)
    {
      idis = Teuchos::rcp(new Epetra_Vector(*interface_->FSICondMap(),true));
    }
    else
    {
      idis  = interface_->ExtractFSICondVector(ExtractDispn());
    }
    break;
  }
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    // d(n)+dt*v(n)
    if (pstype != INPAR::STR::prestress_none && time+dt <= pstime)
      dserror("only constant interface predictor useful for prestressing");

    double dt = GetTimeStepSize();

    idis = interface_->ExtractFSICondVector(ExtractDispn());
    Teuchos::RCP<Epetra_Vector> ivel
      = interface_->ExtractFSICondVector(ExtractVeln());

    idis->Update(dt,* ivel, 1.0);
    break;
  }
  case 4:
  {
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    if (pstype != INPAR::STR::prestress_none && time+dt <= pstime)
      dserror("only constant interface predictor useful for prestressing");

    double dt = GetTimeStepSize();

    idis = interface_->ExtractFSICondVector(ExtractDispn());
    Teuchos::RCP<Epetra_Vector> ivel
      = interface_->ExtractFSICondVector(ExtractVeln());
    Teuchos::RCP<Epetra_Vector> iacc
      = interface_->ExtractFSICondVector(ExtractAccn());

    idis->Update(dt, *ivel, 0.5*dt*dt, *iacc, 1.0);
    break;
  }
  default:
    dserror("unknown interface displacement predictor '%s'",
             DRT::Problem::Instance()->FSIDynamicParams().get<std::string>("PREDICTOR").c_str());
    break;
  }

  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FSIStructureWrapper::ExtractInterfaceDispn()
{
  // prestressing business
  double time = 0.0;
  double pstime = -1.0;
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = GetTime();
    pstime = sdyn.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time <= pstime)
  {
    return Teuchos::rcp(new Epetra_Vector(*interface_->FSICondMap(),true));
  }
  else
  {
    return interface_->ExtractFSICondVector(ExtractDispn());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FSIStructureWrapper::ExtractInterfaceDispnp()
{
  // prestressing business
  double time = 0.0;
  double dt   = 0.0;
  double pstime = -1.0;
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  if (pstype != INPAR::STR::prestress_none)
  {
    time = GetTime();
    dt = GetTimeStepSize();
    pstime = sdyn.get<double>("PRESTRESSTIME");
  }

  if (pstype != INPAR::STR::prestress_none && time+dt <= pstime)
  {
    return Teuchos::rcp(new Epetra_Vector(*interface_->FSICondMap(),true));
  }
  else
  {
    return interface_->ExtractFSICondVector(ExtractDispnp());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// @name Apply interface forces

/// apply interface forces to structural solver
///
/// This prepares a new solve of the structural field within one time
/// step. The middle values are newly created.
///
/// \note This is not yet the most efficient implementation.
void ADAPTER::FSIStructureWrapper::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> fifc = LINALG::CreateVector(*Discretization()->DofRowMap(), true);

  interface_->AddFSICondVector(iforce, fifc);

  SetForceInterface(fifc);

  PreparePartitionStep();

  return;
}



