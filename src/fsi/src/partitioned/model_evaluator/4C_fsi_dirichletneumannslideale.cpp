/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach
       with sliding ALE-structure interfaces



\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_fsi_dirichletneumannslideale.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_mortar_interface.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannSlideale::DirichletNeumannSlideale(const Epetra_Comm& comm)
    : DirichletNeumann(comm)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::setup()
{
  // call setup of base class
  FSI::DirichletNeumann::setup();

  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  set_kinematic_coupling(
      Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") == Inpar::FSI::CoupVarPart::disp);

  Inpar::FSI::SlideALEProj aletype = Core::UTILS::IntegralValue<Inpar::FSI::SlideALEProj>(
      Global::Problem::Instance()->FSIDynamicParams(), "SLIDEALEPROJ");

  slideale_ = Teuchos::rcp(new FSI::UTILS::SlideAleUtils(structure_field()->discretization(),
      MBFluidField()->discretization(), structure_fluid_coupling_mortar(), true, aletype));

  islave_ = Teuchos::rcp(new Epetra_Vector(*structure_fluid_coupling_mortar().SlaveDofMap(), true));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::Remeshing()
{
  // dispn and dispnp of structure, used for surface integral and velocity of the fluid in the
  // interface
  Teuchos::RCP<Epetra_Vector> idisptotal = structure_field()->extract_interface_dispnp();

  slideale_->Remeshing(*structure_field(), MBFluidField()->discretization(), idisptotal, islave_,
      structure_fluid_coupling_mortar(), Comm());

  // Evaluate solid/fluid Mortar coupling
  slideale_->EvaluateMortar(
      structure_field()->extract_interface_dispnp(), islave_, structure_fluid_coupling_mortar());
  // Evaluate solid/ale Mortar coupling
  slideale_->EvaluateFluidMortar(idisptotal, islave_);

  Teuchos::RCP<Epetra_Vector> unew =
      slideale_->InterpolateFluid(MBFluidField()->extract_interface_velnp());
  MBFluidField()->apply_interface_values(islave_, unew);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannSlideale::fluid_op(
    Teuchos::RCP<Epetra_Vector> idispcurr, const FillType fillFlag)
{
  FSI::Partitioned::fluid_op(idispcurr, fillFlag);

  if (fillFlag == User)
  {
    FOUR_C_THROW("not implemented");
    // SD relaxation calculation
    return fluid_to_struct(MBFluidField()->RelaxationSolve(struct_to_fluid(idispcurr), Dt()));
  }
  else
  {
    // normal fluid solve

    // the displacement -> velocity conversion at the interface
    const Teuchos::RCP<Epetra_Vector> ivel = interface_velocity(idispcurr);

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    // new Epetra_Vector for aledisp in interface
    Teuchos::RCP<Epetra_Vector> iale =
        Teuchos::rcp(new Epetra_Vector(*(structure_fluid_coupling_mortar().MasterDofMap()), true));

    Teuchos::RCP<Epetra_Vector> idispn = structure_field()->extract_interface_dispn();

    iale->Update(1.0, *idispcurr, 0.0);

    // iale reduced by old displacement dispn and instead added the real last displacements
    iale->Update(1.0, *ft_stemp_, -1.0, *idispn, 1.0);

    MBFluidField()->nonlinear_solve(struct_to_fluid(iale), struct_to_fluid(ivel));

    MBFluidField()->SetItemax(itemax);

    return fluid_to_struct(MBFluidField()->extract_interface_forces());
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannSlideale::struct_op(
    Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag)
{
  FSI::Partitioned::struct_op(iforce, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    return structure_field()->RelaxationSolve(iforce);
  }
  else
  {
    // normal structure solve
    structure_field()->apply_interface_forces(iforce);
    structure_field()->Solve();
    return structure_field()->extract_interface_dispnp();
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannSlideale::initial_guess()
{
  if (get_kinematic_coupling())
  {
    // real displacement of slave side at time step begin on master side --> for calcualtion of
    // FluidOp
    ft_stemp_ = fluid_to_struct(islave_);
    // predict displacement
    return structure_field()->predict_interface_dispnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (Core::UTILS::IntegralValue<int>(fsipart, "PREDICTOR") != 1)
    {
      FOUR_C_THROW(
          "unknown interface force predictor '%s'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return interface_force();
  }
}

FOUR_C_NAMESPACE_CLOSE
