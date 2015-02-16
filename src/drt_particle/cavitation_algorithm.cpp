/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_adapter/adapter_particle.H"
#include "particle_node.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/cavitationfluid.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"

#include "../drt_inpar/inpar_cavitation.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  coupalgo_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::CouplingStrategyOverFields>(params,"COUPALGO")),
  void_frac_strategy_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::VoidFractionCalculation>(params,"VOID_FRACTION_CALC")),
  gauss_rule_per_dir_(params.get<int>("NUM_GP_VOID_FRACTION")),
  approxelecoordsinit_((bool)DRT::INPUT::IntegralValue<int>(params,"APPROX_ELECOORDS_INIT")),
  simplebubbleforce_((bool)DRT::INPUT::IntegralValue<int>(params,"SIMPLIFIED_BUBBLE_FORCES")),
  timestepsizeratio_(params.get<int>("TIME_STEP_SIZE_RATIO")),
  fluiddis_(Teuchos::null),
  fluid_(Teuchos::null),
  ele_volume_(Teuchos::null),
  fluidfracn_(Teuchos::null),
  fluidfracnp_(Teuchos::null)
{
  // setup fluid time integrator
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),DRT::Problem::Instance()->FluidDynamicParams(),"fluid",false));
  fluid_ = fluid->FluidField();

  // validate input file

  // check whether gravity acceleration for fluid and particles match
  if(gravity_acc_.Norm2() > 0.0)
  {
    std::vector<DRT::Condition*> condition;
    fluiddis_->GetCondition("VolumeNeumann", condition);

    if(condition.size() != 1)
      dserror("exactly one VOL NEUMANN boundary condition expected to represent body forces in fluid");
    const std::vector<int>*    onoff = condition[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = condition[0]->Get<std::vector<double> >("val"  );

    const int dim = 3;
    for(int i=0; i<dim; ++i)
    {
      if(gravity_acc_(i) != (*val)[i])
        dserror("body force for particles does not match body force for fluid");
      if(gravity_acc_(i) != 0.0 and (*onoff)[i] == 0)
        dserror("body force for %d. dof deactivated in VOL NEUMANN bc for fluid although body force acts on particles in this direction.", i);
    }

    // check whether an initial pressure field is set due to the gravity load
    const int startfuncno = DRT::Problem::Instance()->FluidDynamicParams().get<int>("STARTFUNCNO");
    if(startfuncno < 0)
      dserror("pressure field needs to be initialized due to gravity load");
  }

  if(!simplebubbleforce_)
  {
    // check for solver for L2 projection of velocity gradient
    if(DRT::Problem::Instance()->FluidDynamicParams().get<int>("VELGRAD_PROJ_SOLVER") < 0)
      dserror("no solver for L2 projection of velocity gradient specified: check VELGRAD_PROJ_SOLVER");
  }

  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull)
  {
    // check for correct time integration scheme of fluid
    if(fluid_->TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
      dserror("two way full coupled cavitation problem only works with TIMEINTEGR = Af_Gen_Alpha");

    // check for correct physical type of fluid
    if(fluid_->PhysicalType() != INPAR::FLUID::loma)
      dserror("two way full coupled cavitation problem only works with PHYSICAL_TYPE = Loma");

    // check fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
    if (id == -1)
      dserror("no cavitation fluid material specified");

    // check for solver for L2 projection
    if(DRT::Problem::Instance()->CavitationParams().get<int>("VOIDFRAC_PROJ_SOLVER") < 0)
      dserror("no solver for L2 projection of fluid fraction specified: check VOIDFRAC_PROJ_SOLVER");
  }
  else
  {
    // check for correct time integration scheme of fluid
    if(fluid_->TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
      dserror("momentum coupled or one-way coupled cavitation problem does not work with TIMEINTEGR = Af_Gen_Alpha");

    // check for correct physical type of fluid
    if(fluid_->PhysicalType() != INPAR::FLUID::incompressible)
      dserror("two way momentum and one way coupled cavitation problems only works with PHYSICAL_TYPE = Incompressible");

    // check fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id == -1)
      dserror("specify fluid material");
  }

  if(timestepsizeratio_ < 1)
    dserror("fluid time step must be a multiplicative greater or equal unity. Your choice: %d", timestepsizeratio_);

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
{
  const int nstep_particles = NStep()*timestepsizeratio_;
  // time loop
  while (NotFinished() || particles_->StepOld() < nstep_particles)
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // deal with particle inflow
    ParticleInflow();

    // transfer particles into their correct bins
    TransferParticles();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::InitCavitation()
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  Teuchos::RCP<Epetra_Map> fluidelecolmapold = Teuchos::rcp(new Epetra_Map(*fluiddis_->ElementColMap()));
  CreateBins(fluiddis_);

  // gather all fluid coleles in each bin for proper extended ghosting
  std::map<int, std::set<int> > fluideles;
  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcsBasedOnUnderlyingDiscret(fluiddis_, fluideles);

  // read out bubble inflow condition and set bubble inflows in corresponding bins
  // assumption: only row bins are available up to here
  BuildBubbleInflowCondition();

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,BINSTRATEGY::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, fluideles);

  // check whether extended ghosting includes standard ghosting
  for(int i=0; i<fluidelecolmapold->NumMyElements(); ++i)
    if( fluiddis_->ElementColMap()->MyGID(fluidelecolmapold->GID(i)) == false)
      dserror("extended ghosting does not include standard ghosting");

  // assign wall elements based on the fluid discretization to bins initially once
  SetupParticleWalls(fluiddis_);
  AssignWallElesToBins();

  // copy structural dynamic params list and adapt particle specific entries
  const Teuchos::ParameterList& cavitationdyn = DRT::Problem::Instance()->CavitationParams();

  // adapt time step properties for particles in case of independent time stepping
  Teuchos::RCP<Teuchos::ParameterList> adaptedcavitationdyn = Teuchos::rcp(new Teuchos::ParameterList(cavitationdyn));

  adaptedcavitationdyn->set<double>("TIMESTEP", cavitationdyn.get<double>("TIMESTEP") / (double)timestepsizeratio_);
  adaptedcavitationdyn->set<int>("NUMSTEP", timestepsizeratio_ * cavitationdyn.get<int>("NUMSTEP"));
  adaptedcavitationdyn->set<int>("RESTARTEVRY", timestepsizeratio_ * cavitationdyn.get<int>("RESTARTEVRY"));
  adaptedcavitationdyn->set<int>("UPRES", timestepsizeratio_ * cavitationdyn.get<int>("UPRES"));

  // create particle time integrator
  Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(*adaptedcavitationdyn, particledis_));
  particles_ = particles->ParticleField();

  // set cavitation algorithm into time integration
  particles_->SetParticleAlgorithm(Teuchos::rcp(this,false));
  particles_->Init();

  // compute volume of each fluid element and store it
  ele_volume_ = LINALG::CreateVector(*fluiddis_->ElementRowMap(), false);
  int numfluidele = fluiddis_->NumMyRowElements();
  for(int i=0; i<numfluidele; ++i)
  {
    DRT::Element* fluidele = fluiddis_->lRowElement(i);
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));
    double ev = GEO::ElementVolume( fluidele->Shape(), xyze );
    (*ele_volume_)[i] = ev;
  }

  // compute initial fluid fraction
  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
  {
    fluidfracnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), true);
    CalculateFluidFraction();
    // and copy values from n+1 to n
    // leading to an initial zero time derivative
    fluidfracn_ = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
    // set fluid fraction in fluid for computation
    SetFluidFraction();
  }
  else
  {
    // fluid fraction is assumed constant equal unity
    fluidfracnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), false);
    fluidfracnp_->PutScalar(1.0);
    fluidfracn_ = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
    Teuchos::RCP<Epetra_Vector> fluid_fraction = Teuchos::rcp(new Epetra_Vector(*fluiddis_->ElementRowMap()));
    fluid_fraction->PutScalar(1.0);
    // apply fluid fraction to fluid on element level for visualization purpose
    Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->SetFluidFraction(fluid_fraction);
  }

  // determine consistent initial acceleration for the particles
  CalculateAndApplyForcesToParticles();
  particles_->DetermineMassDampConsistAccel();

  // some output
  if (myrank_ == 0)
    IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  if(particles_->StepOld() % timestepsizeratio_ == 0)
  {
    IncrementTimeAndStep();
    PrintHeader();

    fluid_->PrepareTimeStep();
  }

  // apply dirichlet boundary conditions
  particles_->PrepareTimeStep();

  if(structure_ != Teuchos::null)
    structure_->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  if(particles_->StepOld() % timestepsizeratio_ == 0)
  {
    if(coupalgo_ == INPAR::CAVITATION::TwoWayFull || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
    {
      TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateFluidFraction");
      CalculateFluidFraction();
      SetFluidFraction();
    }
  }
  else
  {
    // some output
    if (myrank_ == 0)
      IO::cout << "particle substep no. " << (particles_->StepOld() % timestepsizeratio_)+1 << IO::endl;
  }

  // apply forces and solve particle time step
  PARTICLE::Algorithm::Integrate();

  if(particles_->StepOld() % timestepsizeratio_ == 0)
  {
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::IntegrateFluid");
    fluid_->Solve();
  }

  return;
}


/*----------------------------------------------------------------------*
 | apply fluid fraction to fluid field                     ghamm 04/14  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetFluidFraction()
{
  if(fluid_->PhysicalType() != INPAR::FLUID::loma && myrank_ == 0)
    IO::cout << "Info: Fluid fraction is calculated and can be visualized but it is not "
        "used for the actual calculation" << IO::endl;

  // compute intermediate values for time integration scheme
  Teuchos::RCP<Epetra_Vector> fluidfracaf = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
  Teuchos::RCP<Epetra_Vector> fluidfracam = Teuchos::rcp(new Epetra_Vector(*fluidfracn_));
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->GenAlphaIntermediateValues(fluidfracaf, fluidfracam);
  // compute time derivative of fluid fraction
  const double invdt = 1.0/Dt();
  Teuchos::RCP<Epetra_Vector> fluidfracdtam = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
  fluidfracdtam->Update(-invdt, *fluidfracn_, invdt);

  // set fluid fraction in fluid for computation
  fluid_->SetIterScalarFields(fluidfracaf, fluidfracam, fluidfracdtam, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | calculate fluid forces on particle and apply it         ghamm 01/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateAndApplyForcesToParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateAndApplyForcesToParticles");
  const int dim = 3;

  fluiddis_->ClearState();
  particledis_->ClearState();

  Teuchos::ParameterList p;
  if(!simplebubbleforce_)
  {
    // project velocity gradient of fluid to nodal level via L2 projection and store it in a ParameterList
    Teuchos::RCP<Epetra_MultiVector> projected_velgrad = FLD::UTILS::ComputeL2ProjectedVelGradient(fluiddis_,fluid_->Veln());
    fluiddis_->AddMultiVectorToParameterList(p, "velgradient", projected_velgrad);
  }

  // at the beginning of the coupling step: veln = velnp(previous step) and current velnp contains fluid predictor
  fluiddis_->SetState("veln",fluid_->Veln());
  fluiddis_->SetState("velnm",fluid_->Velnm());

  // state at n+1 contains already dbc values due to PrepareTimeStep(), otherwise n = n+1
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispnp();
  Teuchos::RCP<const Epetra_Vector> bubblevel = particles_->Velnp();
  Teuchos::RCP<const Epetra_Vector> bubbleacc = particles_->Accnp();
  Teuchos::RCP<const Epetra_Vector> bubbleradius = particles_->Radius();

  // vectors to be filled with forces,
  // note: global assemble is needed for fluidforces due to the case with large bins and small fluid eles
  Teuchos::RCP<Epetra_Vector> bubbleforces = LINALG::CreateVector(*particledis_->DofRowMap(),true);
  Teuchos::RCP<Epetra_FEVector> fluidforces = Teuchos::rcp(new Epetra_FEVector(*fluiddis_->DofRowMap()));

  // fluid density and dynamic viscosity
  double rho_l;
  double mu_l;
  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull)
  {
    // get cavitation material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
    if (id == -1)
      dserror("no cavitation fluid material specified");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::CavitationFluid* actmat = static_cast<const MAT::PAR::CavitationFluid*>(mat);
    rho_l = actmat->density_;
    mu_l = actmat->viscosity_;
  }
  else
  {
    // get fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id == -1)
      dserror("no fluid material specified");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
    rho_l = actmat->density_;
    mu_l = actmat->viscosity_;
  }

  // bubble density
  double rho_b = particles_->ParticleDensity();

  // check whether dbc are specified for particles at all
  Teuchos::RCP<const Epetra_Map> dbcmap = particles_->GetDBCMapExtractor()->CondMap();
  const bool haveparticledbc = dbcmap->NumGlobalElements();

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;
  Epetra_SerialDenseVector elevector4;
  Epetra_SerialDenseVector elevector5;

  // only row particles are evaluated
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    DRT::Node* currparticle = particledis_->lRowNode(i);
    // fill particle position
    static LINALG::Matrix<3,1> particleposition(false);
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim; ++d)
      particleposition(d) = (*bubblepos)[posx+d];


    //--------------------------------------------------------------------
    // 1st step: element coordinates of particle position in fluid element
    //--------------------------------------------------------------------

    // variables to store information about element in which the particle is located
    DRT::Element* targetfluidele = NULL;
    static LINALG::Matrix<3,1> elecoord(false);

    // find out in which fluid element the current particle is located
    if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    DRT::Element** fluidelesinbin = currbin->AssociatedFluidEles();
    int numfluidelesinbin = currbin->NumAssociatedFluidEle();

    std::set<int>::const_iterator eleiter;
    // search for underlying fluid element with fast search if desired
    for(int ele=0; ele<numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];
      const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord, approxelecoordsinit_);

      if(insideele == true)
      {
        targetfluidele = fluidele;
        // leave loop over all fluid eles in bin
        break;
      }
    }

    // repeat search for underlying fluid element with standard search in case nothing was found
    if(targetfluidele == NULL and approxelecoordsinit_ == true)
    {
      for(int ele=0; ele<numfluidelesinbin; ++ele)
      {
        DRT::Element* fluidele = fluidelesinbin[ele];
        const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

        //get coordinates of the particle position in parameter space of the element
        bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord, false);

        if(insideele == true)
        {
          targetfluidele = fluidele;
          // leave loop over all fluid eles in bin
          break;
        }
      }
    }


    //--------------------------------------------------------------------
    // 2nd step: forces on this bubble are calculated
    //--------------------------------------------------------------------

    if(targetfluidele == NULL)
    {
      std::cout << "INFO: currparticle with Id: " << currparticle->Id() << " and position: " << particleposition(0) << " "
          << particleposition(1) << " " << particleposition(2) << " " << " does not have an underlying fluid element -> no forces calculated" << std::endl;

      std::vector<double> tmpposition(dim);
      for(int d=0; d<dim; ++d)
        tmpposition[d] = particleposition(d);
      int bubbleBinId = ConvertPosToGid(tmpposition);
      std::cout << "particle is in binId: " << bubbleBinId << " while currbin->Id() is " << currbin->Id() <<
          " . The following number of fluid eles is in this bin:" << numfluidelesinbin << std::endl;

      // do not assemble forces for this bubble and continue with next bubble
      continue;
    }

    // get element location vector and ownerships
    std::vector<int> lm_f;
    std::vector<int> lmowner_f;
    std::vector<int> lmstride;
    targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(dim);
    elevector2.Size(dim);
    elevector3.Size(dim);

    // set action in order to calculate the velocity and material derivative of the velocity
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_mat_deriv_u_and_rot_u);
    params.set<double>("timestep",Dt());
    params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = rot of fluid vel)
    targetfluidele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    if(!simplebubbleforce_)
    {
      // Reshape element matrices and vectors and initialize to zero
      elevector4.Size(dim);
      elevector5.Size(dim);

      // set action in order to calculate the pressure gradient and divergence of the stress tensor
      Teuchos::ParameterList params_surfintegrals(p);
      params_surfintegrals.set<int>("action",FLD::calc_press_grad_and_div_eps);
      params_surfintegrals.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

      // call the element specific evaluate method (elevec4 = pressure gradient; elevec5 = viscous stress term)
      targetfluidele->Evaluate(params_surfintegrals,*fluiddis_,lm_f,elematrix1,elematrix2,elevector4,elevector5,elevector3);
    }

    // get bubble velocity and acceleration
    std::vector<double> v_bub(lm_b.size());
    DRT::UTILS::ExtractMyValues(*bubblevel,v_bub,lm_b);

    // get bubble radius
    const double r_bub = (*bubbleradius)[ particledis_->NodeRowMap()->LID(currparticle->Id()) ];

    // bubble Reynolds number
    static LINALG::Matrix<3,1> v_rel(false);
    for (int d=0; d<dim; ++d)
      v_rel(d) = elevector1[d] - v_bub[d];

    const double v_relabs = v_rel.Norm2();
    const double Re_b = 2.0 * r_bub * v_relabs * rho_l / mu_l;

    bool output = false;
    if(output)
    {
      std::cout << "id_bub: " << currparticle->Id() << " " << std::endl;
      std::cout << "pos_bub: " << particleposition(0) << " " << particleposition(1) << " " << particleposition(2) << " " << std::endl;
      std::cout << "radius_bub: " << r_bub << std::endl;
      std::cout << "v_bub: " << v_bub[0] << " " << v_bub[1] << " " << v_bub[2] << " " << std::endl;
      std::cout << "v_fl: " << elevector1[0] << " " << elevector1[1] << " " << elevector1[2] << " " << std::endl;
      std::cout << "v_rel: " << v_rel(0) << " " << v_rel(1) << " " << v_rel(2) << " " << std::endl;
      std::cout << "v_relabs: " << v_relabs << std::endl;
      std::cout << "bubble Reynolds number: " << Re_b << std::endl;
    }

    // variable to sum forces for the current bubble under observation
    static LINALG::Matrix<3,1> sumforces(false);
    /*------------------------------------------------------------------*/
    //// 2.1) drag force = 0.5 * c_d * rho_l * Pi * r_b^2 * |u-v| * (u-v) or
    //// Stokes law for very small Re: drag force = 6.0 * Pi * mu_l * r_b * (u-v)
    double coeff1 = 0.0;
    if(Re_b < 0.1)
    {
      coeff1 = 6.0 * M_PI * mu_l * r_bub;
    }
    else
    {
      double c_d = 0.0;
      if(Re_b < 1000.0)
        c_d = 24.0 * (1.0 + 0.15 * pow(Re_b,0.687)) / Re_b;
      else
        c_d = 0.44;

      coeff1 = 0.5 * c_d * rho_l * M_PI * r_bub * r_bub * v_relabs;
    }

    static LINALG::Matrix<3,1> dragforce(false);
    dragforce.Update(coeff1, v_rel);
    // assemble
    sumforces.Update(dragforce);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.2) lift force = c_l * rho_l * volume_b * (u-v) x rot_u   with rot_u = nabla x u
    const double c_l = 0.5;
    const double vol_b = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
    static LINALG::Matrix<3,1> rot_u(false);
    for (int d=0; d<dim; ++d)
      rot_u(d) = elevector3(d);

    LINALG::Matrix<3,1> liftforce = GEO::computeCrossProduct(v_rel, rot_u);

    const double coeff2 = c_l * rho_l * vol_b;
    liftforce.Scale(coeff2);
    // assemble
    sumforces.Update(1.0, liftforce, 1.0);
    // store forces for coupling to fluid
    static LINALG::Matrix<3,1> couplingforce(false);
    couplingforce.Update(sumforces);
    /*------------------------------------------------------------------*/

    // material fluid acceleration at bubble position
    static LINALG::Matrix<3,1> Du_Dt(false);
    for (int d=0; d<dim; ++d)
      Du_Dt(d) = elevector2[d];

    if(simplebubbleforce_)
    {
      /*------------------------------------------------------------------*/
      //// 2.3) gravity and buoyancy forces = volume_b * rho_bub * g - volume_b * rho_l * ( g - Du/Dt )

      static LINALG::Matrix<3,1> grav_buoy_force(false);
      grav_buoy_force.Update(rho_b, gravity_acc_);
      grav_buoy_force.Update(-rho_l, gravity_acc_, rho_l, Du_Dt, 1.0);
      grav_buoy_force.Scale(vol_b);
      // assemble
      sumforces.Update(1.0, grav_buoy_force, 1.0);
      /*------------------------------------------------------------------*/
    }
    else
    {
      /*------------------------------------------------------------------*/
      //// 2.3) gravity, pressure gradient and viscous stress term = volume_b * rho_bub * g + volume_b * ( -grad_p + dTau/dx )

      static LINALG::Matrix<3,1> grad_p(false);
      static LINALG::Matrix<3,1> visc_stress(false);
      for (int d=0; d<dim; ++d)
      {
        grad_p(d) = elevector4[d];
        visc_stress(d) = elevector5[d];
      }

      static LINALG::Matrix<3,1> grav_surface_force(false);
      grav_surface_force.Update(rho_b, gravity_acc_);
      grav_surface_force.Update(-1.0, grad_p, 2.0*mu_l, visc_stress, 1.0);
      grav_surface_force.Scale(vol_b);
      // assemble
      sumforces.Update(1.0, grav_surface_force, 1.0);
      /*------------------------------------------------------------------*/
    }

    /*------------------------------------------------------------------*/
    //// 2.4) virtual/added mass = c_VM * rho_l * volume_b * ( Du/Dt - Dv/Dt )
    //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * volume_b * Du/Dt }{ 1 + c_VM * rho_l / rho_b }
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) +         coeff3          * Du/Dt }{          coeff4          }
    const double c_VM = 0.5;
    const double coeff3 = c_VM * rho_l * vol_b;
    static LINALG::Matrix<3,1> bubbleforce(false);
    bool isdbc = false;
    if(haveparticledbc)
    {
      isdbc = dbcmap->MyGID(lm_b[0]);
      bool isdbc1 = dbcmap->MyGID(lm_b[1]);
      bool isdbc2 = dbcmap->MyGID(lm_b[2]);
      if(isdbc != isdbc1 || isdbc1 != isdbc2)
        dserror("one particle can only constrain all or none of the dofs with dbc");
    }

    if(!isdbc) // free flying bubble
    {
      /*------------------------------------------------------------------*/
      //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
      //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * volume_b * Du/Dt }{ 1 + c_VM * rho_l / rho_b }
      //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) +         coeff3          * Du/Dt }{          coeff4          }
      const double coeff4 = 1.0 + c_VM * rho_l / rho_b;
      const double invcoeff4 = 1.0 / coeff4;

      bubbleforce.Update(invcoeff4, sumforces, coeff3*invcoeff4, Du_Dt);
      /*------------------------------------------------------------------*/
    }
    else // dbc controlled bubble
    {
      // extract dbc accelerations
      static LINALG::Matrix<3,1> Dv_Dt(false);
      for (int d=0; d<dim; ++d)
        Dv_Dt(d) = (*bubbleacc)[posx+d];

      // compute bubble force
      const double m_b = vol_b * rho_b;
      bubbleforce.Update(m_b, Dv_Dt);
    }
    /*------------------------------------------------------------------*/

    //--------------------------------------------------------------------
    // 3rd step: assemble bubble/fluid forces
    //--------------------------------------------------------------------

    // assemble of bubble forces (note: row nodes evaluated)
    static Epetra_SerialDenseVector forcecurrbubble(3);
    for(int d=0; d<dim; ++d)
      forcecurrbubble[d] = bubbleforce(d);
    std::vector<int> lmowner_b(lm_b.size(), myrank_);
    LINALG::Assemble(*bubbleforces,forcecurrbubble,lm_b,lmowner_b);

    // coupling forces between fluid and particle only include certain forces
    switch(coupalgo_)
    {
    case INPAR::CAVITATION::TwoWayFull:
    case INPAR::CAVITATION::TwoWayMomentum:
    {
      // calculate added mass force
      static LINALG::Matrix<3,1> addedmassforce(false);
      const double m_b = vol_b * rho_b;
      addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);

      //// coupling force = -(dragforce + liftforce + addedmassforce); actio = reactio --> minus sign
      couplingforce.Update(-1.0, addedmassforce, -1.0);

      // assemble of fluid forces must be done globally because col entries in the fluid can occur
      // although only row particles are evaluated
      const int numnode = targetfluidele->NumNode();
      Epetra_SerialDenseVector funct(numnode);
      // get shape functions of the element; evaluated at the bubble position --> distribution
      DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),targetfluidele->Shape());
      // prepare assembly for fluid forces (pressure degrees do not have to be filled)

      const int numdofperfluidele = numnode*(dim+1);
      double val[numdofperfluidele];
      for(int iter=0; iter<numnode; ++iter)
      {
        for(int d=0; d<dim; ++d)
        {
          val[iter*(dim+1) + d] = funct[iter] * couplingforce(d);
        }
        // no contribution on pressure dof
        val[iter*(dim+1) + 3] = 0.0;
      }
      // do assembly of bubble forces on fluid
      int err = fluidforces->SumIntoGlobalValues(numdofperfluidele, &lm_f[0], &val[0]);
      if (err<0)
        dserror("summing into Epetra_FEVector failed");
      break;
    }
    case INPAR::CAVITATION::OneWay:
    case INPAR::CAVITATION::VoidFracOnly:
    {
      //// coupling force = 0
      couplingforce.PutScalar(0.0);
      break;
    }
    default:
      dserror("coupalgo not available");
      break;
    }

    //--------------------------------------------------------------------
    // 4th step: output
    //--------------------------------------------------------------------
    if(output)
    {
      // gravity
      double m_b = vol_b * rho_b;
      static LINALG::Matrix<3,1> gravityforce(false);
      gravityforce.Update(m_b, gravity_acc_);
      std::cout << "t: " << Time() << " gravity force         : " << gravityforce << std::endl;

      if(simplebubbleforce_)
      {
        static LINALG::Matrix<3,1> buoy_force(false);
        buoy_force.Update(-rho_l,gravity_acc_);
        buoy_force.Scale(vol_b);
        std::cout << "t: " << Time() << " buoy_force          : " << buoy_force << std::endl;

        static LINALG::Matrix<3,1> inertia_force(false);
        inertia_force.Update(rho_l,Du_Dt);
        inertia_force.Scale(vol_b);
        std::cout << "t: " << Time() << " inertia_force       : " << inertia_force << std::endl;
      }
      else
      {
        static LINALG::Matrix<3,1> grad_p(false);
        static LINALG::Matrix<3,1> visc_stress(false);
        for (int d=0; d<dim; ++d)
        {
          grad_p(d) = elevector4[d];
          visc_stress(d) = elevector5[d];
        }

        static LINALG::Matrix<3,1> pressgrad_force(false);
        pressgrad_force.Update(-vol_b,grad_p);
        std::cout << "t: " << Time() << " pressgrad force     : " << pressgrad_force << std::endl;

        static LINALG::Matrix<3,1> viscous_force(false);
        viscous_force.Update(2.0*mu_l*vol_b, visc_stress);
        std::cout << "t: " << Time() << " viscous force       : " << viscous_force << std::endl;
      }

      // added mass force
      static LINALG::Matrix<3,1> addedmassforce(false);
      addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);

      // drag, lift and added mass force
      std::cout << "t: " << Time() << " dragforce force     : " << dragforce << std::endl;
      std::cout << "t: " << Time() << " liftforce force     : " << liftforce << std::endl;
      std::cout << "t: " << Time() << " added mass force    : " << addedmassforce << std::endl;

      // sum over all bubble forces
      std::cout << "t: " << Time() << " particle force      : " << bubbleforce << std::endl;

      // fluid force
      std::cout << "t: " << Time() << " fluid force         : " << couplingforce << std::endl;
    }

  } // end iparticle

  //--------------------------------------------------------------------
  // 5th step: apply forces to bubbles and fluid field
  //--------------------------------------------------------------------
  particles_->SetForceInterface(bubbleforces);

  if(coupalgo_ == INPAR::CAVITATION::OneWay || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
    return; // leave here because nothing to add to fluid

  // call global assemble
  int err = fluidforces->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");

  switch(coupalgo_)
  {
  case INPAR::CAVITATION::TwoWayFull:
  {
    // divide nodal wise fluid forces by fluid fraction
    // due to the special choice of Euler-Lagrange coupling
    const int numnodes = fluid_->Discretization()->NumMyRowNodes();
    for(int i=0; i<numnodes; ++i)
    {
      // fluid fraction is stored in pressure dof
      const double invnodalfraction = 1.0 / (*fluidfracnp_)[i*4+3];

      for(int j=0; j<3; ++j)
      {
        (*(*fluidforces)(0))[i*4+j] *= invnodalfraction;
      }
    }
    // no break here
  }
  case INPAR::CAVITATION::TwoWayMomentum:
  {
    // apply forces to fluid
    fluid_->ApplyExternalForces(fluidforces);
    break;
  }
  default:
    dserror("this case is not yet implemented");
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | particles are inserted into domain                      ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ParticleInflow()
{
  // inflow only once in fluid time step -> special case independent time stepping
  if(particles_->StepOld() % timestepsizeratio_ != 0)
    return;

  std::map<int, std::list<Teuchos::RCP<BubbleSource> > >::const_iterator biniter;

  int timeforinflow = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    // all particles have the same inflow frequency --> it is enough to test one
    // assumption only valid in case of one condition or conditions with identical inflow frequency
    if(biniter->second.size() != 0)
    {
      double inflowtime = 1.0 / biniter->second.front()->inflow_freq_;
      if(Step() % ((int)(inflowtime/Dt())) == 0)
      {
        timeforinflow = 1;
        break;
      }
    }
  }

  int globaltimeforinflow = 0;
  particledis_->Comm().MaxAll(&timeforinflow, &globaltimeforinflow, 1);
  if(globaltimeforinflow == 0)
    return; // no inflow detected


  // initialize bubble id with largest bubble id in use + 1 (on each proc)
  int maxbubbleid = particledis_->NodeRowMap()->MaxAllGID()+1;

  // start filling particles
  int inflowcounter = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::set<Teuchos::RCP<DRT::Node>,BINSTRATEGY::Less> homelessparticles;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
      Teuchos::RCP<DRT::Node> newparticle = Teuchos::rcp(new PARTICLE::ParticleNode(newbubbleid, &inflow_position[0], myrank_));
      PlaceNodeCorrectly(newparticle, &inflow_position[0], homelessparticles);
      if(homelessparticles.size() != 0)
        dserror("New bubble could not be inserted on this proc! Bubble inflow broken.");
    }
    inflowcounter += (int)biniter->second.size();
  }

  std::cout << "Inflow of " << inflowcounter << " bubbles on proc " << myrank_ << std::endl;

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  // insert data for new bubbles into state vectors
  const Epetra_Map* dofrowmap = particledis_->DofRowMap();
  const Epetra_Map* noderowmap = particledis_->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> disn = particles_->WriteAccessDispnp();
  Teuchos::RCP<Epetra_Vector> veln = particles_->WriteAccessVelnp();
  Teuchos::RCP<Epetra_Vector> radiusn = particles_->WriteAccessRadius();
  Teuchos::RCP<Epetra_Vector> massn = particles_->WriteAccessMass();
  const double density = particles_->ParticleDensity();

  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::vector<double> inflow_vel = (*particleiter)->inflow_vel_;
      int inflow_vel_curve = (*particleiter)->inflow_vel_curve_;
      double inflow_radius = (*particleiter)->inflow_radius_;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;

      double curvefac = 1.0;
      // curves are numbered starting with 1 in the input file
      if(inflow_vel_curve > 0)
        curvefac = DRT::Problem::Instance()->Curve(inflow_vel_curve-1).f(Time());

      DRT::Node* currparticle = particledis_->gNode(newbubbleid);
      // get the first gid of a particle and convert it into a LID
      int lid = dofrowmap->LID(particledis_->Dof(currparticle, 0));
      for(int dim=0; dim<3; ++dim)
      {
        (*disn)[lid+dim] = inflow_position[dim];
        (*veln)[lid+dim] = inflow_vel[dim] * curvefac;
      }
      lid = noderowmap->LID(newbubbleid);
      (*radiusn)[lid] = inflow_radius;
      (*massn)[lid] = density * 4.0/3.0 * M_PI * pow(inflow_radius, 3.0);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update()
{
  // here is the transition from n+1 -> n
  PARTICLE::Algorithm::Update();

  if((particles_->StepOld()-1) % timestepsizeratio_ == 0)
  {
    fluid_->Update();

    // update fluid fraction
    fluidfracn_->Update(1.0,*fluidfracnp_,0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  // adapt time step properties for particles in case of independent time stepping
  PARTICLE::Algorithm::ReadRestart(timestepsizeratio_*restart);
  fluid_->ReadRestart(restart);

  // correct time and step in algorithm base
  SetTimeStep(fluid_->Time(),restart);

  // additionally read restart data for fluid fraction
  IO::DiscretizationReader reader(fluid_->Discretization(),restart);
  reader.ReadVector(fluidfracn_,"fluid_fraction");
  return;
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap, std::map<int, std::set<int> >& fluideles)
{
  //--------------------------------------------------------------------
  // 1st and 2nd step
  //--------------------------------------------------------------------

  PARTICLE::Algorithm::SetupGhosting(binrowmap);


  //--------------------------------------------------------------------
  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int> > extendedfluidghosting;
  {
    // do communication to gather all elements for extended ghosting
    const int numproc = fluiddis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many col bins it has
      int numbin = bincolmap_->NumMyElements();
      fluiddis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which col bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* bincolmap = bincolmap_->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = bincolmap[i];
      }
      fluiddis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; ++i)
      {
        sdata[binid[i]].insert(fluideles[binid[i]].begin(),fluideles[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, fluiddis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedfluidghosting = rdata;
      }
    }

    //reduce map of sets to one set and copy to a vector to create fluidcolmap
    std::set<int> redufluideleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedfluidghosting.begin(); iter!= extendedfluidghosting.end(); ++iter)
    {
      redufluideleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> fluidcolgids(redufluideleset.begin(),redufluideleset.end());
    Teuchos::RCP<Epetra_Map> fluidcolmap = Teuchos::rcp(new Epetra_Map(-1,(int)fluidcolgids.size(),&fluidcolgids[0],0,Comm()));

    fluiddis_->ExtendedGhosting(*fluidcolmap,true,true,true,false);

  }

  //--------------------------------------------------------------------
  // 4th step: assign fluid elements to bins
  //--------------------------------------------------------------------
  {
    for(std::map<int, std::set<int> >::const_iterator biniter=extendedfluidghosting.begin(); biniter!=extendedfluidghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(biniter->first));
      for(std::set<int>::const_iterator fluideleiter=biniter->second.begin(); fluideleiter!=biniter->second.end(); ++fluideleiter)
      {
        int fluideleid = *fluideleiter;
        currbin->AddAssociatedFluidEle(fluideleid, fluiddis_->gElement(fluideleid));
//          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" << fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
      }
    }
  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from fluid elements to bins            ghamm 07/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildElementToBinPointers(bool wallpointer)
{
  // first call base class to associate potential particle walls
  PARTICLE::Algorithm::BuildElementToBinPointers(wallpointer);

  // loop over column bins and fill fluid elements
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedFluidEle();
    const int* fluideleids = actbin->AssociatedFluidEleIds();
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for(int iele=0; iele<numfluidele; ++iele)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = fluiddis_->gElement(fluideleid);
    }
    actbin->BuildFluidElePointers(&fluidelements[0]);
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output()
{
  if((particles_->StepOld()-1) % timestepsizeratio_ == 0)
  {
    // call fluid output and add restart data for fluid fraction if necessary
    fluid_->Output();
    const int uprestart = DRT::Problem::Instance()->CavitationParams().get<int>("RESTARTEVRY");
    if( Step()%uprestart == 0 && uprestart != 0 )
      fluid_->DiscWriter()->WriteVector("fluid_fraction", fluidfracn_, IO::DiscretizationWriter::dofvector);
  }

  PARTICLE::Algorithm::Output();

  return;
}


/*----------------------------------------------------------------------*
 | get adjacent bins to corner, where ijk is in 1st octant ghamm 02/13  |
 *----------------------------------------------------------------------*/
std::vector<int> CAVITATION::Algorithm::AdjacentBinstoCorner(int* ijk)
{
  std::vector<int> adjbins;
  adjbins.reserve(8);

  // get all adjacent bins to the current corner, including the bin itself
  for(int i=-1;i<1;i++)
  {
    for(int j=-1;j<1;j++)
    {
      for(int k=-1;k<1;k++)
      {
        int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

        int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
        if(neighborgid != -1)
        {
          adjbins.push_back(neighborgid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return adjbins;
}


/*----------------------------------------------------------------------*
| setup of bubble inflow                                    ghamm 01/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildBubbleInflowCondition()
{
  // build inflow boundary condition
  std::vector<DRT::Condition*> conds;
  particledis_->GetCondition("ParticleInflow", conds);
  // unique bubbleinflow id over all inflow conditions
  int bubbleinflowid = 0;
  for (size_t i=0; i<conds.size(); ++i)
  {
    if(i>0)
      dserror("only taken care of one particle inflow condition so far. "
          "Remedy: bubble_source_ needs to be a vector of the current layout");
    /*
     * inflow condition --> bubble sources
     *
     *  example: num_per_dir = {4, 5, 1}
     *
     *       <-> (dist_x = (vertex2_x-vertex1_x)/(num_per_dir_x-1))
     *
     *   x  x  x  x<-------- vertex2
     *
     *   x  x  x  x
     *
     *   x  x  x  x   ^
     *                | (dist_y = (vertex2_y-vertex1_y)/(num_per_dir_y-1) )
     *   x  x  x  x   ^
     *
     *   x  x  x  x
     *   ^
     *   |
     * vertex1
     *
     */

    // extract data from inflow condition
    const std::vector<double>* vertex1 = conds[i]->Get<std::vector<double> >("vertex1");
    const std::vector<double>* vertex2 = conds[i]->Get<std::vector<double> >("vertex2");
    const std::vector<int>* num_per_dir = conds[i]->Get<std::vector<int> >("num_per_dir");
    const std::vector<double>* inflow_vel = conds[i]->Get<std::vector<double> >("inflow_vel");
    int inflow_vel_curve = conds[i]->GetInt("inflow_vel_curve");
    double inflow_freq = conds[i]->GetDouble("inflow_freq");

    // make sure that a particle material is defined in the dat-file
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if (id==-1)
      dserror("Could not find particle material");

    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
    double initial_radius = actmat->initialradius_;

    double inflowtime = 1.0 / inflow_freq;
    if(std::abs(inflowtime/Dt() - (int)(inflowtime/Dt())) > EPS9)
      dserror("1/inflow_freq with inflow_freq = %f cannot be divided by fluid time step %f", inflowtime, Dt());
/* MUST BE ADDED WHEN PARTICLE CONTACT IS CONSIDERED
    double inflow_vel_mag = sqrt((*inflow_vel)[0]*(*inflow_vel)[0] + (*inflow_vel)[1]*(*inflow_vel)[1] + (*inflow_vel)[2]*(*inflow_vel)[2]);
    if(initial_radius/inflow_vel_mag > inflowtime)
      dserror("Overlap for inflowing bubbles expected: initial_radius/inflow_vel_mag = %f s > inflow_freq = %f s", initial_radius/inflow_vel_mag, inflowtime);
*/
    // loop over all bubble inflow positions and fill them into bin when they are on this proc;
    // up to here, only row bins are available
    std::vector<double> source_pos(3);
    for(int z=0; z<(*num_per_dir)[2]; ++z)
    {
      double dist_z = ((*vertex2)[2] - (*vertex1)[2]) / ((((*num_per_dir)[2]-1)!=0) ? ((*num_per_dir)[2]-1) : 1);
      source_pos[2] = (*vertex1)[2] + z * dist_z;
      for(int y=0; y<(*num_per_dir)[1]; ++y)
      {
        double dist_y = ((*vertex2)[1] - (*vertex1)[1]) / ((((*num_per_dir)[1]-1)!=0) ? ((*num_per_dir)[1]-1) : 1);
        source_pos[1] = (*vertex1)[1] + y * dist_y;
        for(int x=0; x<(*num_per_dir)[0]; ++x)
        {
          double dist_x = ((*vertex2)[0] - (*vertex1)[0]) / ((((*num_per_dir)[0]-1)!=0) ? ((*num_per_dir)[0]-1) : 1);
          source_pos[0] = (*vertex1)[0] + x * dist_x;
          // check whether this source position is on this proc
          int binId = ConvertPosToGid(source_pos);
          bool found = particledis_->HaveGlobalElement(binId);
          if(found == true)
          {
            Teuchos::RCP<BubbleSource> bubbleinflow = Teuchos::rcp(new BubbleSource(
                                                                          bubbleinflowid,
                                                                          source_pos,
                                                                          *inflow_vel,
                                                                          inflow_vel_curve,
                                                                          initial_radius,
                                                                          inflow_freq));
            bubble_source_[binId].push_back(bubbleinflow);
#ifdef DEBUG
            if(particledis_->gElement(binId)->Owner() != myrank_)
              dserror("Only row bins should show up here. Either add additional if-case or move ghosting to a later point in time.");
#endif
          }
          bubbleinflowid++;
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | particle source                                         ghamm 02/13  |
 *----------------------------------------------------------------------*/
CAVITATION::BubbleSource::BubbleSource(
  int bubbleinflowid,
  std::vector<double> inflow_position,
  std::vector<double> inflow_vel,
  int inflow_vel_curve,
  double inflow_radius,
  double inflow_freq
  ) :
  inflowid_(bubbleinflowid),
  inflow_position_(inflow_position),
  inflow_vel_(inflow_vel),
  inflow_vel_curve_(inflow_vel_curve),
  inflow_radius_(inflow_radius),
  inflow_freq_(inflow_freq)
{
}
