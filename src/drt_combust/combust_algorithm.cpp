/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

    detailed description in header file combust_algorithm.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_algorithm.H"
#include "combust_defines.H"
#include "combust_flamefront.H"
#include "combust_reinitializer.H"
#include "combust_utils.H"
#include "combust_fluidimplicitintegration.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_io/io_gmsh.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
  /* Der constructor sollte den gesamten Algorithmus initialisieren.
   * Hier müssen alle Variablen, die den Einzelfeldern übergeordnet sind, initialisiert werden.
   *
   * das heisst:
   * o G-Funktionsvektor (t und ig+1) auf initial value setzen
   * o Geschwindigkeitsvektor (t und iu+1) auf initial value setzen
   * o alle Zähler auf 0 setzen (step_(0), f_giter_(0), g_iter_(0), f_iter_(0))
   * o alle Normen und Grenzwerte auf 0 setzen
   *
   * Zusammenfassend kann an sagen, dass alles was in der Verbrennungsrechnung vor der Zeitschleife
   * passieren soll, hier passieren muss, weil die combust dyn gleich die Zeitschleife ruft.
   *
   * scalar transport velocity field has been initialized in ScaTraFluidCouplingAlgorithm()
  */
COMBUST::Algorithm::Algorithm(Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn,false),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  convtol_(combustdyn.get<double>("CONVTOL")),
  stepbeforereinit_(false),
  stepreinit_(false),
  phireinitn_(Teuchos::null),
//  fgvelnormL2_(?),
//  fggfuncnormL2_(?),
  combusttype_(Teuchos::getIntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  reinitaction_(Teuchos::getIntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
//  reinitaction_(combustdyn.sublist("COMBUSTION GFUNCTION").get<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION")),
  reinitinterval_(combustdyn.sublist("COMBUSTION GFUNCTION").get<int>("REINITINTERVAL")),
  reinitband_(Teuchos::getIntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITBAND")),
  reinitbandwidth_(combustdyn.sublist("COMBUSTION GFUNCTION").get<double>("REINITBANDWIDTH")),
  combustdyn_(combustdyn),
  interfacehandleNP_(Teuchos::null),
  interfacehandleN_(Teuchos::null),
  flamefront_(Teuchos::null)
  {

  if (Comm().MyPID()==0)
  {
    switch(combusttype_)
    {
    case INPAR::COMBUST::combusttype_premixedcombustion:
      std::cout << "COMBUST::Algorithm: this is a premixed combustion problem" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem with kinks in vel and jumps in pres" << std::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflowjump:
      std::cout << "COMBUST::Algorithm: this is a two-phase flow problem with jumps in vel and pres" << std::endl;
      break;
    default: dserror("unknown type of combustion problem");
    }
  }

  if (Teuchos::getIntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") == INPAR::FLUID::timeint_gen_alpha)
    dserror("Generalized Alpha time integration scheme not available for combustion");

  // get pointers to the discretizations from the time integration scheme of each field
  // remark: fluiddis cannot be of type "const Teuchos::RCP<const DRT::Dis...>", because parent
  // class. InterfaceHandle only accepts "const Teuchos::RCP<DRT::Dis...>"              henke 01/09
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  velnpip_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);
  velnpi_ = rcp(new Epetra_Vector(*fluiddis->DofRowMap()),true);

  phinpip_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
  phinpi_ = rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   *
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/
  // construct initial flame front
  flamefront_ = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  // construct interfacehandles using initial flame front
  interfacehandleNP_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));
  interfacehandleN_ = rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefront_));
  // get integration cells according to initial flame front
  interfacehandleNP_->UpdateInterfaceHandle();
  interfacehandleN_->UpdateInterfaceHandle();

//  std::cout << "No initialization of LevelSet" << std::endl;
  if (reinitaction_ != INPAR::COMBUST::reinitaction_none)
  {
    stepreinit_ = true;
    ReinitializeGfunc();
    if (Teuchos::getIntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn_,"TIMEINT") != INPAR::FLUID::timeint_stationary)
    {
      // reset phin vector in ScaTra time integration scheme to phinp vector
      *ScaTraField().Phin() = *ScaTraField().Phinp();
    }
    // pointer not needed any more
    stepreinit_ = false;
  }

  //------------------------
  // set initial fluid field
  //------------------------
  const INPAR::COMBUST::InitialField initfield = Teuchos::getIntegralValue<INPAR::COMBUST::InitialField>(
      combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD");
  const int initfuncno = combustdyn_.sublist("COMBUSTION FLUID").get<int>("INITFUNCNO");
  if (initfield == INPAR::COMBUST::initfield_flame_vortex_interaction)
  {
    // show flame front to fluid time integration scheme
    FluidField().ImportFlameFront(flamefront_);
  }
  FluidField().SetInitialFlowField(initfield, initfuncno);
  //FluidField().Output();
  if (initfield == INPAR::COMBUST::initfield_flame_vortex_interaction)
  {
    // delete fluid's memory of flame front; it should never have seen it in the first place!
    FluidField().ImportFlameFront(Teuchos::null);
  }
  // export interface information to the fluid time integration
  // remark: this is essential here, if DoFluidField() is not called in Timeloop() (e.g. for pure Scatra problems)
  FluidField().ImportInterface(interfacehandleNP_,interfacehandleN_);
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a dynamic combustion problem                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{
  // compute initial volume of minus domain
  const double volume_start = ComputeVolume();

  // get initial field by solving stationary problem first
  if(Teuchos::getIntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == true)
    SolveInitialStationaryProblem();

  // time loop
  while (NotFinished())
  {
    // prepare next time step; update field vectors
    PrepareTimeStep();

    // Fluid-G-function-Interaction loop
    while (NotConvergedFGI())
    {
      // prepare Fluid-G-function iteration
      PrepareFGIteration();

     // TODO @Ursula do G-function field first
//      std::cout << "/!\\ warning === fluid field is solved before G-function field" << std::endl;
//      DoFluidField();

      // TODO In the first iteration of the first time step the convection velocity for the
      //      G-function is zero, if a zero initial fluid field is used.
      //      -> Should the fluid be solved first?
      // solve linear G-function equation
      DoGfuncField();

      // update interface geometry
      UpdateInterface();

      // solve nonlinear Navier-Stokes system
      DoFluidField();

    } // Fluid-G-function-Interaction loop

    // write output to screen and files
    Output();
    //Remark (important for restart): the time level of phi (n+1, n or n-1) used to reconstruct the interface
    //                                conforming to the restart state of the fluid depends on the order
    //                                of Output() and UpdateTimeStep()

    // TODO sollte direkt nach DoGfuncField() gerufen werden
    if (stepreinit_)
    {
      // compute current volume of minus domain
      const double volume_current_before = ComputeVolume();
      // print mass conservation check on screen
      printMassConservationCheck(volume_start, volume_current_before);

      // reinitialize G-function
      ReinitializeGfunc();

      // compute current volume of minus domain
      const double volume_current_after = ComputeVolume();
      // print mass conservation check on screen
      printMassConservationCheck(volume_start, volume_current_after);
    }

    // update all field solvers
    UpdateTimeStep();

    if (!stepreinit_)
    {
      // compute current volume of minus domain
      const double volume_current = ComputeVolume();
      // print mass conservation check on screen
      printMassConservationCheck(volume_start, volume_current);
    }

  } // time loop

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start, volume_end);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary combustion problem                              henke 10/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("--------Stationary-Combustion-------  time step ----------------------------------------\n");
  }
  //-----------------------------
  // prepare stationary algorithm
  //-----------------------------
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // check if ScaTraField().initialvelset == true
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField().MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // compute initial volume of minus domain
  const double volume_start = ComputeVolume();

  //--------------------------------------
  // loop over fluid and G-function fields
  //--------------------------------------
  while (NotConvergedFGI())
  {
    // prepare Fluid-G-function iteration
    PrepareFGIteration();

    // solve nonlinear Navier-Stokes system
    DoFluidField();

    // solve (non)linear G-function equation
    std::cout << "/!\\ warning === G-function field not solved for stationary problems" << std::endl;
    //DoGfuncField();

    //reinitialize G-function
    //if (fgiter_ % reinitinterval_ == 0)
    //{
    //ReinitializeGfunc();
    //}

    // update field vectors
    UpdateInterface();

  } // fluid-G-function loop

  //-------
  // output
  //-------
  // remark: if Output() was already called at initial state, another Output() call will cause an
  //         error, because both times fields are written into the output control file at time and
  //         time step 0.
  //      -> the time and the time step have to be advanced even though this makes no physical sense
  //         for a stationary computation
  //IncrementTimeAndStep();
  //FluidField().PrepareTimeStep();
  //ScaTraField().PrepareTimeStep();

  // write output to screen and files (and Gmsh)
  Output();

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start, volume_end);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: reinitialize G-function                                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::ReinitializeGfunc()
{
  //  // turn on/off screen output for writing process of Gmsh postprocessing file
  //  const bool screen_out = true;
  //  // create Gmsh postprocessing file
  //  const std::string filename1 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_reinitialization_before", Step(), 701, screen_out, ScaTraField().Discretization()->Comm().MyPID());
  //  std::ofstream gmshfilecontent1(filename1.c_str());
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent1 << "View \" " << "Phinp \" {" << endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(ScaTraField().Discretization(),ScaTraField().Phinp(),gmshfilecontent1);
  //    gmshfilecontent1 << "};" << endl;
  //  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent1 << "View \" " << "Convective Velocity \" {" << endl;
  //    // draw vector field 'Convective Velocity' for every element
  //    IO::GMSH::VectorFieldNodeBasedToGmsh(ScaTraField().Discretization(),ScaTraField().ConVel(),gmshfilecontent1);
  //    gmshfilecontent1 << "};" << endl;
  //  }
  //  gmshfilecontent1.close();
  //  if (screen_out) std::cout << " done" << endl;

  //  // get my flame front (boundary integration cells)
  //  // remark: we generate a copy of the flame front here, which is not neccessary in the serial case
  //  std::map<int, GEO::BoundaryIntCells> myflamefront = interfacehandleNP_->GetElementalBoundaryIntCells();
  //
  //#ifdef PARALLEL
  //  // export flame front (boundary integration cells) to all processors
  //  flamefront_->ExportFlameFront(myflamefront);
  //#endif
  //
  //  // reinitialize what will later be 'phin'
  ////  if (stepbeforereinit_)
  ////  {
  ////    cout << "reinitializing phin_" << endl;
  ////    phireinitn_ = LINALG::CreateVector(*ScaTraField().Discretization()->DofRowMap(),true);
  ////    // copy phi vector of ScaTra time integration scheme
  ////    *phireinitn_ = *ScaTraField().Phinp();
  ////    // reinitialize G-function (level set) field
  ////    COMBUST::Reinitializer reinitializer(
  ////        combustdyn_,
  ////        ScaTraField(),
  ////        myflamefront,
  ////        phireinitn_);
  ////  }
  //
  //  // reinitialize what will later be 'phinp'
  //  if (stepreinit_)
  //  {
  //    // reinitialize G-function (level set) field
  //    COMBUST::Reinitializer reinitializer(
  //        combustdyn_,
  //        ScaTraField(),
  //        myflamefront,
  //        ScaTraField().Phinp());
  //
  //    // update flame front according to reinitialized G-function field
  //    flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp(),true);
  //
  //    // update interfacehandle (get integration cells) according to updated flame front
  //    interfacehandleNP_->UpdateInterfaceHandle();
  //  }

  // reinitialize what will later be 'phinp'
  if (stepreinit_)
  {
    // REMARK:
    // maybe there are some modified phi-values in phinp_ and so the current interfacehandle_ is adopted
    // to the modified values
    // we want to reinitialize the orignial phi-values
    // => UpdateFlameFront without modifyPhiVectors

    // update flame front according to reinitialized G-function field
    // the reinitilizer needs the original G-function field
    // ModifyPhiVector uses an alternative modification such that tetgen does not crash, so UpdateFlameFront is called
    // with the boolian true
    flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp(),true);

    // update interfacehandle (get integration cells) according to updated flame front
    interfacehandleNP_->UpdateInterfaceHandle();


    // get my flame front (boundary integration cells)
    // TODO we generate a copy of the flame front here, which is not neccessary in the serial case
    std::map<int, GEO::BoundaryIntCells> myflamefront = interfacehandleNP_->GetElementalBoundaryIntCells();

#ifdef PARALLEL
    // export flame front (boundary integration cells) to all processors
    flamefront_->ExportFlameFront(myflamefront);
#endif

    // reinitialize G-function (level set) field
    COMBUST::Reinitializer reinitializer(
        combustdyn_,
        ScaTraField(),
        myflamefront,
        ScaTraField().Phinp());

    // REMARK:
    // after the reinitialization we update the flamefront in the usual sense
    // this means that we modify phi-values if necessary -> default boolian modifyPhiVectors = true

    // update flame front according to reinitialized G-function field
    flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

    // update interfacehandle (get integration cells) according to updated flame front
    interfacehandleNP_->UpdateInterfaceHandle();
  }

  //  // create Gmsh postprocessing file
  //  const std::string filename2 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_reinitialization_after", Step(), 701, screen_out, ScaTraField().Discretization()->Comm().MyPID());
  //  std::ofstream gmshfilecontent2(filename2.c_str());
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent2 << "View \" " << "Phinp \" {" << endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(ScaTraField().Discretization(),ScaTraField().Phinp(),gmshfilecontent2);
  //    gmshfilecontent2 << "};" << endl;
  //  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent2 << "View \" " << "Convective Velocity \" {" << endl;
  //    // draw vector field 'Convective Velocity' for every element
  //    IO::GMSH::VectorFieldNodeBasedToGmsh(ScaTraField().Discretization(),ScaTraField().ConVel(),gmshfilecontent2);
  //    gmshfilecontent2 << "};" << endl;
  //  }
  //  gmshfilecontent2.close();
  //  if (screen_out) std::cout << " done" << endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: overwrite Navier-Stokes velocity                                        henke 08/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::OverwriteFluidVel()
{
  //----------------------------------------------------------------
  // impose velocity field by function e.g. for level set test cases
  // Navier-Stokes solution velocity field is overwritten
  //----------------------------------------------------------------
  if (Comm().MyPID()==0)
    std::cout << "\n--- overwriting Navier-Stokes solution ... " << std::flush;

  // get fluid (Navier-Stokes) velocity vector in standard FEM configuration (no XFEM dofs)
  const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();
  // velocity function number = 1 (FUNCT1)
  const int velfuncno = 1;

  // loop all nodes on the processor
  for(int lnodeid=0; lnodeid < FluidField().Discretization()->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    DRT::Node*  lnode = FluidField().Discretization()->lRowNode(lnodeid);
    // get standard dofset from fluid time integration
    vector<int> fluidnodedofs = (*(FluidField().DofSet())).Dof(lnode);
    // determine number of space dimensions (numdof - pressure dof)
    const int numdim = ((int) fluidnodedofs.size()) -1;
    if (numdim != 3) dserror("3 components expected for velocity");

    int err = 0;
    // overwrite velocity dofs only
    for(int index=0; index<numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int fgid = fluidnodedofs[index];
      int flid = convel->Map().LID(fgid);
      if (flid < 0) dserror("lid not found in map for given gid");

      // get value of corresponding velocity component
      double value = DRT::Problem::Instance()->Funct(velfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

      // insert velocity value into node-based vector
//        value = 0.0;
      err += convel->ReplaceMyValues(1, &value, &flid);
    }
    if (err != 0) dserror("error overwriting Navier-Stokes solution");
  }

  if (Comm().MyPID()==0)
    std::cout << "done" << std::endl;

  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | protected: compute flame velocity                                                  henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ComputeFlameVel(const Teuchos::RCP<Epetra_Vector>& convel,
    const Teuchos::RCP<const DRT::DofSet>& dofset)
{
  if((Teuchos::getIntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false) and
     (Teuchos::getIntegralValue<INPAR::COMBUST::InitialField>(combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD")
         == INPAR::COMBUST::initfield_zero_field))
  {
    cout << "/!\\ warning === Compute an initial stationary fluid solution to avoid a non-zero initial flame velocity" << endl;
  }
#if(0)
  if(Step()==1)
  {
    cout << "\n--- \t no scalar transport! replace convel by convel:= 0.0 everywere" << endl;
    // get a pointer to the fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();
    // get G-function value vector on fluid NodeColMap
    const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
    const Teuchos::RCP<Epetra_MultiVector> grad_phinp = flamefront_->GradPhi();

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {

      // get the processor local node
      DRT::Node*  lnode      = fluiddis->lRowNode(lnodeid);

      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);
      std::vector<int> lids(3);

      // extract velocity values (no pressure!) from global velocity vector
      for (int icomp=0; icomp<3; ++icomp)
      {
        lids[icomp] = convel->Map().LID(dofids[icomp]);
      }

      LINALG::Matrix<3,1> flvelabs(true);
      // add fluid velocity (Navier Stokes solution) and relative flame velocity
      for (int icomp=0; icomp<3; ++icomp)
      {
        flvelabs(icomp) = 0.0;
        convel->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
      }
    }
  }
#endif

#if(0)
    //-----------------------------------------------------------------------
    // use smoothed Phi-gradient (GradPhi) for calculation of the flame speed
    //-----------------------------------------------------------------------

    // get a pointer to the fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();

    // get G-function value vector on fluid NodeColMap
    const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
    const Teuchos::RCP<Epetra_MultiVector> grad_phinp = flamefront_->GradPhi();

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis->lRowNode(lnodeid);
      // get list of adjacent elements of this node
      DRT::Element** elelist = lnode->Elements();

      LINALG::Matrix<3,1> grad_phi(true);
      // get the global id for current node
      int GID = lnode->Id();

      // get local processor id according to global node id
      const int nodelid = (*grad_phinp).Map().LID(GID);

      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*grad_phinp).Comm().MyPID(),GID);

      const int numcol = (*grad_phinp).NumVectors();
      if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not identically to nsd");

      //----------------------------------------------------
      // compute normal vector at this node for this element
      // n = - grad phi / |grad phi|
      //----------------------------------------------------
      // evaluate gradient of G-function field at this node

      // loop over dimensions (= number of columns in multivector)
      for(int col=0; col< numcol; col++)
      {
        // get columns vector of multivector
        double* globalcolumn = (*grad_phinp)[col];


        // set smoothed gradient entry of phi into column of global multivector
        grad_phi(col) = globalcolumn[nodelid];
      }

      //------------------------
      // get material parameters
      //------------------------
      // get material from first (arbitrary!) element adjacent to this node
      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

      // density burnt domain
      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
      const double rhoplus = matplus->Density();

      // density unburnt domain
      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
      const double rhominus = matminus->Density();

      // laminar flame speed
      const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
      //---------------------------------------------
      // compute relative flame velocity at this node
      //---------------------------------------------
      // get phi value for this node
      const double gfuncval = (*phinp)[nodelid];

      double speedfac = 0.0;
      if (gfuncval >= 0.0){
        // burnt domain -> burnt material
        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
        speedfac = sl * rhominus/rhoplus;
      }
      else{
        // interface or unburnt domain -> unburnt material
        // flame speed factor = laminar flame speed
        speedfac = sl;
      }

      LINALG::Matrix<3,1> flvelrel(true);
      for (int icomp=0; icomp<3; ++icomp)
        flvelrel(icomp) = -1.0* speedfac * grad_phi(icomp)/grad_phi.Norm2();

      //-----------------------------------------------
      // compute (absolute) flame velocity at this node
      //-----------------------------------------------
      LINALG::Matrix<3,1> fluidvel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);
      std::vector<int> lids(3);

      // extract velocity values (no pressure!) from global velocity vector
      for (int icomp=0; icomp<3; ++icomp)
      {
        lids[icomp] = convel->Map().LID(dofids[icomp]);
        fluidvel(icomp) = (*convel)[lids[icomp]];
      }

      LINALG::Matrix<3,1> flvelabs(true);
      // add fluid velocity (Navier Stokes solution) and relative flame velocity
      for (int icomp=0; icomp<3; ++icomp)
      {
        flvelabs(icomp) = fluidvel(icomp) + flvelrel(icomp);
        convel->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
      }
    }
#endif
  // get a pointer to the fluid discretization
  const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField().Discretization();
  // get G-function value vector on fluid NodeColMap
  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();

#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
#endif

#ifdef COMBUST_GMSH_NORMALFIELD
  const std::string filestr = "flamefront_normal_field";
  const std::string name_in_gmsh = "Normal field";

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 500, true, fluiddis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
#endif

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis->lRowNode(lnodeid);
      // get list of adjacent elements of this node
      DRT::Element** elelist = lnode->Elements();

      //cout << "------------------------------------------------------------" << endl;
      //cout << "run for node: " << lnode->Id() << endl;
      //cout << "------------------------------------------------------------" << endl;

      //--------------------------------------------------------
      // compute "average"/"smoothed" normal vector at this node
      //--------------------------------------------------------
      // vector for average normal vector at this node
      LINALG::Matrix<3,1> avnvec(true);

      // loop over adjacent elements of this node
      for(int iele=0; iele<lnode->NumElement(); ++iele)
      {
        const DRT::Element* ele = elelist[iele];
        const int numnode = ele->NumNode();

        //cout << "------------------------------------------------------------" << endl;
        //cout << "run for element: " << ele->Id() << endl;
        //cout << "------------------------------------------------------------" << endl;

        // extract G-function values for nodes of this element
        Epetra_SerialDenseMatrix myphi(numnode,1);
        DRT::UTILS::ExtractMyNodeBasedValues(ele, myphi, *phinp);

        // get node coordinates of this element
        Epetra_SerialDenseMatrix xyze(3,numnode);
        GEO::fillInitialPositionArray<Epetra_SerialDenseMatrix>(ele, xyze);

        // TODO: function should be templated DISTYPE -> LINALG::Matrix<3,DISTYPE>
        Epetra_SerialDenseMatrix deriv(3,numnode);

#ifdef DEBUG
        bool nodefound = false;
#endif
        // find out which node in the element is my local node lnode
        for (int inode=0; inode<numnode; ++inode)
        {
          if (ele->NodeIds()[inode] == lnode->Id())
          {
            // get local (element) coordinates of this node
            LINALG::Matrix<3,1> coord = DRT::UTILS::getNodeCoordinates(inode,ele->Shape());
            // evaluate derivatives of shape functions at this node
            DRT::UTILS::shape_function_3D_deriv1(deriv,coord(0),coord(1),coord(2),ele->Shape());
#ifdef DEBUG
            nodefound = true;
#endif
          }
        }
#ifdef DEBUG
        if (nodefound==false)
          dserror("node was not found in list of elements");
#endif
        //----------------------------------------------------
        // compute normal vector at this node for this element
        // n = - grad phi / |grad phi|
        //----------------------------------------------------
        // evaluate gradient of G-function field at this node
        // remark: grad phi = sum (grad N_i * phi_i)

        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        Epetra_SerialDenseMatrix xjm(3,3);
        xjm.Multiply('N','T',1.0,deriv,xyze,0.0);

        // inverse of jacobian (xjm)
        LINALG::NonSymmetricInverse(xjm,3);
        // alternativ: Epetra_SerialDenseMatrix xji(3,3);
        // alternativ: xjm.Invert(xji);

        // compute global derivates
        Epetra_SerialDenseMatrix derxy(3,numnode);
        // derxy(i,j) = xji(i,k) * deriv(k,j)
        derxy.Multiply('N','N',1.0,xjm,deriv,0.0);
        // alternativ: xji.Multiply(false,deriv,derxy);

        Epetra_SerialDenseMatrix gradphi(3,1);
        derxy.Multiply(false,myphi,gradphi);
        double ngradphi = sqrt(gradphi(0,0)*gradphi(0,0)+gradphi(1,0)*gradphi(1,0)+gradphi(2,0)*gradphi(2,0));

        // normal vector at this node for this element
        LINALG::Matrix<3,1> nvec(true);
        if ((ngradphi < 1.0E-12) and (ngradphi > -1.0E-12)) // 'ngradphi' == 0.0
        {
          // length of normal is zero for this element -> level set must be constant within element (gradient is zero)
          std::cout << "/!\\ warning === no contribution to average normal vector from element " << ele->Id() << std::endl;
          // this element shall not contribute to the average normal vector
          continue;
        }
        else
        {
          for (int icomp=0; icomp<3; ++icomp) nvec(icomp) = -gradphi(icomp,0) / ngradphi;
        }
        //cout << "normal vector for element: " << ele->Id() << " at node: " << lnode->Id() << " is: " << nvec << endl;

        // add normal vector to linear combination (could also be weighted in different ways!)
        for (int icomp=0; icomp<3; ++icomp)
          avnvec(icomp) = avnvec(icomp) + nvec(icomp);
      }
      //---------------------------
      // compute unit normal vector
      //---------------------------
      // compute norm of average normal vector
      double avnorm = sqrt(avnvec(0)*avnvec(0)+avnvec(1)*avnvec(1)+avnvec(2)*avnvec(2));
      // divide vector by its norm to get unit normal vector
      if ((avnorm < 1.0E-12) and (avnorm > -1.0E-12)) // 'avnorm' == 0.0
      {
        // length of average normal is zero at this node -> node must be the tip of a
        // "regular level set cone" (all normals add up to zero normal vector)
        // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
        //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
        std::cout << "/!\\ warning === flame velocity at this node is only the convective velocity" << std::endl;
      }
      else
      {
        for (int icomp=0; icomp<3; ++icomp) avnvec(icomp) /= avnorm;
      }

#ifdef COMBUST_GMSH_NORMALFIELD
      LINALG::SerialDenseMatrix xyz(3,1);
      xyz(0,0) = lnode->X()[0];
      xyz(1,0) = lnode->X()[1];
      xyz(2,0) = lnode->X()[2];

      IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, avnvec, xyz, gmshfilecontent);
#endif

      //------------------------
      // get material parameters
      //------------------------
      // get material from first (arbitrary!) element adjacent to this node
      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

      // density burnt domain
      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
      const double rhoplus = matplus->Density();

      // density unburnt domain
      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
      const double rhominus = matminus->Density();

      // laminar flame speed
      const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
      //---------------------------------------------
      // compute relative flame velocity at this node
      //---------------------------------------------
      // get phi value for this node
      const int lid = phinp->Map().LID(lnode->Id());
      const double gfuncval = (*phinp)[lid];

      double speedfac = 0.0;
      if (gfuncval >= 0.0) // burnt domain -> burnt material
        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
        speedfac = sl * rhominus/rhoplus;
      else // interface or unburnt domain -> unburnt material
        // flame speed factor = laminar flame speed
        speedfac = sl;

      LINALG::Matrix<3,1> flvelrel(true);
      for (int icomp=0; icomp<3; ++icomp)
        flvelrel(icomp) = speedfac * avnvec(icomp);

      //-----------------------------------------------
      // compute (absolute) flame velocity at this node
      //-----------------------------------------------
      LINALG::Matrix<3,1> fluidvel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);
      std::vector<int> lids(3);

      // extract velocity values (no pressure!) from global velocity vector
      for (int icomp=0; icomp<3; ++icomp)
      {
        lids[icomp] = convel->Map().LID(dofids[icomp]);
        fluidvel(icomp) = (*convel)[lids[icomp]];
      }

      LINALG::Matrix<3,1> flvelabs(true);
      // add fluid velocity (Navier Stokes solution) and relative flame velocity
      for (int icomp=0; icomp<3; ++icomp)
      {
        flvelabs(icomp) = fluidvel(icomp) + flvelrel(icomp);
        convel->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
      }
      //cout << "------------------------------------------------------------" << endl;
      //cout << "run for node: " << lnode->Id() << endl;
      //cout << "------------------------------------------------------------" << endl;
      //cout << "convection velocity: " << flvelabs << endl;
    }

#ifdef COMBUST_GMSH_NORMALFIELD
    gmshfilecontent << "};\n";
  }
  gmshfilecontent.close();
  if (true) std::cout << " done" << endl;
#endif

return convel;
}

/*------------------------------------------------------------------------------------------------*
 | protected: FGI iteration converged?                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
  // if (fgiter <= fgitermax and ComputeGfuncNorm() < maxepsg and ComputeFluidNorm() < maxepsf)
  //return (fgiter_ < fgitermax_ and true);

  bool notconverged = true;
  
  // TODO fix this mess; implement a posteriori estimator
  //if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion or
  //    combusttype_ == INPAR::COMBUST::combusttype_twophaseflow)

  if (combusttype_ == INPAR::COMBUST::combusttype_twophaseflow)
  {
      /* at the moment only the convergence of the g-function field is checked
       * to check the convergence of the fluid field, uncomment the corresponding lines
       */
      if (fgiter_ == 0)
      {
         // G-function field and fluid field haven't been solved, yet

         // store old solution vectors
//        if (velnpip_->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
//          dserror("vectors must have the same length 1");
//        velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
         phinpip_->Update(1.0,*(ScaTraField().Phinp()),0.0);

         if (fgiter_ == fgitermax_)
         notconverged = false;

      }
      else if (fgiter_ > 0)
      {
//          double velnormL2 = 1.0;
          double gfuncnormL2 = 1.0;

          // store new solution vectors and compute L2-norm
//        velnpi_->Update(1.0,*velnpip_,0.0);
          phinpi_->Update(1.0,*phinpip_,0.0);
//        velnpip_->Update(1.0,*(FluidField().ExtractInterfaceVeln()),0.0);
//        velnpip_->Norm2(&velnormL2);
          phinpip_->Update(1.0,*(ScaTraField().Phinp()),0.0);
          phinpip_->Norm2(&gfuncnormL2);

//          if (velnormL2 < 1e-5) velnormL2 = 1.0;
          if (gfuncnormL2 < 1e-5) gfuncnormL2 = 1.0;

//        fgvelnormL2_ = 0.0;
          fggfuncnormL2_ = 0.0;

          // compute increment and L2-norm of increment
//        Teuchos::RCP<Epetra_Vector> incvel = rcp(new Epetra_Vector(velnpip_->Map()),true);
//        if (incvel->MyLength() != FluidField().ExtractInterfaceVeln()->MyLength())
//          dserror("vectors must have the same length 2");
//        incvel->Update(1.0,*velnpip_,-1.0,*velnpi_,0.0);
//        incvel->Norm2(&fgvelnormL2_);
          Teuchos::RCP<Epetra_Vector> incgfunc = rcp(new Epetra_Vector(*ScaTraField().Discretization()->DofRowMap()),true);
          incgfunc->Update(1.0,*phinpip_,-1.0,*phinpi_,0.0);
          incgfunc->Norm2(&fggfuncnormL2_);

          if (Comm().MyPID()==0)
         {
             printf("\n|+------------------------ FGI ------------------------+|");
             printf("\n|iter/itermax|----tol-[Norm]--|-fluid inc--|-g-func inc-|");
             printf("\n|   %2d/%2d    | %10.3E[L2] | ---------- | %10.3E |",fgiter_,fgitermax_,convtol_,fggfuncnormL2_/gfuncnormL2);
             printf("\n|+-----------------------------------------------------+|\n");
          }

          if (fggfuncnormL2_/gfuncnormL2 <= convtol_) //((fgvelnormL2_/velnormL2 <= convtol_) and (fggfuncnormL2_/gfuncnormL2 <= convtol_))
          {
             notconverged = false;
          }
          else
          {
             if (fgiter_ == fgitermax_)
             {
                 notconverged = false;
                 if (Comm().MyPID()==0)
                 {
                    printf("|+---------------- not converged ----------------------+|");
                    printf("\n|+-----------------------------------------------------+|\n");
                 }
             }
          }
      }
  }
  else // INPAR::COMBUST::combusttype_premixedcombustion
  {
    if (fgiter_ < fgitermax_)
    {
    }
    else
    {
      // added by me 21/10/09
      notconverged = false;
    }
  }

  return notconverged;
}

/*------------------------------------------------------------------------------------------------*
 | protected: do a stationary first time step for combustion algorithm               schott 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveInitialStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("==============================================================================================\n");
    printf("----------------Stationary timestep prepares instationary algorithm---------------------------\n");
    printf("==============================================================================================\n");
  }
  //-----------------------------
  // prepare stationary algorithm
  //-----------------------------
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // check if ScaTraField().initialvelset == true
  /* remark: initial velocity field has been transfered to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // modify time and timestep for stationary timestep
  SetTimeStep(0.0,0); // algorithm timestep

  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }

  FluidField().PrepareTimeStep();
  
  // compute initial volume of minus domain
  const double volume_start = ComputeVolume();

  //-------------------------------------
  // solve nonlinear Navier-Stokes system
  //-------------------------------------
  DoFluidField();

  // update field vectors
  UpdateInterface();

  //-------
  // output
  //-------
  // remark: if Output() was already called at initial state, another Output() call will cause an
  //         error, because both times fields are written into the output control file at time and
  //         time step 0.
  //      -> the time and the time step have to be advanced even though this makes no physical sense
  //         for a stationary computation
  //IncrementTimeAndStep();
  //FluidField().PrepareTimeStep();
  //ScaTraField().PrepareTimeStep();

  // write output to screen and files (and Gmsh)
  Output();

  // compute final volume of minus domain
  const double volume_end = ComputeVolume();
  // print mass conservation check on screen
  printMassConservationCheck(volume_start, volume_end);
  
  return;
}



/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  fgiter_ = 0;
  // fgnormgfunc = large value, determined in Input file
  fgvelnormL2_ = 1.0;
  // fgnormfluid = large value
  fggfuncnormL2_ = 1.0;

  // TODO @Florian clarify if this parameter is still needed
//  stepbeforereinit_ = false;
//  if (Step()>0 and Step() % reinitinterval_ == 0) stepbeforereinit_ = true;
  stepreinit_ = false;
  if (Step() % reinitinterval_ == 0) stepreinit_ = true; //Step()>1 and

  if (Comm().MyPID()==0)
  {
    //cout<<"---------------------------------------  time step  ------------------------------------------\n";
    printf("----------------------Combustion-------  time step %2d ----------------------------------------\n",Step());
    printf("TIME: %11.4E/%11.4E  DT = %11.4E STEP = %4d/%4d \n",Time(),MaxTime(),Dt(),Step(),NStep());
  }

  FluidField().PrepareTimeStep();
  // TODO @Martin: Kommentar einfuegen wofuer und warum
  interfacehandleN_->UpdateInterfaceHandle();

  // prepare time step
  // remark: initial velocity field has been transferred to scalar transport field in constructor of
  //         ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such
  //         as the one-step-theta scheme, are thus initialized correctly.
  ScaTraField().PrepareTimeStep();

  // synchronicity check between combust algorithm and base algorithms
  if (FluidField().Time() != Time())
    dserror("Time in Fluid time integration differs from time in combustion algorithm");
  if (ScaTraField().Time() != Time())
    dserror("Time in ScaTra time integration  differs from time in combustion algorithm");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
  fgiter_ += 1;
  if (Comm().MyPID()==0)
  {
    //cout<<"\n---------------------------------------  FGI loop  -------------------------------------------\n";
    printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
  }
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  if (Comm().MyPID()==0)
  {
    std:: cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------" << std::endl;
  }

  // TODO do we want an input parameter for this?
  if (true) // INPAR::XFEM::timeintegration_semilagrangian
  {
    // show flame front to fluid time integration scheme
    FluidField().ImportFlameFront(flamefront_);
  }
  // export interface information to the fluid time integration
  FluidField().ImportInterface(interfacehandleNP_,interfacehandleN_);
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  // solve nonlinear Navier-Stokes equations
  FluidField().NonlinearSolve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a G-function time integration step                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncField()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n---------------------------------------  G-FUNCTION SOLVER  ----------------------------------\n";
  }
  // assign the fluid velocity field to the G-function as convective velocity field
  switch(combusttype_)
  {
  case INPAR::COMBUST::combusttype_twophaseflow:
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    // for two-phase flow, the fluid velocity field is continuous; it can be directly transferred to
    // the scalar transport field

    ScaTraField().SetVelocityField(
      //OverwriteFluidVel(),
      FluidField().ExtractInterfaceVeln(),
      Teuchos::null,
      FluidField().DofSet(),
      FluidField().Discretization()
    );

    // Transfer history vector only for subgrid-velocity
    //ScaTraField().SetVelocityField(
    //    FluidField().ExtractInterfaceVeln(),
    //    FluidField().Hist(),
    //    FluidField().DofSet(),
    //    FluidField().Discretization()
    //);
    break;
  }
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    // for combustion, the velocity field is discontinuous; the relative flame velocity is added

    // extract convection velocity from fluid solution
    const Teuchos::RCP<Epetra_Vector> convel = FluidField().ExtractInterfaceVeln();

#if 0
    //std::cout << "convective velocity is transferred to ScaTra" << std::endl;
    Epetra_Vector convelcopy = *convel;
    //    *copyconvel = *convel;

    const Teuchos::RCP<Epetra_Vector> xfemvel = ComputeFlameVel(convel,FluidField().DofSet());
    if((convelcopy).MyLength() != (*xfemvel).MyLength())
      dserror("length is not the same!");

    const int dim = (convelcopy).MyLength();
    int counter = 0;

    for(int idof=0; idof < dim ;++idof)
    {
      if ((convelcopy)[idof] == (*xfemvel)[idof])
        counter++;
    }
    // number of identical dofs in convection velocity and flame velocity vector
    cout << "number of identical velocity components: " << counter << endl;
#endif

    ScaTraField().SetVelocityField(
        //OverwriteFluidVel(),
        //FluidField().ExtractInterfaceVeln(),
        ComputeFlameVel(convel,FluidField().DofSet()),
        Teuchos::null,
        FluidField().DofSet(),
        FluidField().Discretization()
    );
    break;
  }
  default:
    dserror("unknown type of combustion problem");
  }

  // TODO @Martin: Besprechung 18.8.2010, fuer mehr als eine FGI Iteration muss hier vielleicht ein
  //               Vektor umgesetzt werden ( flamefront_->SetOldPhiVector(); )

  //solve convection-diffusion equation
  ScaTraField().Solve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateInterface()
{
  // update flame front according to evolved G-function field
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  // update interfacehandle (get integration cells) according to updated flame front
  interfacehandleNP_->UpdateInterfaceHandle();

  // update the Fluid and the FGI vector at the end of the FGI loop
  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField().Update();

  if (stepreinit_)
  {
    ScaTraField().UpdateReinit();
  }
  else
  {
     ScaTraField().Update();
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  //------------------------------------------------------------------------------------------------
  // this hack is necessary for the visualization of disconituities in Gmsh             henke 10/09
  //------------------------------------------------------------------------------------------------
  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  FluidField().Output();
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  // causes error in DEBUG mode (trueresidual_ is null)
  //FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


void COMBUST::Algorithm::printMassConservationCheck(const double volume_start, const double volume_end)
{
  if (Comm().MyPID() == 0)
  {
    // compute mass loss
    if (volume_start == 0.0)
      dserror(" there is no 'minus domain'! -> division by zero checking mass conservation");
    double massloss = -(volume_start - volume_end) / volume_start *100;
    // 'isnan' seems to work not reliably; error occurs in line above
    if (std::isnan(massloss))
      dserror("NaN detected in mass conservation check");

    std::cout << "---------------------------------------" << endl;
    std::cout << "           mass conservation           " << endl;
    std::cout << " initial mass: " << volume_start << endl;
    std::cout << " final mass:   " << volume_end   << endl;
    std::cout << " mass loss:    " << massloss << "%" << endl;
    std::cout << "---------------------------------------" << endl;
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | compute volume on all processors                                                   henke 02/10 |
 *------------------------------------------------------------------------------------------------*/
double COMBUST::Algorithm::ComputeVolume()
{
  // compute negative volume of discretization on this processor
  double myvolume = interfacehandleNP_->ComputeVolumeMinus();

  double sumvolume = 0.0;

  // sum volumes on all processors
  // remark: ifndef PARALLEL sumvolume = myvolume
  Comm().SumAll(&myvolume,&sumvolume,1);

  return sumvolume;
}


/* -------------------------------------------------------------------------------*
 * Restart (fluid is solved before g-func)                               rasthofer|
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Restart(int step)
{
  if (Comm().MyPID()==0)
    std::cout << "Restart of combustion problem" << std::endl;

  // restart of scalar transport (G-function) field
  ScaTraField().ReadRestart(step);

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  //--------------------------
  // write output to Gmsh file
  //--------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_after_restart", Step(), 701, true, gfuncdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phin \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phin(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinm \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinm(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Convective Velocity \" {" << endl;
    // draw vector field 'Convective Velocity' for every element
    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField().ConVel(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  gmshfilecontent.close();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  Teuchos::RCP<COMBUST::FlameFront> flamefrontOld = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));

  // export phi n-1 vector from scatra dof row map to fluid node column map
  const Teuchos::RCP<Epetra_Vector> phinrow = rcp(new Epetra_Vector(*fluiddis->NodeRowMap()));
  if (phinrow->MyLength() != ScaTraField().Phin()->MyLength())
    dserror("vectors phinrow and phin must have the same length");
  *phinrow = *ScaTraField().Phin();
  const Teuchos::RCP<Epetra_Vector> phincol = rcp(new Epetra_Vector(*fluiddis->NodeColMap()));
  LINALG::Export(*phinrow,*phincol);

  // reconstruct old flame front
  //flamefrontOld->ProcessFlameFront(combustdyn_,ScaTraField().Phin());
  flamefrontOld->ProcessFlameFront(combustdyn_,phincol);

  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOldNP =
      rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOldN =
        rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  interfacehandleOldNP->UpdateInterfaceHandle();
  interfacehandleOldN->UpdateInterfaceHandle();
  FluidField().ImportInterface(interfacehandleOldNP,interfacehandleOldN);

  // restart of fluid field
  FluidField().ReadRestart(step);

  // reset interface for restart
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());

  interfacehandleNP_->UpdateInterfaceHandle();
  interfacehandleN_->UpdateInterfaceHandle();
  //-------------------
  // write fluid output
  //-------------------
  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  FluidField().Output();
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  SetTimeStep(FluidField().Time(),step);

  UpdateTimeStep();

  return;
}

/* -------------------------------------------------------------------------------*
 * Restart (g-func is solved before fluid)                               rasthofer|
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::RestartNew(int step)
{
  if (Comm().MyPID()==0)
    std::cout << "Restart of combustion problem" << std::endl;

  // restart of scalar transport (G-function) field
  ScaTraField().ReadRestart(step);

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField().Discretization();

  //--------------------------
  // write output to Gmsh file
  //--------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_after_restart", Step(), 701, true, gfuncdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phin \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phin(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinm \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField().Phinm(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Convective Velocity \" {" << endl;
    // draw vector field 'Convective Velocity' for every element
    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField().ConVel(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  gmshfilecontent.close();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  Teuchos::RCP<COMBUST::FlameFront> flamefrontOld = rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));

  // export phi n-1 vector from scatra dof row map to fluid node column map
  const Teuchos::RCP<Epetra_Vector> phinprow = rcp(new Epetra_Vector(*fluiddis->NodeRowMap()));
  if (phinprow->MyLength() != ScaTraField().Phinp()->MyLength())
    dserror("vectors phinrow and phin must have the same length");
  *phinprow = *ScaTraField().Phinp();
  const Teuchos::RCP<Epetra_Vector> phinpcol = rcp(new Epetra_Vector(*fluiddis->NodeColMap()));
  LINALG::Export(*phinprow,*phinpcol);

  // reconstruct old flame front
  flamefrontOld->ProcessFlameFront(combustdyn_,phinpcol);

  // build interfacehandle using old flame front
  // TODO @Martin Test + Kommentar
  // remark: interfacehandleN = interfacehandleNP, weil noch aeltere Information nicht vorhanden

  // TODO remove old code when new code tested
  //Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOld =
  //  rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  //interfacehandleOld->UpdateInterfaceHandle();
  //FluidField().ImportInterface(interfacehandleOld);

  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOldNP =
      rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandleOldN =
        rcp(new COMBUST::InterfaceHandleCombust(fluiddis,gfuncdis,flamefrontOld));
  interfacehandleOldNP->UpdateInterfaceHandle();
  interfacehandleOldN->UpdateInterfaceHandle();
  FluidField().ImportInterface(interfacehandleOldNP,interfacehandleOldN);

  // restart of fluid field
  FluidField().ReadRestart(step);

  //-------------------
  // write fluid output
  //-------------------
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField().Phin(), ScaTraField().Phinp());
  interfacehandleNP_->UpdateInterfaceHandle();
  interfacehandleN_->UpdateInterfaceHandle();
  // show flame front to fluid time integration scheme
  FluidField().ImportFlameFront(flamefront_);
  FluidField().Output();
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  FluidField().ImportFlameFront(Teuchos::null);

  SetTimeStep(FluidField().Time(),step);

  UpdateTimeStep();

  return;
}

#endif // #ifdef CCADISCRET
