/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all structural time integration strategies.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_TIMINT_BASE_HPP
#define BACI_STRUCTURE_NEW_TIMINT_BASE_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_structure_new.hpp"
#include "baci_io_every_iteration_writer.hpp"
#include "baci_structure_new_timint_basedataglobalstate.hpp"
#include "baci_structure_new_timint_basedataio.hpp"
#include "baci_structure_new_timint_basedatasdyn.hpp"

// forward declaration
class Epetra_Vector;
class Epetra_Map;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class BlockSparseMatrixBase;
}  // namespace CORE::LINALG
namespace ADAPTER
{
  class StructureTimeAda;
}
namespace STR
{
  class ModelEvaluator;
  class Dbc;
  class Integrator;
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TIMINT
  {
    /** \brief Abstract class for all time integration strategies
     *
     *  \author Michael Hiermeier */
    class Base : public ADAPTER::StructureNew, IO::EveryIterationWriterInterface
    {
      friend class ADAPTER::StructureTimeAda;

     public:
      /// constructor
      Base();


      /// initialize (all already existing) class variables
      virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
          const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate);

      /// setup of the new class variables
      void Setup() override;

      /// tests if there are more time steps to do
      [[nodiscard]] bool NotFinished() const override;

      /// reset everything (needed for biofilm simulations)
      void Reset() override;

      /** \brief reset step configuration after time step
       *
       *  This function is supposed to reset all variables which are directly related
       *  to the current new step n+1. To be more precise all variables ending with "Np"
       *  have to be reset. */
      void ResetStep() override;

      /// wrapper for things that should be done before PrepareTimeStep is called
      void PrePredict() override {}

      /// wrapper for things that should be done before solving the nonlinear iterations
      void PreSolve() override {}

      /// wrapper for things that should be done after convergence of Newton scheme
      void PostOutput() override {}

      /// things that should be done after the actual time loop is finished
      void PostTimeLoop() override;

      /// @name General access methods
      ///@{
      /// Access discretization (structure only)
      Teuchos::RCP<DRT::Discretization> Discretization() override;

      /// Access to pointer to DoF row map of the discretization (structure only)
      const Epetra_Map* DofRowMapView() override
      {
        CheckInit();
        return dataglobalstate_->DofRowMapView();
      }

      /// DoF map of structural vector of unknowns
      Teuchos::RCP<const Epetra_Map> DofRowMap() override
      {
        CheckInit();
        return dataglobalstate_->DofRowMap();
      }

      //! DoF map of vector of unknowns
      //! Alternative method capable of multiple DoF sets
      Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) override
      {
        CheckInit();
        return dataglobalstate_->DofRowMap(nds);
      }

      /// Access linear structural solver
      Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override
      {
        CheckInit();
        return datasdyn_->GetLinSolvers()[INPAR::STR::model_structure];
      }

      /// Return MapExtractor for Dirichlet boundary conditions
      Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override;
      [[nodiscard]] Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() const;

      //! Return locsys manager
      Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() override;

      //! Return the desired model evaluator (read-only)
      [[nodiscard]] const STR::MODELEVALUATOR::Generic& ModelEvaluator(
          INPAR::STR::ModelType mtype) const override;

      //! Return the desired model evaluator (read and write)
      STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype) override;

      ///@}

      /// Return domain map of the mass matrix (implicit and explicit)
      [[nodiscard]] const Epetra_Map& GetMassDomainMap() const override;

      /// @name Coupled problem routines
      /// @{
      /// wrapper for things that should be done before updating
      void PreUpdate() override {}

      /// Update routine for coupled problems with monolithic approach
      void Update() override;

      /// Update routine for coupled problems with monolithic approach with time adaptivity
      void Update(double endtime) override = 0;

      /// Update time and step counter
      virtual void UpdateStepTime();

      /// wrapper for things that should be done after solving the update
      void PostUpdate() override;
      /// @}

      /// @name Access global state from outside via adapter (needed for coupled problems)
      ///@{
      /// unknown displacements at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> DispNp() const override
      {
        CheckInit();
        return dataglobalstate_->GetDisNp();
      }

      /* \brief write access to displacements at \f$t^{n+1}\f$
       *
       * Calling this method makes only sense if state is supposed
       * to be manipulated. We must not forget to synchronize the
       * manipulated state with the NOX group.
       * Otherwise, the manipulations will be overwritten by NOX.
       * Therefore, we set the flag state_is_insync_with_noxgroup_
       * to false.
       * This will be checked:
       * See \ref CheckStateInSyncWithNOXGroup
       *
       * See also \ref ADAPTER::StructureNew::SetState
       */
      Teuchos::RCP<Epetra_Vector> WriteAccessDispNp() override
      {
        CheckInit();
        SetStateInSyncWithNOXGroup(false);
        return dataglobalstate_->GetDisNp();
      }

      /// known displacements at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> DispN() const override
      {
        CheckInit();
        return dataglobalstate_->GetDisN();
      }

      /// write access to displacements at \f$t^{n}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessDispN() override
      {
        CheckInit();
        return dataglobalstate_->GetDisN();
      }

      /// unknown velocities at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelNp() const override
      {
        CheckInit();
        return dataglobalstate_->GetVelNp();
      }

      /// write access to velocities at \f$t^{n+1}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessVelNp() override
      {
        CheckInit();
        return dataglobalstate_->GetVelNp();
      }

      /// unknown velocities at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelN() const override
      {
        CheckInit();
        return dataglobalstate_->GetVelN();
      }

      /// write access to velocities at \f$t^{n}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessVelN() override
      {
        CheckInit();
        return dataglobalstate_->GetVelN();
      }

      /// known velocities at \f$t_{n-1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelNm() const override
      {
        CheckInit();
        return dataglobalstate_->GetVelNm();
      }

      /// unknown accelerations at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> AccNp() const override
      {
        CheckInit();
        return dataglobalstate_->GetAccNp();
      }

      //! known accelerations at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> AccN() const override
      {
        CheckInit();
        return dataglobalstate_->GetAccN();
      }
      ///@}

      /// @name access and modify model evaluator stuff via adapter
      /// @{
      /// are there any algebraic constraints?
      bool HaveConstraint() override
      {
        CheckInitSetup();
        return datasdyn_->HaveModelType(INPAR::STR::model_lag_pen_constraint);
      }

      /// do we need a semi-smooth Newton-type plasticity algorithm
      virtual bool HaveSemiSmoothPlasticity()
      {
        CheckInitSetup();
        return datasdyn_->HaveEleTech(INPAR::STR::EleTech::plasticity);
      }

      /// FixMe get constraint manager defined in the structure
      Teuchos::RCP<CONSTRAINTS::ConstrManager> GetConstraintManager() override
      {
        dserror("Not yet implemented!");
        return Teuchos::null;
      }

      /// FixMe get contact/meshtying manager
      Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge() override
      {
        dserror("Not yet implemented!");
        return Teuchos::null;
      }

      /// do we have this model
      bool HaveModel(INPAR::STR::ModelType model) override
      {
        return datasdyn_->HaveModelType(model);
      }

      /// Add residual increment to Lagrange multipliers stored in Constraint manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) override
      {
        dserror("Not yet implemented!");
      }

      /// Add residual increment to pressures stored in Cardiovascular0D manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr) override
      {
        dserror("Not yet implemented!");
      }
      /// @}

      /// @name Time step helpers
      ///@{
      /// Return current time \f$t_{n}\f$ (derived)
      [[nodiscard]] double GetTimeN() const override
      {
        CheckInit();
        return dataglobalstate_->GetTimeN();
      }

      /// Sets the current time \f$t_{n}\f$ (derived)
      void SetTimeN(const double time_n) override
      {
        CheckInit();
        dataglobalstate_->GetTimeN() = time_n;
      }

      /// Return target time \f$t_{n+1}\f$ (derived)
      [[nodiscard]] double GetTimeNp() const override
      {
        CheckInit();
        return dataglobalstate_->GetTimeNp();
      }

      /// Sets the target time \f$t_{n+1}\f$ of this time step (derived)
      void SetTimeNp(const double time_np) override
      {
        CheckInit();
        dataglobalstate_->GetTimeNp() = time_np;
      }

      /// Get upper limit of time range of interest (derived)
      [[nodiscard]] double GetTimeEnd() const override
      {
        CheckInit();
        return datasdyn_->GetTimeMax();
      }

      /// Get upper limit of time range of interest (derived)
      void SetTimeEnd(double timemax) override
      {
        CheckInit();
        datasdyn_->GetTimeMax() = timemax;
      }

      /// Get time step size \f$\Delta t_n\f$
      [[nodiscard]] double GetDeltaTime() const override
      {
        CheckInit();
        return (*dataglobalstate_->GetDeltaTime())[0];
      }

      /// Set time step size \f$\Delta t_n\f$
      void SetDeltaTime(const double dt) override
      {
        CheckInit();
        (*dataglobalstate_->GetDeltaTime())[0] = dt;
      }

      /// Return time integration factor
      [[nodiscard]] double TimIntParam() const override;

      /// Return current step number \f$n\f$
      [[nodiscard]] int GetStepN() const override
      {
        CheckInit();
        return dataglobalstate_->GetStepN();
      }

      /// Sets the current step \f$n\f$
      void SetStepN(int step_n) override
      {
        CheckInit();
        dataglobalstate_->GetStepN() = step_n;
      }

      /// Return current step number $n+1$
      [[nodiscard]] int GetStepNp() const override
      {
        CheckInit();
        return dataglobalstate_->GetStepNp();
      }

      /// Sets the current step number \f$n+1\f$
      void SetStepNp(int step_np) override
      {
        CheckInitSetup();
        dataglobalstate_->GetStepNp() = step_np;
      }

      //! Get number of time steps
      [[nodiscard]] int GetStepEnd() const override
      {
        CheckInit();
        return datasdyn_->GetStepMax();
      }

      /// Sets number of time steps
      void SetStepEnd(int step_end) override
      {
        CheckInitSetup();
        datasdyn_->GetStepMax() = step_end;
      }

      //! Get divcont type
      [[nodiscard]] virtual enum INPAR::STR::DivContAct GetDivergenceAction() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivergenceAction();
      }

      //! Get number of times you want to halve your timestep in case nonlinear solver diverges
      [[nodiscard]] virtual int GetMaxDivConRefineLevel() const
      {
        CheckInitSetup();
        return datasdyn_->GetMaxDivConRefineLevel();
      }

      //! Get random factor for time step adaption
      [[nodiscard]] virtual double GetRandomTimeStepFactor() const
      {
        CheckInitSetup();
        return datasdyn_->GetRandomTimeStepFactor();
      }

      //! Set random factor for time step adaption
      virtual double SetRandomTimeStepFactor(double rand_tsfac)
      {
        CheckInitSetup();
        return datasdyn_->GetRandomTimeStepFactor() = rand_tsfac;
      }

      //! Get random factor for time step adaption
      [[nodiscard]] virtual int GetDivConRefineLevel() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivConRefineLevel();
      }

      //! Set random factor for time step adaption
      virtual int SetDivConRefineLevel(int divconrefinementlevel)
      {
        CheckInitSetup();
        return datasdyn_->GetDivConRefineLevel() = divconrefinementlevel;
      }

      //! Get random factor for time step adaption
      [[nodiscard]] virtual int GetDivConNumFineStep() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivConNumFineStep();
      }

      //! Set random factor for time step adaption
      virtual int SetDivConNumFineStep(int divconnumfinestep)
      {
        CheckInitSetup();
        return datasdyn_->GetDivConNumFineStep() = divconnumfinestep;
      }

      /// set evaluation action
      void SetActionType(const DRT::ELEMENTS::ActionType& action) override;

      // group id in nested parallelity
      [[nodiscard]] int GroupId() const;
      ///@}

      /// @name Structure with ale specific methods
      ///@{
      /// FixMe set/apply material displacements to structure field (structure with ale)
      void SetDispMatNp(Teuchos::RCP<Epetra_Vector> dispmatnp) override
      {
        dserror("Not supported at the moment!");
      }

      /// FixMe write access to material displacements (strutcure with ale) at \f$t^{n+1}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessDispMatNp() override
      {
        CheckInitSetup();
        dserror("Not yet supported!");
        return Teuchos::null;
      }
      ///@}


      /// Time adaptivity (derived pure virtual functionality)
      /// @{
      /// Resize MStep Object due to time adaptivity in FSI (derived)
      void ResizeMStepTimAda() override;

      /// @}

      /// Output writer related routines (file and screen output)
      /// @{
      /// Access output object
      Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() override
      {
        return DataIO().GetOutputPtr();
      }

      /// Calculate all output quantities depending on the constitutive model
      /// (and, hence, on a potential material history)
      void PrepareOutput(bool force_prepare_timestep) override;

      /// output results (implicit and explicit)
      virtual void Output() { Output(false); }
      void Output(bool forced_writerestart) override;

      /// output error norms
      virtual void OutputErrorNorms();

      /// Write Gmsh output for structural field
      void WriteGmshStrucOutputStep() override;

      /// FixMe Check if there are any elements with the micro material definition.
      /// Maybe the detection can be moved to the element loop in the ad_str_structure_new.cpp.
      /// There is already one.
      bool HaveMicroMat() override
      {
        dserror("Not yet considered!");
        return false;
      }

      /// create result test for encapsulated structure algorithm
      Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override;

      /** \brief Get data that is written during restart
       *
       *  This routine is only for simple structure problems!
       *  \date 06/13
       *  \author biehler */
      void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
          Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
          Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
          Teuchos::RCP<std::vector<char>> nodedata) override;

      /** Read restart values
       *
       * \param stepn (in): restart step at \f${n}\f$
       */
      void ReadRestart(const int stepn) override;

      /// Set restart values (deprecated)
      void SetRestart(int stepn,                        //!< restart step at \f${n}\f$
          double timen,                                 //!< restart time at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> disn,             //!< restart displacements at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> veln,             //!< restart velocities at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> accn,             //!< restart accelerations at \f$t_{n}\f$
          Teuchos::RCP<std::vector<char>> elementdata,  //!< restart element data
          Teuchos::RCP<std::vector<char>> nodedata      //!< restart element data
          ) override;
      /// @}

      /// Biofilm related stuff
      /// @{
      /// FixMe set structure displacement vector due to biofilm growth
      void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) override
      {
        dserror("Currently unsupported!");
      }
      /// @}

      /// @name Pure virtual adapter functions (have to be implemented in the derived classes)
      /// @{
      /// integrate the current step (implicit and explicit)
      virtual int IntegrateStep() = 0;
      /// right-hand-side of Newton's method (implicit only)
      Teuchos::RCP<const Epetra_Vector> RHS() override { return GetF(); };
      [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> GetF() const = 0;
      /// @}

     public:
      /// @name External accessors for the class variables
      ///@{
      /// Get the indicator if we are currently restarting the simulation
      [[nodiscard]] inline const bool& IsRestarting() const { return isrestarting_; }

      /// Get the indicator if we need to restart the initial state
      [[nodiscard]] inline bool IsRestartingInitialState() const
      {
        return datasdyn_->IsRestartingInitialState();
      }

      /// Get TimIntBase data for global state quantities (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataGlobalState> GetDataGlobalStatePtr() const
      {
        CheckInit();
        return dataglobalstate_;
      }

      /// Get TimIntBase data for global state quantities (read & write access)
      Teuchos::RCP<BaseDataGlobalState>& GetDataGlobalStatePtr()
      {
        CheckInit();
        return dataglobalstate_;
      }

      [[nodiscard]] const BaseDataGlobalState& GetDataGlobalState() const
      {
        CheckInit();
        return *dataglobalstate_;
      }

      /// Get TimIntBase data for io quantities (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataIO> GetDataIOPtr() const
      {
        CheckInit();
        return dataio_;
      }

      [[nodiscard]] const BaseDataIO& GetDataIO() const
      {
        CheckInit();
        return *dataio_;
      }

      /// Get TimIntBase data or struct dynamics quantitites (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataSDyn> GetDataSDynPtr() const
      {
        CheckInit();
        return datasdyn_;
      }

      [[nodiscard]] const BaseDataSDyn& GetDataSDyn() const
      {
        CheckInit();
        return *datasdyn_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (read access)
      [[nodiscard]] const STR::Dbc& GetDBC() const
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (write access)
      STR::Dbc& GetDBC()
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a pointer to the Dirichlet Boundary Condition handler (read access)
      [[nodiscard]] Teuchos::RCP<const STR::Dbc> GetDBCPtr() const
      {
        CheckInitSetup();
        return dbc_ptr_;
      }

      /// return the integrator (read-only)
      [[nodiscard]] const STR::Integrator& Integrator() const
      {
        CheckInitSetup();
        return *int_ptr_;
      }

      /// Get the global state
      const BaseDataGlobalState& DataGlobalState() const
      {
        CheckInit();
        return *dataglobalstate_;
      }

      /// Get internal TimIntBase data for structural dynamics quantities (read and write access)
      BaseDataSDyn& DataSDyn()
      {
        CheckInit();
        return *datasdyn_;
      }

      /// return a pointer to the Dirichlet Boundary Condition handler (read and write access)
      const Teuchos::RCP<STR::Dbc>& DBCPtr()
      {
        CheckInitSetup();
        return dbc_ptr_;
      }

      [[nodiscard]] bool HasFinalStateBeenWritten() const override;

      /// get the indicator state
      [[nodiscard]] inline const bool& IsInit() const { return isinit_; }

      /// get the indicator state
      [[nodiscard]] inline const bool& IsSetup() const { return issetup_; }

      //! @name Attribute access functions
      //@{

      //! Provide Name
      virtual enum INPAR::STR::DynamicType MethodName() const = 0;

      //! Provide title
      std::string MethodTitle() const { return INPAR::STR::DynamicTypeString(MethodName()); }

      //! Return true, if time integrator is implicit
      virtual bool IsImplicit() const = 0;

      //! Return true, if time integrator is explicit
      virtual bool IsExplicit() const = 0;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      virtual int MethodSteps() const = 0;

      //! Give order of accuracy
      int MethodOrderOfAccuracy() const
      {
        return std::min(MethodOrderOfAccuracyDis(), MethodOrderOfAccuracyVel());
      }

      //! Give local order of accuracy of displacement part
      virtual int MethodOrderOfAccuracyDis() const = 0;

      //! Give local order of accuracy of velocity part
      virtual int MethodOrderOfAccuracyVel() const = 0;

      //! Return linear error coefficient of displacements
      virtual double MethodLinErrCoeffDis() const = 0;

      //! Return linear error coefficient of velocities
      virtual double MethodLinErrCoeffVel() const = 0;

      //@}

      ///@}
     protected:
      /// Check if Init() and Setup() have been called, yet.
      inline void CheckInitSetup() const
      {
        dsassert(IsInit() and IsSetup(), "Call Init() and Setup() first!");
      }

      /// Check if Init() has been called
      inline void CheckInit() const { dsassert(IsInit(), "Call Init() first!"); }

      /// Get the global state
      BaseDataGlobalState& DataGlobalState()
      {
        CheckInit();
        return *dataglobalstate_;
      }

      /// Get the pointer to global state
      const Teuchos::RCP<BaseDataGlobalState>& DataGlobalStatePtr() const
      {
        CheckInit();
        return dataglobalstate_;
      }

      /// Get internal TimIntBase data for io quantities (read and write access)
      BaseDataIO& DataIO()
      {
        CheckInit();
        return *dataio_;
      }

      /// return a pointer to the input/output data container (read and write access)
      const Teuchos::RCP<BaseDataIO>& DataIOPtr()
      {
        CheckInit();
        return dataio_;
      }

      /// return a pointer to the structural dynamic data container (read and write access)
      const Teuchos::RCP<BaseDataSDyn>& DataSDynPtr()
      {
        CheckInit();
        return datasdyn_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (read and write access)
      STR::Dbc& DBC()
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a reference to the integrator (read and write access)
      STR::Integrator& Integrator()
      {
        CheckInitSetup();
        return *int_ptr_;
      }

      /// return a pointer to the integrator (read and write access)
      const Teuchos::RCP<STR::Integrator>& IntegratorPtr()
      {
        CheckInitSetup();
        return int_ptr_;
      }

      /** \brief Output to file
       *
       *  This routine always prints the last converged state, i.e.
       *  \f$D_{n}, V_{n}, A_{n}\f$.
       *
       *  \date 03/07
       *  \author mwgee (originally) */
      void OutputStep(bool forced_writerestart);

     private:
      /*! \brief Create a new input/output step in the output writer
       *
       * New step is created only once per time step. This is controlled by \c datawritten.
       * Do nothing if data has already been written in this time step.
       *
       * \param[in,out] Indicator whether data has already been written in this time step (true) or
       *                not (false)
       */
      void NewIOStep(bool& datawritten);

      /// output of the current state
      void OutputState();

      /** \brief output of the current state */
      void OutputState(IO::DiscretizationWriter& iowriter, bool write_owner) const;

      /** \brief output of the debug state */
      void OutputDebugState(IO::DiscretizationWriter& iowriter, bool write_owner) const override;

      /// output during runtime
      void RuntimeOutputState();

      /// output reaction forces
      void OutputReactionForces();

      /// output element volumes
      void OutputElementVolume(IO::DiscretizationWriter& iowriter) const;

      /// output stress and/or strain state
      void OutputStressStrain();

      /// output energy
      void OutputEnergy() const;

      /// output optional quantity
      void OutputOptionalQuantity();

      /// write restart information
      void OutputRestart(bool& datawritten);

      /// add restart information to output state
      void AddRestartToOutputState();

      /** \brief set the number of nonlinear iterations of the last time step
       *
       *  \pre UpdateStepTime() must be called beforehand, otherwise the wrong
       *  step-id is considered.
       *
       *  \author hiermeier \date 11/17 */
      void SetNumberOfNonlinearIterations();

      /** \brief decide which contributions to the total system energy shall be
       *         computed and written during simulation
       *
       *  \author grill */
      void SelectEnergyTypesToBeWritten();

      /** \brief initialize file stream for energy values and write all the
       *         column headers for the previously selected energy contributions
       *         to be written separately
       *
       *  \author grill */
      void InitializeEnergyFileStreamAndWriteHeaders();

     protected:
      /// flag indicating if Init() has been called
      bool isinit_;

      /// flag indicating if Setup() has been called
      bool issetup_;

      /// flag indicating that the simulation is currently restarting
      bool isrestarting_;

      /// flag indicating that displacement state was manipulated
      /// but NOX group has not been informed.
      bool state_is_insync_with_noxgroup_;

     protected:
      inline void SetStateInSyncWithNOXGroup(const bool insync)
      {
        state_is_insync_with_noxgroup_ = insync;
      }

      inline void ThrowIfStateNotInSyncWithNOXGroup() const
      {
        if (!state_is_insync_with_noxgroup_)
        {
          dserror(
              " state has been requested but the manipulated state has\n"
              "not been communicated to NOX.\n"
              "Manipulations made in the state vector will have no effect.\n"
              "Call SetState(x) to synchronize the states stored in the global\n"
              "state object and in the NOX group!");
        }
      }

     private:
      /// pointer to the different data containers
      Teuchos::RCP<BaseDataIO> dataio_;
      Teuchos::RCP<BaseDataSDyn> datasdyn_;
      Teuchos::RCP<BaseDataGlobalState> dataglobalstate_;

      /// pointer to the integrator (implicit or explicit)
      Teuchos::RCP<STR::Integrator> int_ptr_;

      /// pointer to the dirichlet boundary condition handler
      Teuchos::RCP<STR::Dbc> dbc_ptr_;
    };  // class Base
  }     // namespace TIMINT
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_TIMINT_BASE_H
