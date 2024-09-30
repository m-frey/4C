/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled monolithic solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_MONOLITHIC_TWOWAY_HPP
#define FOUR_C_POROMULTIPHASE_MONOLITHIC_TWOWAY_HPP

#include "4C_config.hpp"

#include "4C_inpar_poromultiphase.hpp"
#include "4C_poromultiphase_monolithic.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class Solver;
  class Equilibration;
  enum class EquilibrationMethod;
}  // namespace Core::LinAlg

namespace Core::LinearSolver
{
  enum class SolverType;
}

namespace Core::Conditions
{
  class LocsysManager;
}

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithicTwoWay : public PoroMultiPhaseMonolithic
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithicTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    /// setup
    void setup_system() override;

    /// time step of coupled problem
    void time_step() override;

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> extractor() const override
    {
      return blockrowdofmap_;
    }

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    void evaluate(Teuchos::RCP<const Core::LinAlg::Vector> sx,
        Teuchos::RCP<const Core::LinAlg::Vector> fx, const bool firstcall) override;

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables) public for access from monolithic scatra problem
    void update_fields_after_convergence(Teuchos::RCP<const Core::LinAlg::Vector>& sx,
        Teuchos::RCP<const Core::LinAlg::Vector>& fx) override;

    // access to monolithic rhs vector
    Teuchos::RCP<const Core::LinAlg::Vector> rhs() const override { return rhs_; }

    // access to monolithic block system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() const override
    {
      return systemmatrix_;
    }

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> combined_dbc_map() const override { return combinedDBCMap_; };

   protected:
    //! Newton output to screen
    virtual void newton_output();

    //! Newton error check after loop
    virtual void newton_error_check();

    //! build the combined dirichletbcmap
    virtual void build_combined_dbc_map();

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map();

    virtual void setup_rhs();

    virtual Teuchos::RCP<Core::LinAlg::Vector> setup_structure_partof_rhs();

    //! build block vector from field vectors, e.g. rhs, increment vector
    void setup_vector(Core::LinAlg::Vector& f,        //!< vector of length of all dofs
        Teuchos::RCP<const Core::LinAlg::Vector> sv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Core::LinAlg::Vector> fv   //!< vector containing only fluid dofs
    );

    //! extract the field vectors from a given composed vector x.
    /*!
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations,
     and 1D artery pressure)
     */
    virtual void extract_field_vectors(Teuchos::RCP<const Core::LinAlg::Vector> x,
        Teuchos::RCP<const Core::LinAlg::Vector>& sx, Teuchos::RCP<const Core::LinAlg::Vector>& fx);

    //! extract only the structure and fluid field vectors from a given composed vector x.
    /*!
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
     */
    void extract_structure_and_fluid_vectors(Teuchos::RCP<const Core::LinAlg::Vector> x,
        Teuchos::RCP<const Core::LinAlg::Vector>& sx, Teuchos::RCP<const Core::LinAlg::Vector>& fx);

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix() { setup_system_matrix(*systemmatrix_); }

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat);

    /// setup composed system matrix from field solvers
    virtual void setup_maps();

    // Setup solver for monolithic system
    bool setup_solver() override;

    //! build the block null spaces
    void build_block_null_spaces(Teuchos::RCP<Core::LinAlg::Solver>& solver) override;

    //! Evaluate mechanical-fluid system matrix
    virtual void apply_str_coupl_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
    );

    //! Evaluate fluid-mechanical system matrix
    virtual void apply_fluid_coupl_matrix(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
    );

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void evaluate(Teuchos::RCP<const Core::LinAlg::Vector> iterinc);

    //! return structure fluid coupling sparse matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> struct_fluid_coupling_matrix();

    //! return fluid structure coupling sparse matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fluid_struct_coupling_matrix();

    //! Solve the linear system of equations
    void linear_solve();

    //! Create the linear solver
    virtual void create_linear_solver(const Teuchos::ParameterList& solverparams,
        const Core::LinearSolver::SolverType solvertype);

    //! Setup Newton-Raphson
    void setup_newton();

    //! Print Header to screen
    virtual void print_header();

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables)
    void update_fields_after_convergence();

    //! build norms for convergence check
    virtual void build_convergence_norms();

    void poro_fd_check();

    // check for convergence
    bool converged();

    /// Print user output that structure field is disabled
    void print_structure_disabled_info();

    //! convergence tolerance for increments
    double ittolinc_;
    //! convergence tolerance for residuals
    double ittolres_;
    //! maximally permitted iterations
    int itmax_;
    //! minimally necessary iterations
    int itmin_;
    //! current iteration step
    int itnum_;
    //! @name Global vectors
    Teuchos::RCP<Core::LinAlg::Vector> zeros_;  //!< a zero vector of full length

    Teuchos::RCP<Core::LinAlg::Vector> iterinc_;  //!< increment between Newton steps k and k+1
    //!< \f$\Delta{x}^{<k>}_{n+1}\f$

    Teuchos::RCP<Core::LinAlg::Vector> rhs_;  //!< rhs of Poroelasticity system

    Teuchos::RCP<Core::LinAlg::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ?
    bool solveradapttol_;                        //!< adapt solver tolerance


    //@}

    //! @name Global matrixes

    //! block systemmatrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    //! structure-fluid coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_sf_;
    //! fluid-structure coupling matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs_;

    //@}

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<Core::LinAlg::Equilibration> equilibration_;

    //! equilibration method applied to system matrix
    Core::LinAlg::EquilibrationMethod equilibration_method_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Epetra_Map> combinedDBCMap_;

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_fluid_;   //!< tolerance residual increment for fluid
    double tolfres_fluid_;  //!< tolerance force residual for fluid

    double normrhs_;  //!< norm of residual forces

    double normrhsfluid_;  //!< norm of residual forces (fluid )
    double normincfluid_;  //!< norm of residual unknowns (fluid )

    double normrhsstruct_;  //!< norm of residual forces (structure)
    double normincstruct_;  //!< norm of residual unknowns (structure)

    double normrhsart_;       //!< norm of residual (artery)
    double normincart_;       //!< norm of residual unknowns (artery)
    double arterypressnorm_;  //!< norm of artery pressure

    double maxinc_;  //!< maximum increment
    double maxres_;  //!< maximum residual

    enum Inpar::POROMULTIPHASE::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum Inpar::POROMULTIPHASE::VectorNorm vectornorminc_;   //!< type of norm for increments

    Teuchos::Time timernewton_;  //!< timer for measurement of solution time of newton iterations
    double dtsolve_;             //!< linear solver time
    double dtele_;               //!< time for element evaluation + build-up of system matrix

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<Core::Conditions::LocsysManager> locsysman_;

    //! flag for finite difference check
    Inpar::POROMULTIPHASE::FdCheck fdcheck_;

  };  // PoroMultiPhasePartitioned

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithicTwoWayArteryCoupling : public PoroMultiPhaseMonolithicTwoWay
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithicTwoWayArteryCoupling(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! extract the field vectors from a given composed vector.
    /*!
     x is the sum of all increments up to this point.
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations,
     and 1D artery pressure)
     */
    void extract_field_vectors(Teuchos::RCP<const Core::LinAlg::Vector> x,
        Teuchos::RCP<const Core::LinAlg::Vector>& sx,
        Teuchos::RCP<const Core::LinAlg::Vector>& fx) override;

    //! build norms for convergence check
    void build_convergence_norms() override;

   protected:
    /// setup composed system matrix from field solvers
    void setup_maps() override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    /// setup global rhs
    void setup_rhs() override;

    //! build the combined dirichletbcmap
    void build_combined_dbc_map() override;

    //! Create the linear solver
    void create_linear_solver(const Teuchos::ParameterList& solverparams,
        const Core::LinearSolver::SolverType solvertype) override;

    //! build the block null spaces
    void build_artery_block_null_space(
        Teuchos::RCP<Core::LinAlg::Solver>& solver, const int& arteryblocknum) override;

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_artporo_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockrowdofmap_artporo_;
  };


}  // namespace POROMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif
