// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_method_direct.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method_projector.hpp"

#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <magic_enum/magic_enum.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! persistent SuperLU_DIST process grid, stored behind the opaque handle in DirectSolver
  struct SuperluDistGrid
  {
    SLUD::gridinfo_t grid;
  };
}  // namespace

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::DirectSolver::DirectSolver(Core::LinearSolver::SolverType solvertype)
    : solvertype_(solvertype), factored_(false), solver_(nullptr), projector_(nullptr)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::DirectSolver::~DirectSolver()
{
  if (superlu_grid_ != nullptr)
  {
    auto* grid_holder = static_cast<SuperluDistGrid*>(superlu_grid_);
    SLUD::superlu_gridexit(&grid_holder->grid);
    delete grid_holder;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::DirectSolver::setup(std::shared_ptr<Core::LinAlg::SparseOperator> matrix,
    std::shared_ptr<Core::LinAlg::MultiVector<double>> b, const bool refactor, const bool reset,
    std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector)
{
  // 1. project the linear system if close to being singular and set the final matrix and vectors
  projector_ = projector;
  std::shared_ptr<Core::LinAlg::MultiVector<double>> rhs = b;
  if (projector_ != nullptr)
  {
    FOUR_C_ASSERT_ALWAYS(b->num_vectors() == 1,
        "Expecting only one solution vector during projector call! Got {} vectors.",
        b->num_vectors());
    LinAlg::Vector<double> projected_b = projector_->to_reduced(b->get_vector(0));
    rhs = std::make_shared<Core::LinAlg::MultiVector<double>>(projected_b.as_multi_vector());

    matrix = projector_->to_reduced(*matrix);
  }


  // 2. merge the block system matrix into a standard sparse matrix if necessary
  auto crsA = std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(matrix);
  if (!crsA)
  {
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> Ablock =
        std::dynamic_pointer_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

    int matrixDim = Ablock->full_range_map().num_global_elements();
    if (matrixDim > 50000 and Communication::my_mpi_rank(matrix->domain_map().get_comm()) == 0)
      std::cout << "\n WARNING: Direct linear solver is merging matrix, this is very expensive! \n";

    crsA = Ablock->merge();
  }

  b_ = rhs;
  a_ = crsA;

  // SuperLU_DIST is driven directly through its native interface, see solve_superlu_dist()
  if (solvertype_ == SolverType::Superlu)
  {
    factored_ = false;
    return;
  }

  // 3. create linear solver
  if (reset or refactor or not is_factored())
  {
    std::string solver_type;
    Teuchos::ParameterList params("Amesos2");

    switch (solvertype_)
    {
      case SolverType::KLU2:
      {
        solver_type = "KLU2";
        auto& klu_params = params.sublist(solver_type);
        klu_params.set("IsContiguous", false, "Are GIDs Contiguous");
        break;
      }
      case SolverType::MUMPS:
      {
        solver_type = "MUMPS";
        auto& mumps_params = params.sublist(solver_type);
        mumps_params.set("IsContiguous", false, "Are GIDs Contiguous");
        break;
      }
      case SolverType::UMFPACK:
      {
        solver_type = "Umfpack";
        auto& umfpack_params = params.sublist(solver_type);
        umfpack_params.set("IsContiguous", false, "Are GIDs Contiguous");
        break;
      }
      default:
        FOUR_C_THROW("Unsupported solver type {}!", magic_enum::enum_name(solvertype_));
    }

    FOUR_C_ASSERT_ALWAYS(Amesos2::query(solver_type),
        "Requested direct solver {} is not available in Amesos2! Check your Amesos2 installation "
        "and choose an available solver!",
        solver_type);

    solver_ = Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>(
        solver_type, Teuchos::rcpFromRef(a_->epetra_matrix()));
    solver_->setB(Teuchos::rcpFromRef(b_->get_epetra_multi_vector()));

    solver_->setParameters(Teuchos::make_rcp<Teuchos::ParameterList>(std::move(params)));

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int Core::LinearSolver::DirectSolver::solve_superlu_dist(Core::LinAlg::MultiVector<double>& x)
{
  FOUR_C_ASSERT_ALWAYS(x.num_vectors() == 1,
      "The SuperLU_DIST interface supports a single right-hand side, got {} vectors.",
      x.num_vectors());

  const Epetra_CrsMatrix& matrix = a_->epetra_matrix();
  const auto* mpi_comm = dynamic_cast<const Epetra_MpiComm*>(&matrix.Comm());
  FOUR_C_ASSERT_ALWAYS(mpi_comm != nullptr, "SuperLU_DIST requires an MPI communicator.");

  const int num_procs = mpi_comm->NumProc();
  const int num_global_rows = matrix.NumGlobalRows();

  // 1. process grid spanning all ranks (created once per solver object)
  if (superlu_grid_ == nullptr)
  {
    int nprow = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(num_procs))));
    while (num_procs % nprow != 0) --nprow;
    const int npcol = num_procs / nprow;

    auto* grid_holder = new SuperluDistGrid;
    SLUD::superlu_gridinit(mpi_comm->Comm(), nprow, npcol, &grid_holder->grid);
    superlu_grid_ = grid_holder;
  }
  SLUD::gridinfo_t& grid = static_cast<SuperluDistGrid*>(superlu_grid_)->grid;

  // 2. redistribute the system to the contiguous block-row layout SuperLU_DIST requires
  const Epetra_Map linear_map(num_global_rows, 0, *mpi_comm);
  const Epetra_Import to_linear(linear_map, matrix.RowMap());

  Epetra_CrsMatrix linear_matrix(::Copy, linear_map, 0);
  linear_matrix.Import(matrix, to_linear, ::Insert);
  linear_matrix.FillComplete(linear_map, linear_map);

  Epetra_MultiVector linear_b(linear_map, 1);
  linear_b.Import(b_->get_epetra_multi_vector(), to_linear, ::Insert);

  const int m_loc = linear_map.NumMyElements();
  const int fst_row = m_loc > 0 ? linear_map.MinMyGID() : 0;
  const int nnz_loc = linear_matrix.NumMyNonzeros();

  // SuperLU_DIST takes ownership of these arrays (freed by Destroy_CompRowLoc_Matrix_dist)
  SLUD::int_t* rowptr = SLUD::intMalloc_dist(m_loc + 1);
  SLUD::int_t* colind = SLUD::intMalloc_dist(std::max(nnz_loc, 1));
  double* nzval = SLUD::D::doubleMalloc_dist(std::max(nnz_loc, 1));

  {
    std::vector<std::pair<int, double>> row_entries;
    int pos = 0;
    rowptr[0] = 0;
    for (int lrow = 0; lrow < m_loc; ++lrow)
    {
      int len = 0;
      double* values = nullptr;
      int* local_cols = nullptr;
      linear_matrix.ExtractMyRowView(lrow, len, values, local_cols);

      row_entries.clear();
      for (int k = 0; k < len; ++k)
        row_entries.emplace_back(linear_matrix.ColMap().GID(local_cols[k]), values[k]);
      std::sort(row_entries.begin(), row_entries.end());

      for (const auto& [gcol, value] : row_entries)
      {
        colind[pos] = gcol;
        nzval[pos] = value;
        ++pos;
      }
      rowptr[lrow + 1] = pos;
    }
  }

  SLUD::SuperMatrix slu_matrix;
  SLUD::D::dCreate_CompRowLoc_Matrix_dist(&slu_matrix, num_global_rows, num_global_rows, nnz_loc,
      m_loc, fst_row, nzval, colind, rowptr, SLUD::SLU_NR_loc, SLUD::SLU_D, SLUD::SLU_GE);

  // 3. solve with equilibration, MC64 row permutation and double iterative refinement --
  //    the combination needed to get accurate solutions for ill-conditioned systems
  SLUD::amesos2_superlu_dist_options_t options;
  SLUD::set_default_options_dist(&options);
  options.IterRefine = SLUD::SLU_DOUBLE;
  options.PrintStat = SLUD::NO;

  SLUD::D::dScalePermstruct_t scale_perm;
  SLUD::D::dScalePermstructInit(num_global_rows, num_global_rows, &scale_perm);
  SLUD::D::dLUstruct_t lu;
  SLUD::D::dLUstructInit(num_global_rows, &lu);
  SLUD::D::dSOLVEstruct_t solve_struct;
  SLUD::SuperLUStat_t stat;
  SLUD::PStatInit(&stat);

  std::vector<double> rhs_and_solution(std::max(m_loc, 1));
  for (int i = 0; i < m_loc; ++i) rhs_and_solution[i] = linear_b[0][i];

  double berr = 0.0;
  int info = 0;
  SLUD::D::pdgssvx(&options, &slu_matrix, &scale_perm, rhs_and_solution.data(), std::max(m_loc, 1),
      1, &grid, &lu, &solve_struct, &berr, &stat, &info);

  SLUD::PStatFree(&stat);
  SLUD::D::dScalePermstructFree(&scale_perm);
  SLUD::D::dDestroy_LU(num_global_rows, &grid, &lu);
  SLUD::D::dLUstructFree(&lu);
  if (options.SolveInitialized) SLUD::D::dSolveFinalize(&options, &solve_struct);
  SLUD::Destroy_CompRowLoc_Matrix_dist(&slu_matrix);

  FOUR_C_ASSERT_ALWAYS(info == 0, "SuperLU_DIST pdgssvx failed with info = {}.", info);

  // 4. map the solution back to the original layout
  for (int i = 0; i < m_loc; ++i) linear_b[0][i] = rhs_and_solution[i];
  x.get_epetra_multi_vector().Export(linear_b, to_linear, ::Insert);

  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int Core::LinearSolver::DirectSolver::solve(Core::LinAlg::MultiVector<double>& x)
{
  std::unique_ptr<LinAlg::Vector<double>> projected_x = nullptr;
  if (projector_ != nullptr)
  {
    FOUR_C_ASSERT_ALWAYS(x.num_vectors() == 1,
        "Expecting only one solution vector during projector call! Got {} vectors.",
        x.num_vectors());

    // Create empty x in reduced space
    projected_x = std::make_unique<LinAlg::Vector<double>>(a_->domain_map(), true);
  }
  Core::LinAlg::MultiVector<double>& solution =
      projector_ != nullptr ? projected_x->as_multi_vector() : x;

  if (solvertype_ == SolverType::Superlu)
  {
    solve_superlu_dist(solution);
    factored_ = true;
  }
  else
  {
    solver_->setX(Teuchos::rcpFromRef(solution.get_epetra_multi_vector()));

    if (not is_factored())
    {
      solver_->symbolicFactorization();
      solver_->numericFactorization();

      factored_ = true;
    }

    solver_->solve();
  }

  if (projector_ != nullptr)
  {
    x.get_vector(0) = projector_->to_full(*projected_x);
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE