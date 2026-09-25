# Linear buckling analysis for 4C — implementation plan

**Target:** an Abaqus `*BUCKLE`-equivalent eigenvalue buckling procedure for 3D solids (`SOLID HEX8/20/27, TET4/10, WEDGE6, PYRAMID5`), including RVEs with periodic BCs.

**Formulation:** the classical one used by Abaqus and by the FEniCSx benchmark (verification case C). `K_Δ` is the **pure initial-stress (geometric) stiffness** of the perturbation stress. There is no finite-difference tangent.

**Verification:** cases A–D in [`verification/buckling_verification_cases.pdf`](verification/buckling_verification_cases.pdf). Reference values and how they were produced are in [`verification/reference_values.md`](verification/reference_values.md).

Paths are relative to the repo root; `SN` = `src/structure_new/src`.

### Design choices that must match Abaqus / FEniCSx (they determine the results)

| # | Design choice | Consequence |
|---|---|---|
| D1 | `K_Δ` = initial-stress stiffness of the **perturbation stress only**. The initial-displacement matrix and the change of the material tangent are not included. | Same λ as Abaqus/FEniCSx. Case C must match to 1e-4. |
| D2 | Perturbation stress from a **linear perturbation solve at the base state** (`K₀ Δu = Q`) | λ independent of the magnitude of Q; λ(−Q) = −λ(Q) |
| D3 | **Preload = separate base state:** the full nonlinear tangent `K₀ = K_T(u₀)` includes the preload stress; the preload is *not* scaled by λ | critical load = preload + λ·Q (Abaqus general step + `*BUCKLE`) |
| D4 | Dead loads only in `K_Δ` (no follower-load stiffness in the first version) | Same as FEniCSx. Abaqus adds load stiffness for pressure loads, so this is a documented difference of the first version |
| D5 | Modes normalized to max. displacement component = 1; negative λ reported | Same output convention as Abaqus |

**Free implementation choices (4C-native, not copied from anywhere):** how `K_Δ` is computed in 4C, the eigensolver, Dirichlet handling, input format, output. Section 2 describes the current proposal. It can be changed as long as D1–D5 hold and verification cases A–D pass.

### Rejected alternative: finite-difference `K_Δ`

An earlier draft computed `K_Δ ≈ [K_T(u₀ + εΔu) − K_T(u₀)]/ε`. It works with every element technology without new element code, but it also contains the initial-displacement matrix. It therefore violates D1: `λ_i` shifts by roughly the pre-buckling strain `λ_i/E`, i.e. up to 9e-3 for mode 8 of case C, and cannot match Abaqus/FEniCSx. It may still serve as a debugging aid in unit tests, but it is not part of the feature.

---

## 1. Formulation (identical to Abaqus `*BUCKLE` and FEniCSx)

**Base state u₀.** Either `u₀ = 0` (`NUMSTEP: 0`) or the converged result of the nonlinear static preload steps, like the Abaqus general step before `*BUCKLE`.

**Base stiffness.** `K₀ = K_T(u₀)` is the full consistent tangent at the base state: material part, geometric part of the *preload* stress, and all model evaluators (periodic-BC penalty, springs). At `u₀ = 0` it equals the linear stiffness `K` of FEniCSx.

**Perturbation load pattern Q** (Abaqus: the loads defined inside the buckle step):
- Neumann part: `Q_N = f_ext(t_pert) − f_ext(t_base)`.
- Dirichlet part: `Δu_D = u_D(t_pert) − u_D(t_base)`. This is needed for RVE macro-strain loading.

**Linear perturbation solve at the base state** (Abaqus: "the perturbation stresses are computed from a linear perturbation analysis"; FEniCSx: `u = problem.solve()`):

    K₀ Δu = Q_N − K₀ Δu_D   on free DOFs,   Δu = Δu_D on Dirichlet DOFs

**Perturbation stress** at every Gauss point (total Lagrangian):

    ΔE = B(u₀) Δd_e              (linearized Green–Lagrange strain; B = strain gradient at F₀)
    ΔS = ℂ(u₀) : ΔE              (ℂ = material tangent dS/dE at the base state)

At `u₀ = 0` this is exactly `σ₀ = ℂ : ε(Δu)` of the FEniCSx tour.

**Differential stiffness.** This is the Abaqus "differential initial stress matrix" and the FEniCSx `K_G`:

    K_Δ = Σ_e Σ_gp  w · N_X · ΔS · N_Xᵀ        (only the stress term — no ΔF·ℂ·F terms, no Δℂ)

**Eigenproblem.**

    (K₀ + λ K_Δ) φ = 0 ,    critical load = base load + λ_i · Q

Negative `λ` means the load pattern is reversed, as in Abaqus.

**Dirichlet DOFs** (FEniCSx convention): 1 on the diagonal of `K₀`, 0 on rows and columns of `K_Δ`. These DOFs then give `λ = ∞` and never pollute the spectrum. With this, `K₀ = K_T(0)` and the same K_Δ reproduce FEniCSx case C up to the quadrature rule.

**Not included** (the last two are omitted by Abaqus and FEniCSx as well; the first is the documented limitation D4):
- follower-load stiffness (pressure load stiffness);
- the initial-displacement matrix;
- the change of `ℂ` with the perturbation.

---

## 2. Architecture

Analysis type `DYNAMICTYPE: "LinearBuckling"`. The static steps before the buckling step form the preload. The buckling itself is a procedure, not a model evaluator, but it **uses** all model evaluators to assemble `K₀`. Only the structure evaluator contributes `K_Δ`: periodic constraints and springs are linear and have no initial-stress term.

```
caldyn_drt()                              src/structure/4C_structure_dyn_nln_drt.cpp:40
  └─ case LinearBuckling → Solid::linear_buckling_drt()                 (new)
       ├─ Adapter::build_structure_algorithm (as dyn_nlnstructural_drt), Statics integrator
       ├─ integrate()                      → preload (0 steps allowed) → u₀
       └─ Solid::Buckling::LinearBucklingAnalysis::run()                (new)
            1  K₀ = K_T(u₀)                  model evaluator manager (all models)
            2  Q_N, Δu_D                     f_ext / Dbc at t_pert vs t_base
            3  Δu = K₀⁻¹(...)                Core::LinAlg::Solver, factorization kept
            4  K_Δ                           ModelEvaluator::Structure::evaluate_differential_stiffness(Δu)
                                               → element action struct_calc_differential_stiffness
            5  eigenpairs                    Core::LinAlg::SubspaceEigenSolver (reuses factorization of K₀)
            6  post                          sort, normalize, table, CSV, VTU, result test
```

### 2.1 Element level (new, core of the feature)

| File | Change |
|---|---|
| `src/core/legacy_enum_definitions/4C_legacy_enum_definitions_element_actions.hpp` | new `ActionType::struct_calc_differential_stiffness` + string conversion |
| `src/solid_ele/4C_solid_ele_evaluate.cpp` (switch at ~93) | new case → `std::visit(... interface->evaluate_differential_stiffness(...))`, result in `elemat1` |
| `src/solid_ele/4C_solid_ele_calc.hpp/.cpp` | new `SolidEleCalc::evaluate_differential_stiffness(ele, material, discretization, lm, params, K_delta)` |
| `src/solid_ele/4C_solid_ele_calc_lib.hpp` | small helper `evaluate_perturbation_gl_strain(Bop, Δd_e)`; reuse `evaluate_strain_gradient` (l. 777) and `add_geometric_stiffness_matrix` (l. 1189) |

Gauss-point loop, mirroring `evaluate_nonlinear_force_stiffness_mass` (`4C_solid_ele_calc.cpp:139`):

```cpp
// states: "displacement" = u₀, "buckling perturbation displacement" = Δu
for_each_gauss_point(nodal_coordinates, element_properties_, stiffness_matrix_integration_,
  [&](xi, shape_functions, jacobian_mapping, integration_factor, gp) {
    evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
             history_data_, gp, [&](const auto& F0, const auto& E0, const auto& linearization) {
      const Stress<celltype> stress0 = evaluate_material_stress<celltype>(solid_material, ..., F0, E0, ...);
      const auto Bop = evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping(F0));
      const auto dE  = Bop * delta_d_e;                         // Voigt, strain-like
      const auto dS  = stress0.cmat_ : dE;                      // perturbation PK2
      add_geometric_stiffness_matrix(jacobian_mapping, dS, integration_factor, K_delta);
    });
  });
```

**Element technology coverage:**
- **First version (M1):** `ElementTechnology::none` with `KINEM nonlinear`, all cell types. `if constexpr` dispatch on the formulation type; other formulations throw
  `"Linear buckling: differential stiffness not implemented for element technology <x>"`.
- **M4:** EAS (`eas_mild/eas_full`). `ΔE` also contains the enhanced part `Δα = −K_αα⁻¹ K_αd Δd_e`, using the condensation matrices already built in `4C_solid_ele_calc_eas_helpers.hpp`. This matters because `HEX8 + eas_full` is the cheap element for thin plates.
- **Later:** F-bar (`Bop`, `Hop`, `fbar_factor` exist in `4C_solid_ele_calc_fbar.hpp`), shell-ANS, MULF.
- `KINEM linear` is rejected (no geometric stiffness exists in that formulation).

The material must provide `cmat_` at `u₀`. That is standard for `So3Material::evaluate`. `evaluate_material_stress` must not update material history: it is called with the converged `u₀` and without an update step.

### 2.2 Structure / driver level

| File | Change |
|---|---|
| `SN/input/4C_structure_new_input.hpp/.cpp` | `DynamicType::LinearBuckling` (+ selection string at ~86); group `STRUCTURAL DYNAMIC/LINEAR BUCKLING` (next to `/GENALPHA` ~600) |
| `SN/4C_structure_new_timint_factory.cpp` (69, 113), `SN/4C_structure_new_factory.cpp` (67) | map `LinearBuckling` → `TimeInt::Implicit` + `IMPLICIT::Statics` |
| `src/structure/4C_structure_dyn_nln_drt.cpp` | case → `linear_buckling_drt()` |
| `SN/model_evaluator/4C_structure_new_model_evaluator_structure.hpp/.cpp` | `evaluate_differential_stiffness(const Vector& delta_u, SparseMatrix& K_delta)`: action param, `set_state` of u₀ and Δu, `discret().evaluate(...)`, `complete()`; `evaluate_neumann` reused for `Q_N` |
| **new** `SN/buckling/4C_structure_new_buckling_analysis.hpp/.cpp` | `Solid::Buckling::LinearBucklingAnalysis` (steps 1–6) |
| **new** `SN/buckling/4C_structure_new_buckling_input.hpp/.cpp` | parameter struct + reader |
| **new** `src/core/linalg/src/sparse/4C_linalg_subspace_eigen_solver.hpp/.cpp` | generic eigensolver (§3) |
| `src/core/linalg/src/dense/4C_linalg_utils_densematrix_eigen.hpp/.cpp` | `symmetric_generalized_eigen(A, B, λ, V)` = LAPACK `SYGV` wrapper |
| `SN/utils/4C_structure_new_resulttest.cpp` | `SPECIAL` quantities `buckling_eigenvalue_<i>` |
| **new** `SN/buckling/CMakeLists.txt` | `four_c_auto_define_module()` (sources are picked up automatically, as in `SN/implicit/`) |

**Model-evaluator guard.** Allowed:
- `model_structure`;
- `model_constraints` (periodic RVE, penalty);
- `model_springdashpot`.

Anything else (contact, meshtying, beams, 0D, …) throws a clear error.

### 2.3 Input

```yaml
STRUCTURAL DYNAMIC:
  DYNAMICTYPE: "LinearBuckling"
  NUMSTEP: 0              # >0: nonlinear static preload (Abaqus: general step before *BUCKLE)
  MAXTIME: 0.0            # t_base
  TIMESTEP: 1.0
  LINEAR_SOLVER: 1
STRUCTURAL DYNAMIC/LINEAR BUCKLING:
  NUM_MODES: 10           # Abaqus: number of eigenvalues requested
  PERTURBATION_TIME: 1.0  # Q = loads(t_pert) - loads(t_base), Neumann + Dirichlet
  SHIFT: 0.0              # Abaqus: minimum eigenvalue of interest / FEniCSx shift
  TOLERANCE: 1.0e-10
  MAX_ITER: 300
  SUBSPACE_SIZE: 0        # 0 = automatic min(2n, n+8)
  MODE_NORMALIZATION: max_displacement   # Abaqus default | stiffness (phi^T K0 phi = 1)
  WRITE_CSV: true
RESULT DESCRIPTION:
  - STRUCTURE: {DIS: "structure", SPECIAL: true, QUANTITY: "buckling_eigenvalue_1",
                VALUE: 0.167963169, TOLERANCE: 1.0e-5}
```

**Output:**
- A screen table, like the Abaqus `.dat` eigenvalue output: mode, `λ_i`, critical-load factor, `negative` flag, dominant component.
- `<prefix>_buckling.csv`.
- VTU with one step per mode (`time = step = i`, field `displacement` = `φ_i`, normalized so that max |component| = 1). These appear as frames in ParaView, like the Abaqus buckle frames.

---

## 3. Eigensolver `Core::LinAlg::SubspaceEigenSolver` (no new dependency)

It solves `K_Δ φ = μ (−K₀) φ` with `μ = 1/(λ − σ)`, wanting the largest `|μ|`. With shift σ, `K₀` is replaced by `K₀ + σ K_Δ`, which is the same spectral transformation as the FEniCSx tour.

Bathe subspace iteration, about 300 lines:
- **Vectors:** `Epetra`-backed `MultiVector`, so it runs MPI-parallel.
- **Linear solve:** one factorization via `Core::LinAlg::Solver` with `refactor=false` after the first solve, then `solve_with_multi_vector`.
- **Projected problem:** dense `SYGV`, B-orthonormalization, convergence on the relative change of `μ_i`.
- **Start vectors:** random, with Dirichlet rows zeroed.

The same class can later serve natural-frequency analysis (`K φ = ω² M φ`).

---

## 4. Implementation steps

### M1 — core formulation + FEniCSx code-to-code (the strictest test first)
1. Element action + `SolidEleCalc::evaluate_differential_stiffness` (displacement-based).
2. **Unit test** (`unittests/solid_ele/`, the existing location of the `solid_ele` unit tests): one `HEX8`/`HEX27` element at `u₀ = 0` under a homogeneous Δu (uniaxial strain). Check that `ΔS = ℂ:ΔE` exactly and that `K_Δ` equals `add_geometric_stiffness_matrix` with that stress. Also check symmetry and zero rigid-translation energy.
3. `ModelEvaluator::Structure::evaluate_differential_stiffness`.
4. `SubspaceEigenSolver` + `SYGV` wrapper, with unit tests in `src/core/linalg/tests/` (extend `4C_linalg_utils_densematrix_eigen_test.cpp`; new `4C_linalg_subspace_eigen_solver_test.cpp` + `.np2` variant) against dense LAPACK on a random SPD/indefinite pair.
5. Input section, `DynamicType`, factories, driver (Neumann pattern only, `u₀ = 0`), screen table, `buckling_eigenvalue_i` result test.
6. **Test C** `buckling_fenics_benchmark_hex27.4C.yaml`:
   - `50×5×5 HEX27`, `E = 1000`, `ν = 0.3`;
   - clamped at x=0, `u_y = u_z = 0` at x=L, traction `(−1,0,0)`;
   - **modes 1–8 within 1e-4 of the FEniCSx values** (0.167963169, 0.496955439, 0.987888956, 1.500093336, …).
7. **Test B** cantilever `HEX20`: λ₁ = 103.605 N etc. within 1 % and 2 %, mode order weak/strong/weak/strong/weak/weak.

### M2 — output, Dirichlet pattern, plate and periodic tests
1. VTU mode output, CSV, normalization, negative-eigenvalue flag.
2. Dirichlet perturbation pattern `Δu_D`.
3. **Test A1/A2** SS plate `HEX20`, 2 layers, mid-plane line supports:
   - A1: λ₁ = 75.877 N/mm (Mindlin), 6 modes;
   - A2: 8 closely spaced modes in exact `(m,n)` order;
   - load-reversal check.
4. **Test D** periodic lattice D1 (1×1) / D2 (2×2):
   - identities I1–I6;
   - plausibility λ ≈ 10.0 / 2.5.

### M3 — preload, shift, parallel
1. `NUMSTEP > 0` base state; history-safe material evaluation at `u₀`.
2. **Test B-preload:** `P_pre = 51.80 N` ⇒ `P_pre + λ₁ = 103.605` within 1e-3.
3. `SHIFT`.
4. `NP 2` variants of A/B/C (D stays serial, an existing RVE limitation).
5. User documentation page (`doc/`).

### M4 — element technologies
EAS differential stiffness, plus test A1 with `HEX8 eas_full` (4 layers). Then F-bar.

### M5 — imperfections
Write modes to a file, plus an input option to add `Σ aᵢ φᵢ` to the reference geometry for a follow-up nonlinear run (Abaqus `*IMPERFECTION`).

**Work cycle per step:**
1. code;
2. incremental build;
3. new tests + neighbouring tests (`rve3d_periodic_bcs*`, `solid_ele_*`);
4. `pre-commit`/clang-format;
5. commit and push to `dev-linear-buckling` (see [`development_workflow.md`](development_workflow.md)).

---

## 5. Test ↔ verification-document mapping

| 4C test file (planned) | Doc case | Hard check |
|---|---|---|
| `buckling_fenics_benchmark_hex27.4C.yaml` | C | modes 1–8 vs FEniCSx 50×5×5 within 1e-4 |
| `buckling_fenics_benchmark_hex27_fine.4C.yaml` (nightly) | C | vs FEniCSx 80×8×8 within 1e-4 |
| `buckling_cantilever_hex20.4C.yaml` | B | λ₁…λ₆ within 1 %/2 %, mode order |
| `buckling_cantilever_hex20_preload.4C.yaml` | B | preload consistency 1e-3 |
| `buckling_plate_ss_square_hex20.4C.yaml` | A1 | 6 modes within 1 %/2 % of Mindlin |
| `buckling_plate_ss_long_hex20.4C.yaml` | A2 | 8 modes, `(m,n)` order |
| `buckling_rve_lattice_1x1.4C.yaml`, `..._2x2.4C.yaml` | D | I1–I6 |

Regression tolerances are tightened to about 1e-8 on the 4C values once the first verified results exist.

---

## 6. Build / test environment

See [`development_workflow.md`](development_workflow.md): build in the CI dependency image, run the targeted tests locally, full regression via a PR inside the fork `m-frey/4C` before any upstream PR.

---

## 7. Code style and mergeability (upstream 4C `main`)

**Goal:** the feature reads like existing 4C code and can be merged upstream without restructuring. The rules below come from the repository itself:
- `CONTRIBUTING.md`;
- `doc/documentation/src/developer_guide/codingguidelines.rst` and `testing.rst`;
- `.clang-format`, `.clang-tidy`, `.pre-commit-config.yaml`;
- `utilities/code_checks/commit-msg`.

When in doubt, the neighbouring code in `src/structure_new` and `src/solid_ele` is the model.

### 7.1 Enforced automatically (pre-commit + CI) — must pass before every push
- **Setup:** `./utilities/set_up_dev_env.sh` installs the hooks. Every commit then runs through `pre-commit`.
- **Formatting:**
  - C++: `clang-format` (Google base style, Allman braces, 100 columns, `&`/`*` bound to the type, regrouped includes).
  - Also: `cmake-format`, `yamlfmt` for the test input files, `black` for Python, trailing whitespace, no non-ASCII characters, `typos`.
- **Headers:**
  - LGPL license header (inserted by the `insert-license` hook);
  - header guard `FOUR_C_<PATH>_HPP` (checked);
  - `#include "4C_config.hpp"` first;
  - code between `FOUR_C_NAMESPACE_OPEN` / `FOUR_C_NAMESPACE_CLOSE`.
- **File names:** `4C_<module>_<name>.{hpp,cpp}` (`check-filenames`). For example:
  - `4C_structure_new_buckling_analysis.cpp`;
  - `4C_linalg_subspace_eigen_solver.hpp`;
  - unit tests `<file>_test.cpp` in the module's `tests/` directory (`<file>_test.np2.cpp` for MPI).
- **Input files:** every new `tests/input_files/*.4C.yaml` must be registered in `tests/list_of_tests.cmake` (`check-test-files`). No orphaned companion files.
- **Includes:** no include cycles (`check-includes-for-cycles`). Respect the declared module dependencies (`four_c_add_internal_dependency`): `solid_ele` must not include `structure_new`, and the new `buckling/` directory gets its own `CMakeLists.txt` with `four_c_auto_define_module()`.
- **Commit messages** (`commit-msg` hook):
  - subject ≤ 50 characters, capitalized, imperative mood, no trailing period;
  - blank line before the body;
  - no single-word subjects;
  - style like the recent history, e.g. `Add differential stiffness to solid elements`.

### 7.2 Coding guidelines (reviewed by maintainers)
- **Naming** (`.clang-tidy`):
  - `CamelCase` for namespaces, classes, structs, enums, concepts;
  - `snake_case` functions;
  - private members with trailing `_`;
  - `lower_case` value template parameters, `CamelCase` type template parameters;
  - no single-letter variable names except loop indices.
- **Modern C++20:** `std::shared_ptr`/`std::unique_ptr` only where ownership is involved. Otherwise pass by (const) reference. No new `Teuchos::RCP`. Const-correct interfaces.
- **No new `#define` flags.** Debug checks go behind `FOUR_C_ENABLE_ASSERTIONS` (`FOUR_C_ASSERT`); user errors use `FOUR_C_THROW` with an actionable message.
- **Parameters:**
  - The buckling settings are read once from the input into a plain `struct` (`Solid::Buckling::Parameters`), following the guideline "prefer parameter container classes over `Teuchos::ParameterList`".
  - The `ParameterList` is only used where the existing element-evaluate interface requires it (action + states).
- **Enums:**
  - `enum class` everywhere (`ModeNormalization`, …).
  - New input selections use `parameter<EnumType>(...)` with `EnumTools`. **Not** `deprecated_selection` (marked deprecated in `4C_io_input_spec_builders.hpp`).
  - `"LinearBuckling"` is added to the existing `DYNAMICTYPE` selection without restructuring it.
- **Header hygiene:** forward declarations instead of header-in-header includes. Expensive templates go into `.templates.hpp`/`.inst.hpp` only if the existing element code does the same.
- **Doxygen** on every new public class and function (`/*! \brief ... */` or `///`, as in the surrounding files). Short comments only where the reason is not obvious.

### 7.3 Fitting into existing architecture (no parallel infrastructure)
- **Element:**
  - The new action is added to the existing `Core::Elements::ActionType` list and dispatched in `4C_solid_ele_evaluate.cpp` exactly like `struct_calc_nlnstiff` (`std::visit` on `solid_calc_variant_`).
  - The Gauss loop reuses `for_each_gauss_point`, `evaluate`, `evaluate_material_stress`, `evaluate_strain_gradient` and `add_geometric_stiffness_matrix`.
  - No copies of existing kinematics code.
- **Unsupported technologies:** they are rejected with `if constexpr` + `FOUR_C_THROW`, like other formulation-specific features in `solid_ele`. They are not silently ignored.
- **Driver:**
  - reuses `Adapter::build_structure_algorithm`, the `Statics` integrator, `Solid::Dbc`, the model evaluator manager, `Core::LinAlg::Solver` and `DiscretizationVisualizationWriterMesh`;
  - `caldyn_drt` gets one additional `case`, no refactoring of unrelated code.
- **Eigensolver:** it lives in `core/linalg` as a generic, documented utility with its own unit test, so it can be reused (e.g. modal analysis). No buckling specifics inside it.
- **Scope:** no changes to unrelated modules. Existing behaviour and all existing tests stay unchanged. The only visible changes are the new `DYNAMICTYPE` value and the new input section.

### 7.4 Tests and documentation (required for review)
- **Unit tests** next to the tested code, following the neighbours: `unittests/solid_ele/` (legacy location still used by `solid_ele`) and `src/core/linalg/tests/` (auto-registered by `four_c_auto_define_tests()`):
  - the element differential stiffness;
  - the eigensolver.
- **Input-file tests A–D** following `testing.rst`: a meaningful `TITLE`, small meshes, a `RESULT DESCRIPTION` with explicit tolerances, `NP` set in `list_of_tests.cmake`. Large meshes (e.g. the FEniCSx 80×8×8 case) go to the nightly tests or are omitted.
- **Documentation:**
  - A short page in the Sphinx analysis guide (`doc/documentation/src/analysis_guide_templates`) covering theory (D1–D5), input and output.
  - The input-parameter reference is generated automatically from the `InputSpec` descriptions, so every new parameter needs a precise `.description`.

### 7.5 Pull-request strategy
- **Commits:**
  - small, self-contained and reviewable in isolation;
  - each builds and passes tests;
  - the history is kept clean (upstream merges with merge commits, no squash).
  - Suggested order:
    1. dense `SYGV` wrapper;
    2. subspace eigensolver + tests;
    3. element differential stiffness + unit test;
    4. structure model evaluator hook;
    5. input + driver;
    6. input-file tests;
    7. documentation.
- **Upstream PRs:** possibly split in two:
  1. generic infrastructure (eigensolver, element action);
  2. the buckling analysis itself.
  This makes review easier for the maintainers.
- **AI declaration:** the PR description states that the code was developed with AI assistance, as recommended in `CONTRIBUTING.md`. Commits carry an attribution line. The contributing author reviews and takes responsibility for the code.
- **Validation before opening upstream:**
  - full CI run inside the fork `m-frey/4C`;
  - `clang-tidy` clean on the changed files;
  - no new warnings in the assertion build.
