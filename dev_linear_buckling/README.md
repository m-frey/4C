# Linear buckling analysis for 4C — development folder

This folder holds the planning, verification and workflow material for adding an **Abaqus `*BUCKLE`-equivalent linear (eigenvalue) buckling analysis** for 3D solid elements to 4C. Periodic RVE boundary conditions are supported as well. The analysis applies a load pattern and returns the buckling eigenvalues (load factors) and buckling modes.

> This folder is development material for the fork branch `dev-linear-buckling`. It is **not** meant for upstream 4C. The upstream pull request(s) will be built from the implementation commits only (see [`development_workflow.md`](development_workflow.md), section 1).

## Status

Planning and verification design are complete. Implementation has not started.

- [x] Feasibility study of the existing code (structure_new, solid_ele, constraint framework, linear solvers)
- [x] Design choices fixed (D1–D5, identical results to Abaqus / FEniCSx)
- [x] Verification cases A–D defined, reference values computed
- [x] Code-style and mergeability rules collected from the repository
- [ ] **M1** element differential stiffness, eigensolver, driver, input; tests C (FEniCSx code-to-code) and B (cantilever)
- [ ] **M2** VTU/CSV output, Dirichlet load pattern; tests A1/A2 (plate) and D (periodic lattice)
- [ ] **M3** nonlinear preload, shift, MPI variants, user documentation
- [ ] **M4** EAS (then F-bar) differential stiffness; test A1 with `HEX8 eas_full`
- [ ] **M5** imperfection seeding for follow-up nonlinear runs
- [ ] Full CI in the fork, then upstream PR(s)

## Contents

| File | What it is |
|---|---|
| [`implementation_plan.md`](implementation_plan.md) | The plan: design choices D1–D5, formulation, architecture (files to add/change), input format, eigensolver, milestones, test mapping, **code style and mergeability rules (section 7)** |
| [`development_workflow.md`](development_workflow.md) | Branch policy, build environment (CI Docker image), build and test commands, pre-commit, commit-message rules |
| [`verification/buckling_verification_cases.pdf`](verification/buckling_verification_cases.pdf) (`.tex` source alongside) | Verification document: one chapter per test case (setup, derivation) and one chapter with all **hard results to match** and acceptance tolerances |
| [`verification/reference_values.md`](verification/reference_values.md) | Provenance of every reference number, embedded scripts and the FEniCSx reproduction procedure with raw output |

## How to continue (read this first)

For a developer or an AI coding session picking this up:

1. Read the **design choices D1–D5** at the top of `implementation_plan.md`. They are fixed and determine the results. Everything else is an implementation choice and may be changed if D1–D5 still hold and cases A–D pass.
2. Read **section 7 of the plan** (code style and mergeability) and `development_workflow.md`. New code must look like the neighbouring 4C code and pass all pre-commit hooks and CI checks.
3. Set up the environment and do a first full build (`development_workflow.md`, sections 2–3).
4. Start **M1, step 1**: the element action `struct_calc_differential_stiffness` in `src/solid_ele` with its unit test. Then continue through M1 in the listed order. The first input-file test to get passing is case C, which must match FEniCSx to 1e-4.
5. After each step: build, run the new and neighbouring tests, run `pre-commit`, commit (one logical change, 4C commit-message rules), push to `dev-linear-buckling`, and tick the checklist above.

## Verification cases at a glance

| Case | Structure | Reference | First eigenvalue to match |
|---|---|---|---|
| A1 | Simply supported square plate, 100 × 100 × 1 mm, uniaxial compression | analytical (Kirchhoff + Mindlin) | 75.877241 N/mm (≤ 1 %) |
| A2 | Simply supported plate 300 × 100 × 1 mm | analytical | 8 closely spaced modes in exact order |
| B | Cantilever column 100 mm, 3 × 2 mm section | analytical (Euler + Timoshenko) | 103.604771 N (≤ 1 %) |
| C | FEniCSx "Linear Buckling Analysis of a 3D solid", 50 × 5 × 5 `HEX27` | published benchmark, rerun | 0.167963169 (≤ 1e-4) |
| D | Periodic square-lattice cell 1×1 and 2×2 | exact identities + frame estimate | 10.0145 / 2.5036 (± 15 %) |

## Prerequisites and open points

- **Build image:** the CI dependency image lives on `ghcr.io`. In Claude Code cloud sessions the host `pkg-containers.githubusercontent.com` must be allowed in the environment's network settings to pull it.
- **Fork CI:** GitHub Actions must be enabled on the fork for the full regression run.
- **Existing limitations that are inherited:** periodic RVE conditions in 4C are serial-only and penalty-only (`src/constraint_framework/4C_constraint_framework_submodelevaluator_mpc.cpp`).
- **Decide later:** one or two upstream PRs. Two would split the generic infrastructure (eigensolver + element action) from the buckling analysis itself.
