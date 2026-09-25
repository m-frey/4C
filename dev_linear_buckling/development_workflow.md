# Development workflow

How to build, test and commit on the `dev-linear-buckling` branch so that the result can go upstream without rework. Coding conventions are in section 7 of [`implementation_plan.md`](implementation_plan.md).

## 1. Branches

| Branch | Purpose |
|---|---|
| `dev-linear-buckling` (fork `m-frey/4C`) | Development branch: planning folder `dev_linear_buckling/` + implementation |
| `main` (fork) | Kept in sync with upstream `4C-multiphysics/4C:main` |
| later: clean PR branch(es) | Created from upstream `main` for the upstream pull request(s), containing only the implementation commits, **without** `dev_linear_buckling/` |

Rules:
- Merge upstream `main` into `dev-linear-buckling` regularly (merge commit, no rebase of pushed history).
- Keep planning changes (files in `dev_linear_buckling/`) and code changes in **separate commits**. The implementation commits can then be cherry-picked onto a clean branch for the upstream PR, and the planning material never enters upstream history.

## 2. Environment

4C is built inside the CI dependency image (Trilinos, SuiteSparse, …). This is the same environment CI uses, so local results match CI:

```bash
docker pull ghcr.io/4c-multiphysics/4c-dependencies-ubuntu24.04:eeb1c043   # tag as in .github/workflows/buildtest.yml
docker run -it --rm -v "$PWD":/home/user/4C -w /home/user/4C \
  --env OMPI_ALLOW_RUN_AS_ROOT=1 --env OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  ghcr.io/4c-multiphysics/4c-dependencies-ubuntu24.04:eeb1c043 bash
```

Notes:
- The image tag changes when upstream updates its dependencies. Always take it from `.github/workflows/buildtest.yml` after merging upstream `main`.
- **Claude Code cloud sessions:** Docker works in the container, but pulling from `ghcr.io` also needs the host `pkg-containers.githubusercontent.com` allowed under *Network access* in the environment settings. Docker Hub (`dolfinx/dolfinx`) already works.
- Resources seen in the cloud container: 4 cores, 15 GB RAM. A full 4C build takes on the order of an hour there; incremental builds take minutes. The container is temporary, so push work regularly.

## 3. Build

Same preset as the CI assertion build:

```bash
mkdir -p ../4C-build && cd ../4C-build
cmake /home/user/4C --fresh --preset=docker_assertions
cmake --build . --target full -- -j "$(nproc)"
```

For faster iterations, build only what is needed: the `4C` executable and the unit-test targets of the touched modules. List them with `cmake --build . --target help | grep -i test`.

## 4. Tests

```bash
cd ../4C-build
# new buckling tests
ctest -R buckling --output-on-failure
# neighbouring tests that must stay green
ctest -R "solid_ele|rve3d_periodic_bcs|structure_new" -j "$(nproc)" --output-on-failure
# unit tests of touched modules (named unittests_<module>..., list them with: ctest -N | grep unittests)
ctest -R "unittests_(solid_ele|core_linalg)" --output-on-failure
```

- Input-file tests live in `tests/input_files/*.4C.yaml` and are registered in `tests/list_of_tests.cmake` with `four_c_test(TEST_FILE <name>.4C.yaml NP <n>)`. Unregistered files fail the `check-test-files` hook.
- The reference values and tolerances for the new tests are in [`verification/`](verification/). Start with the loose verification tolerances. Once a result is verified, freeze the 4C value in the `RESULT DESCRIPTION` with a tight regression tolerance.
- Full regression: open a PR **inside the fork** (`dev-linear-buckling` → fork `main`). This triggers `buildtest` and `checkcode` on GitHub runners, provided Actions are enabled for the fork.

## 5. Code checks before every commit

```bash
./utilities/set_up_dev_env.sh          # once: creates utilities/python-venv and installs the git hooks
source utilities/python-venv/bin/activate
pre-commit run --files <changed files>  # or: pre-commit run --all-files
```

- The hooks format C++, CMake, YAML and Python, insert license headers, and check file names, header guards, includes, include cycles, typos, non-ASCII characters in code, test registration and the commit message.
- CI (`checkcode.yml`) runs `pre-commit run --all-files` on pull requests to `main`, so everything in the branch must pass. That includes this folder: Markdown is typo-checked, and Python files are only allowed in the locations the `check-python-files` hook accepts. That is why the scripts are embedded in Markdown here.
- `clang-tidy` also runs in CI with the preset `docker_codeclimate`. Keep changed files free of new findings.

## 6. Commits

Enforced by `utilities/code_checks/commit-msg`:
- subject ≤ 50 characters, capitalized, imperative mood, no trailing period;
- a blank line between subject and body;
- more than one word in the subject.

Style of the existing history, e.g. `Add differential stiffness to solid elements`, `Add subspace eigensolver`, `Add linear buckling input section`.

One logical change per commit; each commit builds and passes its tests. Suggested sequence: section 7.5 of the plan.

AI assistance: declare it in the PR description, as recommended in `CONTRIBUTING.md`, and add an attribution line to AI-generated commits. The author reviews and takes responsibility for every change.

## 7. Regenerating reference values

See [`verification/reference_values.md`](verification/reference_values.md). Analytical values need Python with `numpy`/`scipy`. The FEniCSx values need Docker with `dolfinx/dolfinx:v0.9.0`, about 1 min for the published mesh and about 6 min for `80×8×8` on 4 cores.

Rebuild the verification PDF after editing the `.tex`:

```bash
cd dev_linear_buckling/verification
pdflatex buckling_verification_cases.tex && pdflatex buckling_verification_cases.tex
rm -f *.aux *.log *.out *.toc
```
