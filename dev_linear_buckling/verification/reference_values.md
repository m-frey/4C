# Reference values: provenance and reproduction

This file documents where every number in [`buckling_verification_cases.pdf`](buckling_verification_cases.pdf) comes from and how to regenerate it. The scripts are embedded as code blocks rather than committed as `.py` files. The repository's `check-python-files` hook only allows Python files in `utilities/four_c_python/`, `tests/input_files/` and module directories.

Units: cases A, B and D use mm, N, MPa. Case C uses the dimensionless units of the FEniCSx tour.

## Summary of hard values

| Case | Quantity | Value | Source |
|---|---|---|---|
| A1 | λ₁ … λ₆ (N/mm), Mindlin | 75.877241, 118.458034, 210.295972, 302.996602, 341.190800, 355.100634 | analytical, script 1 |
| A2 | λ₁ … λ₈ (N/mm), Mindlin | 75.877241, 82.314032, 89.064318, 97.411294, 118.458034, 144.519114, 175.197311, 210.822928 | analytical, script 1 |
| B | λ₁ … λ₆ (N), Timoshenko | 103.604771, 233.037442, 930.569806, 2087.883909, 2574.572274, 5016.053142 | analytical, script 1 |
| C | λ₁ … λ₈, FEniCSx 50×5×5 | 0.167963169, 0.496955439, 0.987888956, 1.500093336, 1.642490988, 2.455332621, 3.429242204, 4.392907156 | FEniCSx run, section 3 |
| C | λ₁ … λ₈, FEniCSx 80×8×8 | 0.168078306, 0.496866498, 0.987687419, 1.500692451, 1.640922961, 2.451987909, 3.421467271, 4.392598990 | FEniCSx run, section 3 |
| D | λ₁ D1 / D2, frame estimate (rigid joints) | 10.0145 / 2.5036 | closed form, script 3 |

## 1. Cases A and B (analytical)

Run with `python3` (needs only the standard library):

```python
from math import pi, sqrt

E, nu = 210000.0, 0.3
G = E / (2 * (1 + nu))


def plate_modes(a, b, t, count, kappa=5 / 6):
    """Simply supported plate under uniaxial N_x: Kirchhoff and Mindlin values."""
    D = E * t**3 / (12 * (1 - nu**2))
    rows = []
    for m in range(1, 15):
        for n in range(1, 5):
            k = (m * b / a + n**2 * a / (m * b)) ** 2
            n_kirchhoff = k * pi**2 * D / b**2
            wave2 = (m * pi / a) ** 2 + (n * pi / b) ** 2
            n_mindlin = n_kirchhoff / (1 + D * wave2 / (kappa * G * t))
            rows.append((n_kirchhoff, n_mindlin, m, n, k))
    return sorted(rows)[:count]


def cantilever_modes(length, width, height, count):
    """Clamped-free column: Euler and Engesser (Timoshenko) values, both axes."""
    area = width * height
    kappa = 10 * (1 + nu) / (12 + 11 * nu)  # Cowper
    rows = []
    for axis, inertia in (("weak", width * height**3 / 12), ("strong", height * width**3 / 12)):
        for n in range(1, 5):
            p_euler = (2 * n - 1) ** 2 * pi**2 * E * inertia / (4 * length**2)
            p_timo = p_euler / (1 + p_euler / (kappa * G * area))
            rows.append((p_euler, p_timo, axis, n))
    return sorted(rows)[:count]


for row in plate_modes(100.0, 100.0, 1.0, 6):
    print("A1", row)
for row in plate_modes(300.0, 100.0, 1.0, 8):
    print("A2", row)
for row in cantilever_modes(100.0, 3.0, 2.0, 6):
    print("B", row)
```

## 2. Case C, beam-theory comparison values (analytical)

Needs `numpy` and `scipy`:

```python
import numpy as np
from scipy.optimize import brentq

E, nu = 1.0e3, 0.3
G = E / (2 * (1 + nu))
length, width, height = 1.0, 0.01, 0.03  # width along y, height along z
area = width * height
kappa = 10 * (1 + nu) / (12 + 11 * nu)

# clamped-pinned column: tan(alpha) = alpha
roots = [brentq(lambda a: np.tan(a) - a, n * np.pi + 1e-9, (n + 0.5) * np.pi - 1e-9) for n in range(1, 9)]
rows = []
for axis, inertia in (("weak (y)", height * width**3 / 12), ("strong (z)", width * height**3 / 12)):
    for n, alpha in enumerate(roots, 1):
        p_euler = alpha**2 * E * inertia / length**2
        p_timo = p_euler / (1 + p_euler / (kappa * G * area))
        rows.append((p_euler / area, p_timo / area, axis, n, alpha / np.pi))
for row in sorted(rows)[:8]:
    print(row)
```

The exact roots are α/π = 1.430297, 2.459024, 3.470890, 4.477409, … The published tour quotes 2.4509 for the second root, which is a typo.

## 3. Case C, FEniCSx reference run (numerical)

Source: J. Bleyer, *Numerical Tours of Computational Mechanics with FEniCSx*, tour "Linear Buckling Analysis of a 3D solid", <https://github.com/bleyerj/comet-fenicsx>, directory `tours/eigenvalue_problems/buckling_3d_solid`, commit `5d06010c` (4 March 2025). The tour is licensed CC BY-SA 4.0 and is therefore **not** copied into this repository. Reproduce it as follows:

1. Get the tour and the container:
   ```bash
   git clone https://github.com/bleyerj/comet-fenicsx.git
   cd comet-fenicsx/tours/eigenvalue_problems/buckling_3d_solid
   docker pull dolfinx/dolfinx:v0.9.0
   ```
2. Apply these minimal changes. Nothing that affects the numerics is changed.
   - `eigenvalue_solver.py`: replace `vr.vector, vi.vector` with `vr.x.petsc_vec, vi.x.petsc_vec` (API change in dolfinx 0.9).
   - `buckling_3d_solid.py`:
     - read `Nx, Ny, Nz` from the command line instead of the fixed `50, 5, 5`;
     - set `N_eig = 12`;
     - limit the tour's own comparison loop to 9 entries, since its `alpha` array has 9 entries;
     - remove the `pyvista` plotting block;
     - append:
     ```python
     print("mesh", Nx, Ny, Nz, "dofs", V.dofmap.index_map.size_global * 3)
     for i in range(min(10, len(eigval))):
         vec = eigvec_r[i].x.array.reshape(-1, 3)
         norm = np.linalg.norm(vec)
         parts = [np.linalg.norm(vec[:, k]) / norm for k in (0, 1, 2)]
         print(f"MODE {i+1} {eigval[i].real:.9f} ux={parts[0]:.3f} uy={parts[1]:.3f} uz={parts[2]:.3f}")
     ```
3. Run:
   ```bash
   docker run --rm -v "$PWD":/w -w /w dolfinx/dolfinx:v0.9.0 python3 buckling_3d_solid.py 50 5 5
   ```

Solver settings (unchanged from the tour): SLEPc Krylov–Schur, shift-and-invert with shift 0.01, tolerance 1e-12, Q2 Lagrange hexahedra (27 nodes, i.e. 4C `SOLID HEX27`).

Raw output. The component columns are the fractions of the mode norm per displacement component:

```text
mesh 50 5 5 dofs 36663                      (published mesh, primary reference)
MODE 1 0.167963169 ux=0.011 uy=1.000 uz=0.000
MODE 2 0.496955439 ux=0.019 uy=1.000 uz=0.000
MODE 3 0.987888956 ux=0.027 uy=1.000 uz=0.001
MODE 4 1.500093336 ux=0.033 uy=0.000 uz=0.999
MODE 5 1.642490988 ux=0.034 uy=0.999 uz=0.001
MODE 6 2.455332621 ux=0.042 uy=0.999 uz=0.002
MODE 7 3.429242204 ux=0.049 uy=0.999 uz=0.003
MODE 8 4.392907156 ux=0.056 uy=0.000 uz=0.998
MODE 9 4.556151336 ux=0.057 uy=0.998 uz=0.004
MODE 10 5.840050132 ux=0.064 uy=0.998 uz=0.005

mesh 100 5 5 dofs 72963
MODE 1 0.168117833 ux=0.011 uy=1.000 uz=0.000
MODE 2 0.496878159 ux=0.019 uy=1.000 uz=0.000
MODE 3 0.987769722 ux=0.027 uy=1.000 uz=0.001
MODE 4 1.500819381 ux=0.033 uy=0.000 uz=0.999
MODE 5 1.640838561 ux=0.034 uy=0.999 uz=0.001
MODE 6 2.451837817 ux=0.042 uy=0.999 uz=0.002
MODE 7 3.420793724 ux=0.049 uy=0.999 uz=0.003
MODE 8 4.392497519 ux=0.056 uy=0.000 uz=0.998
MODE 9 4.541413933 ux=0.057 uy=0.998 uz=0.004
MODE 10 5.813782178 ux=0.064 uy=0.998 uz=0.005

mesh 80 8 8 dofs 139587                     (converged reference)
MODE 1 0.168078306 ux=0.011 uy=1.000 uz=0.000
MODE 2 0.496866498 ux=0.018 uy=1.000 uz=0.000
MODE 3 0.987687419 ux=0.026 uy=1.000 uz=0.001
MODE 4 1.500692451 ux=0.032 uy=0.000 uz=0.999
MODE 5 1.640922961 ux=0.033 uy=0.999 uz=0.001
MODE 6 2.451987909 ux=0.041 uy=0.999 uz=0.002
MODE 7 3.421467271 ux=0.048 uy=0.999 uz=0.003
MODE 8 4.392598990 ux=0.054 uy=0.000 uz=0.999
MODE 9 4.542564030 ux=0.055 uy=0.998 uz=0.003
MODE 10 5.816138067 ux=0.062 uy=0.998 uz=0.004
```

Notes:
- A `100×10×10` run needs more than 16 GB of memory with the tour's direct solver and was not used.
- A `50×10×10` run failed in SLEPc (fewer converged pairs than requested) and was not investigated.
- The FE eigenvalues are not monotone under refinement: the discrete pre-stress enters `K_σ`, so they are not bounds.

## 4. Case D, frame-theory estimate (closed form)

The x-struts carry the load. The struts have plane-strain bending stiffness `E w³ / (12 (1 − ν²))` per unit depth and free length `ℓ`.

```python
from math import pi

nu, strut_width, strain_ref = 0.3, 1.0, 1.0e-3
for free_length in (20.0, 19.0):  # point joints: L; rigid joints of size w: L - w
    eps_d1 = 4 * pi**2 * strut_width**2 / (12 * (1 - nu**2) * free_length**2)  # clamped-clamped
    eps_d2 = eps_d1 / 4  # alternating sway (clamped-guided)
    print(free_length, eps_d1 / strain_ref, eps_d2 / strain_ref)
```

Output: `20.0 → 9.0381 / 2.2595`, `19.0 → 10.0145 / 2.5036`.

A periodic Euler–Bernoulli frame FE model (16 elements per strut, rigid joint zones) reproduced these values to four significant digits: 10.0147 and 2.5036. In the 2×2 spectrum it also contained the 1×1 eigenvalues, as required by identity I1.
