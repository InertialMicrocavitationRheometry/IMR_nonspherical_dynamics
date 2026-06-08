# IMR Nonspherical Dynamics

MATLAB code for simulating small-amplitude nonspherical perturbations of an
inertial microcavitation bubble in a soft material.

The current repository is organized around one runnable driver script,
`s_basic_simulation.m`, and a set of shared solver, numerical, and visualization
functions in `common/`. The driver computes a radial bubble history, evolves a
nonspherical spherical-harmonic perturbation, compares irrotational and
rotational perturbation models, and exports a strain-field snapshot PDF.

## Current Capabilities

- Radial bubble history through an external IMR radial solver interface.
- Nonspherical perturbation evolution for selected spherical harmonic modes.
- Rotational and irrotational perturbation model comparison.
- Finite-difference and trapezoidal-integration utilities for the transformed
  radial domain.
- Eulerian/current-domain visualization of strain, strain rate, and stress
  fields with `make_axisym_displacement_movie_all_fields`.
- Hypergeometric-function support through the included `mexhyp2f1` MEX binary
  and rebuild sources.

## Repository Layout

```text
IMR_nonspherical_dynamics/
|-- s_basic_simulation.m              Main runnable MATLAB driver
|-- common/                           Shared MATLAB functions and MEX support
|   |-- compute_rotational_perturbation_evolution.m
|   |-- f_*.m                         Numerical kernels and model terms
|   |-- make_axisym_displacement_movie_all_fields.m
|   |-- hyp2f1.m
|   |-- mexhyp2f1.mexw64              Windows MEX binary
|   |-- mexhyp2f1.mexa64              Linux MEX binary
|   `-- make_hyp2f1/                  Sources and build script for mexhyp2f1
|-- LICENSE
`-- README.md
```

The `data/` directory and `data.zip` are ignored by git. The current
`s_basic_simulation.m` script does not load files from `data/`, but the directory
may contain local experimental or simulation outputs used by related workflows.

## Requirements

- MATLAB, preferably a recent release.
- A MATLAB-supported C compiler if `mexhyp2f1` must be rebuilt.
- The external IMR radial solver repository available as a sibling directory:

```text
../IMRv2/src/forward_solver/
```

`common/f_call_IMRv2.m` adds that path and calls `f_imr_fd`. Run MATLAB from the
root of this repository so that relative path resolves correctly.

The included `mexhyp2f1.mexw64` and `mexhyp2f1.mexa64` files support Windows and
Linux. On another platform, or after a MATLAB/MEX compatibility change, rebuild
the MEX file from `common/make_hyp2f1/`.

## Running the Basic Simulation

From the repository root in MATLAB:

```matlab
addpath common
s_basic_simulation
```

The script:

1. Defines radial, material, and perturbation parameters.
2. Calls `f_call_IMRv2` to compute the spherical radial trajectory.
3. Calls `compute_rotational_perturbation_evolution` for irrotational and
   rotational perturbation evolution.
4. Plots the radius and perturbation amplitude histories.
5. Writes a snapshot PDF named `strain_test.pdf` in the current working
   directory.

The default driver uses mode `n = 8`, `xN = 256` radial grid points, and
`tsteps = 5000`.

## Rebuilding `mexhyp2f1`

Only rebuild this if MATLAB cannot load the included MEX binary.

```matlab
cd common/make_hyp2f1
make_hyp2f1
copyfile(['mexhyp2f1.' mexext], '..')
cd ../..
```

This requires a working MEX compiler configuration. In MATLAB, use `mex -setup`
if no compiler has been configured.

## Main Files

- `s_basic_simulation.m`: top-level example/driver.
- `common/f_call_IMRv2.m`: wrapper around the external IMR radial solver.
- `common/compute_rotational_perturbation_evolution.m`: main perturbation
  evolution solver.
- `common/f_*.m`: numerical stencils, transformed-domain mappings, matrix
  assembly, source terms, and constitutive/model coefficients.
- `common/make_axisym_displacement_movie_all_fields.m`: snapshot/movie
  visualization for strain, strain rate, and stress fields.
- `common/hyp2f1.m`: MATLAB wrapper for the compiled hypergeometric MEX routine.

## Notes

- Run from the repository root unless you also update the relative paths in the
  scripts.
- Generated figures, PDFs, movies, and large `.mat` datasets should usually stay
  out of version control.
- `LICENSE` currently specifies the GNU General Public License v3.0.

