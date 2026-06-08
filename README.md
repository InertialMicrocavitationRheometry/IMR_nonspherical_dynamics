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
|-- examples/                         Ready-to-run parameter permutations
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

See [IMRv2 Compatibility](#imrv2-compatibility) for setup details.

The included `mexhyp2f1.mexw64` and `mexhyp2f1.mexa64` files support Windows and
Linux. On another platform, or after a MATLAB/MEX compatibility change, rebuild
the MEX file from `common/make_hyp2f1/`.

## IMRv2 Compatibility

This repository does not include the IMRv2 radial bubble solver. Instead,
`common/f_call_IMRv2.m` is a wrapper that adds the IMRv2 forward-solver folder to
the MATLAB path and calls `f_imr_fd`.

The default expected folder layout is:

```text
Code/
|-- IMR_nonspherical_dynamics/
`-- IMRv2/
    `-- src/
        `-- forward_solver/
            `-- f_imr_fd.m
```

When MATLAB is run from the `IMR_nonspherical_dynamics` root,
`common/f_call_IMRv2.m` uses:

```matlab
addpath ../IMRv2/src/forward_solver/
```

If IMRv2 is stored somewhere else, either add the correct IMRv2 forward-solver
folder before running:

```matlab
addpath('C:\path\to\IMRv2\src\forward_solver')
```

or update the `addpath` line in `common/f_call_IMRv2.m`.

To check compatibility from MATLAB:

```matlab
cd path/to/IMR_nonspherical_dynamics
addpath common
which f_call_IMRv2
addpath ../IMRv2/src/forward_solver
which f_imr_fd
```

`which f_imr_fd` should return the IMRv2 solver path. If it returns
`f_imr_fd not found`, the radial solver path is not configured.

The wrapper currently assumes the IMRv2 `f_imr_fd` interface accepts name-value
inputs such as `radial`, `bubtherm`, `tvector`, `vapor`, `medtherm`,
`masstrans`, `method`, `stress`, `collapse`, `mu`, `g`, `alphax`, `surft`,
`r0`, `req`, `kappa`, `t8`, `rho8`, `pa`, `omega`, and `wave_type`, and returns
radial displacement, velocity, and acceleration in the positions used here:

```matlab
[t, R, Rd, ~, ~, ~, ~, Rdd] = f_imr_fd(...);
```

If a newer IMRv2 version changes option names, output order, or units, update
`common/f_call_IMRv2.m` before running LIC or ultrasound examples.

IMRv2 is required for:

- `s_basic_simulation.m`
- `examples/run_lic_*.m`
- `examples/run_ultrasound_*.m`

IMRv2 is not required for:

- `examples/run_free_*.m`, which use a prescribed constant-radius radial history.

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

## Example Permutations

The `examples/` folder contains scripts for free, laser-induced cavitation
(`lic`), and ultrasound-forced radial histories across viscous, elastic, and
viscoelastic material presets. Each case is available as both a single-mode and
multimode perturbation example.

For example:

```matlab
run('examples/run_lic_viscoelastic_single_mode.m')
```

The examples share `examples/run_nonspherical_example.m` and
`examples/example_config.m`, so changes to the common workflow can be made in
one place.

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
- `examples/`: simulation permutations built from the main driver workflow.

## Notes

- Run from the repository root unless you also update the relative paths in the
  scripts.
- Generated figures, PDFs, movies, and large `.mat` datasets should usually stay
  out of version control.
- `LICENSE` currently specifies the GNU General Public License v3.0.
