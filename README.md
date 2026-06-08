# IMR Nonspherical Dynamics

`IMR_nonspherical_dynamics` is a MATLAB-based computational tool for modeling nonspherical bubble dynamics in soft materials, with applications to nonspherical Inertial Microcavitation Rheometry (IMR). The code simulates the coupled evolution of a radially oscillating cavitation bubble and small-amplitude nonspherical surface perturbations, enabling comparison between reduced-order models, full field simulations, and experimental observations from high-speed videography.

IMR relates the observed kinematics of a cavitation bubble to the pressure, stress, strain, and strain-rate fields generated in the surrounding material. This repository extends that framework to nonspherical bubble oscillations, allowing the stability, damping, stiffness, and mode evolution of nonspherical perturbations to be studied in viscous, elastic, and viscoelastic media.

`IMR_nonspherical_dynamics` is actively developed by Sawyer Remillard at Brown University.

For questions, contact Sawyer Remillard or Mauro Rodriguez, or request to join the IMR Slack workspace.

---

## Features

`IMR_nonspherical_dynamics` provides several tools for modeling and analyzing nonspherical inertial microcavitation.

### Nonspherical bubble dynamics simulation

The repository implements numerical solvers for nonspherical perturbations to a spherical cavitation bubble. The models are solved using centered finite differences, trapezoidal integration, and MATLAB-based time-integration routines.

The code supports studies of:

- Radial bubble dynamics
- Nonspherical surface perturbation growth and decay
- Spherical harmonic mode evolution
- Viscous, elastic, and viscoelastic material response
- Stress, strain, and strain-rate fields surrounding the bubble
- Comparison between potential-flow, reduced-order, and full-field formulations

Relevant background for the broader IMR framework includes:

- Estrada et al., "High Strain-rate Soft Material Characterization via Inertial Cavitation," *Journal of the Mechanics and Physics of Solids*, 2017.
- Warnez and Johnsen, "Numerical modeling of bubble dynamics in viscoelastic media with relaxation," *Physics of Fluids*, 2015.

### Nonspherical IMR analysis

The code is designed to support nonspherical extensions of IMR by connecting measurable bubble-shape evolution to the mechanical response of the surrounding material. This includes analysis of perturbation damping, perturbation stiffness, and mode-dependent material effects.

### MATLAB integration

The repository provides MATLAB scripts and functions for:

- Running nonspherical bubble simulations
- Controlling model parameters
- Generating figures and movies
- Visualizing bubble-shape evolution
- Plotting stress, strain, and strain-rate fields
- Comparing different model assumptions

### Experimental validation support

The repository includes tools intended to facilitate comparison between nonspherical bubble simulations and experimental observations. In particular, the code can be used to interpret high-speed imaging data in terms of spherical harmonic perturbation amplitudes when suitable bubble-surface tracking data are available.

---

## Getting Started

### Prerequisites

`IMR_nonspherical_dynamics` is implemented in MATLAB.

Recommended dependencies:

- MATLAB R2021a or newer
- Optimization Toolbox, recommended for fitting and parameter studies
- Image Processing Toolbox, optional, for workflows involving image or video analysis

Additional toolboxes may be required depending on which scripts or visualization routines are used.

---

## Installation

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/InertialMicrocavitationRheometry/IMR_nonspherical_dynamics.git
cd IMR_nonspherical_dynamics
```

Add the repository to your MATLAB path:

```matlab
addpath(genpath('IMR_nonspherical_dynamics'))
savepath
```

If the repository name or GitHub organization differs, replace the URL and folder name above with the correct repository path.

---

## Running an Example

To run a basic nonspherical bubble dynamics simulation, open MATLAB from the repository root and run one of the example or driver scripts.

For example:

```matlab
s_basic_simulation
```

or, if using a script from the examples directory:

```matlab
run('examples/run_example.m')
```

The exact example scripts may change as the repository develops. See the `examples/` directory for available cases.

---

## Repository Structure

A typical organization is:

```text
IMR_nonspherical_dynamics/
├── common/              # Shared utilities and helper functions
├── data/                # Small example datasets or parameter files
├── examples/            # Example simulation scripts
├── figures/             # Figure-generation scripts
├── src/                 # Core model and solver routines
├── tests/               # Verification and validation scripts
├── README.md            # Project overview
└── LICENSE              # License information
```

Large simulation outputs, movies, and raw experimental data should generally not be committed directly to the repository. Instead, they should be regenerated locally or stored externally.

---

## Documentation

Documentation is currently under development. For now, users should begin with the example scripts and the comments in the main solver and plotting routines.

Suggested starting points:

- `examples/`
- `src/`
- `common/`
- Main driver scripts in the repository root

---

## Citation

If you use `IMR_nonspherical_dynamics` in your research, please cite the associated publication once available.

For the broader IMR framework, please cite:

```bibtex
@article{Estrada_2017,
  title   = {High Strain-rate Soft Material Characterization via Inertial Cavitation},
  author  = {J. B. Estrada and C. Barajas and D. L. Henann and E. Johnsen and C. Franck},
  journal = {Journal of the Mechanics and Physics of Solids},
  year    = {2017},
  volume  = {112},
  pages   = {291--317},
  doi     = {10.1016/j.jmps.2017.12.006}
}
```

Additional relevant references include:

```bibtex
@article{Warnez_Johnsen_2015,
  title   = {Numerical modeling of bubble dynamics in viscoelastic media with relaxation},
  author  = {Warnez, M. T. and Johnsen, E.},
  journal = {Physics of Fluids},
  year    = {2015}
}
```

A project-specific citation should be added here once the nonspherical bubble dynamics manuscript is published or available as a preprint.

---

## License

`IMR_nonspherical_dynamics` is released under the GNU General Public License v3.0 unless otherwise specified in the repository.

See `LICENSE` for details.

---

## Acknowledgments

This work builds on the broader Inertial Microcavitation Rheometry framework developed by researchers at Brown University, the University of Michigan, the University of Wisconsin-Madison, Georgia Tech, and collaborators.

Development of the nonspherical bubble dynamics model has been supported by research efforts in soft-material characterization, cavitation dynamics, and high strain-rate mechanics.

---

## Contact

For issues, bug reports, or feature requests, please open an issue on GitHub.

For research questions, contact:

- Sawyer Remillard, Brown University
- Mauro Rodriguez, Brown University

Users interested in broader IMR development may also request to join the IMR Slack workspace.
