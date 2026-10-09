# A Fully Discrete, Mass-Conservative Adaptive Discontinuous Spectral-Element Method for the One-Dimensional Richards Equation in Mixed Form

This repository provides the Julia code for reproducing the numerical experiments and figures presented in the following manuscript:

```bibtex
@unpublished{MontoyaSpiteriAdaptiveRichards2026,
  title={A Fully Discrete, Mass-Conservative Adaptive Discontinuous
         Spectral-Element Method for the One-Dimensional Richards Equation
         in Mixed Form},
  author={Montoya, Tristan and Spiteri, Raymond J.},
  year={2026},
  note={Submitted to SIAM Journal on Scientific Computing.},
}
```

When using or adapting this code or the associated solver in [HydroTrixi.jl](https://github.com/uofs-simlab/HydroTrixi.jl) for your research, please cite the manuscript above.

## Abstract

We present and analyze a spatially and temporally adaptive high-order method for the one-dimensional Richards equation describing the vertical transport of water in unsaturated soil, with an emphasis on resolving sharp wetting fronts while satisfying a fully discrete water mass balance. The spatial discretization uses an arbitrary-order local discontinuous Galerkin spectral-element formulation with collocated Legendre-Gauss-Lobatto quadrature, and the time integration is based on a high-order Rosenbrock-Wanner method for stiff problems. We focus primarily on the mixed formulation of the Richards equation, which advances water content as the conserved variable while evaluating the nonlinear diffusion operator in terms of pressure head; the resulting fully discrete scheme satisfies the water mass balance in exact arithmetic, independent of spatial resolution and time step size, including under spatial and temporal adaptivity. In contrast, the pressure-head formulation advances pressure head through a nonlinear capacity function and satisfies the corresponding conservation identity only at the semi-discrete level. Numerical experiments examine convergence through spatial and temporal refinement studies, with mass conservation assessed under spatial and temporal adaptivity using the recently proposed solver flux output method (SFOM). The SFOM results confirm machine-precision mass conservation for the mixed formulation, whereas the pressure-head formulation exhibits mass balance errors that grow as the temporal error tolerance is relaxed. Moreover, in the cases considered, the mixed formulation yields more accurate solutions than the pressure-head formulation for a prescribed time step size. The combination of high-order accuracy, spatial and temporal adaptivity, and fully discrete mass conservation makes the proposed mixed-form discontinuous spectral-element method a promising approach for simulating unsaturated flow in next-generation land models.

## Installation

The scripts use syntax compatible with Julia 1.10 and newer; use Julia 1.13.1 to reproduce the recorded environment in `code/Manifest.toml`.

Clone this repository and instantiate the environment from `code/`:

```sh
git clone https://github.com/tristanmontoya/paper-2026-adaptive-richards.git
cd paper-2026-adaptive-richards/code
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The environment fixes HydroTrixi.jl to tag [`v0.1.1`](https://github.com/uofs-simlab/HydroTrixi.jl/tree/v0.1.1), including the solver, diagnostics, and visualization methods used by these scripts. For repeated experiments and plotting, start a persistent REPL from `code/`:

```sh
OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=.
```

## Run the experiments

From the Julia REPL in `code/`, run both suites and retain the returned directory paths:

```julia
include("convergence/run_richards_convergence.jl")
convergence = RichardsConvergence.run_convergence()

include("conservation/run_richards_conservation.jl")
conservation = RichardsConservation.run_conservation()
```

This performs all 60 convergence solves and 24 conservation solves. Numerical tables are saved in `results/richards_convergence/RUN_ID/data/` and `results/richards_conservation/RUN_ID/data/`, where `RUN_ID` is assigned automatically when each runner starts as a UTC timestamp, formatted as `yyyymmddTHHMMSSsssZ` (for example, `20261008T213045123Z`). Conservation snapshots are saved in `RUN_ID/snapshots/`. The runner scripts do not generate figures.

### Individual studies and command-line use

From `code/`, run a selected study:

```sh
OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=. \
  convergence/run_richards_convergence.jl DN space 4

OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=. \
  conservation/run_richards_conservation.jl haverkamp mixed water_content
```

For the convergence studies (Section 5.2 in the paper), the first argument (`DD` or `DN`) specifies the boundary conditions, the second argument specifies the type of study (`space` or `time`), and the third argument represents the polynomial degree of the spatial discretization (`3` or `4` for the spatial convergence studies, and `3` for the temporal convergence studies). Both the mixed and pressure-head formulations are run for each convergence study.

For the conservation studies (Section 5.3 in the paper), the first argument (`haverkamp` or `new_mexico`) specifies the infiltration problem, the second argument specifies the formulation (`mixed` or `pressure_head`), and the third argument specifies the solution transfer variable (`water_content` or `pressure_head`). Each conservation study runs at relative tolerances `1e-5`, `1e-7`, and `1e-9`.

In the REPL, use `Symbol`s such as `:space` and `:mixed` for named study options, and pass a vector of study specifications as `NamedTuple`s:

```julia
include("convergence/run_richards_convergence.jl")
convergence = RichardsConvergence.run_convergence(
    [(; bc = :dirichlet_neumann, kind = :space, N = 4)])

include("conservation/run_richards_conservation.jl")
conservation = RichardsConservation.run_conservation(
    [(; benchmark = :haverkamp, formulation = :mixed, transfer = :pressure_head)])
```

More specifications can be added to the vector to select multiple studies. Calling either runner without a study vector runs its full suite. Either runner accepts the keyword argument `output_directory = "NEW_RUN_DIRECTORY"` to choose a user-specified directory rather than the default timestamped directory.

## Generate the paper figures

Saved numerical tables, conservation snapshots, and reference PDFs are included in the [convergence results](results/richards_convergence/final/) and [conservation results](results/richards_conservation/final/) directories. Use `final` in place of `RUN_ID` in the plotting commands below to regenerate the figures without rerunning the simulations.

From the Julia REPL in `code/`, generate figures from the available saved tables. Missing solution snapshots are reported and skipped. Replace `RUN_ID` in each path with the saved run directory name for that suite; convergence and conservation runs generally have different `RUN_ID` values:

```julia
include("convergence/plot_richards_convergence.jl")
convergence_pdfs = RichardsConvergencePlots.plot_convergence(
    "../results/richards_convergence/RUN_ID")

include("conservation/plot_richards_conservation.jl")
conservation_pdfs = RichardsConservationPlots.plot_conservation(
    "../results/richards_conservation/RUN_ID")
```

Each call returns the PDF paths and prints its new `results/richards_convergence/RUN_ID/figures/TIMESTAMP/` or `results/richards_conservation/RUN_ID/figures/TIMESTAMP/` directory. Repeated plotting preserves earlier figures and numerical tables as it generates a new timestamp. Command-line plotting also runs from `code/`:

```sh
julia --project=. convergence/plot_richards_convergence.jl \
  ../results/richards_convergence/RUN_ID

julia --project=. conservation/plot_richards_conservation.jl \
  ../results/richards_conservation/RUN_ID
```

Both plotters accept additional directories so that they can plot results that span multiple `RUN_ID`s; figures are written within the first directory listed. An error is returned if the same `.dat` filename appears in multiple listed directories. Conservation figures use whichever of the three supported tolerance tables are available for each study. The convergence plotter reads every `.dat` file in each input directory's `data/` subdirectory and generates a corresponding PDF, replacing the `.dat` extension with `.pdf` in the file name. The input tables and output panels selected for the manuscript are listed below.

### Figure 1: Spatial convergence

Input tables from `results/richards_convergence/RUN_ID/data/`:

- `richards_manufactured_solution_dirichlet_dirichlet_space_N3.dat`
- `richards_manufactured_solution_dirichlet_dirichlet_space_N4.dat`
- `richards_manufactured_solution_dirichlet_neumann_space_N3.dat`
- `richards_manufactured_solution_dirichlet_neumann_space_N4.dat`

Manuscript output panels:

- **(a)** `richards_manufactured_solution_dirichlet_dirichlet_space_N3.pdf`
- **(b)** `richards_manufactured_solution_dirichlet_dirichlet_space_N4.pdf`
- **(c)** `richards_manufactured_solution_dirichlet_neumann_space_N3.pdf`
- **(d)** `richards_manufactured_solution_dirichlet_neumann_space_N4.pdf`

### Figure 2: Temporal convergence

Input tables from `results/richards_convergence/RUN_ID/data/`: 

- `richards_manufactured_solution_dirichlet_dirichlet_time_N3.dat`
- `richards_manufactured_solution_dirichlet_neumann_time_N3.dat`

Manuscript output panels:
  
- **(a)** `richards_manufactured_solution_dirichlet_dirichlet_time_N3.pdf`
- **(b)** `richards_manufactured_solution_dirichlet_neumann_time_N3.pdf`

### Figure 3: Infiltration profiles and meshes

Snapshot tables from `results/richards_conservation/RUN_ID/snapshots/`:

- `richards_celia_haverkamp_mixed_rtol1e-7_half_solution_mesh.dat`
- `richards_celia_haverkamp_mixed_rtol1e-7_full_solution_mesh.dat`
- `richards_celia_new_mexico_mixed_rtol1e-7_half_solution_mesh.dat`
- `richards_celia_new_mexico_mixed_rtol1e-7_full_solution_mesh.dat`

The following conservation history tables must also be present in `results/richards_conservation/RUN_ID/data/`:

- `richards_celia_haverkamp_mixed_rtol1e-7.dat`
- `richards_celia_new_mexico_mixed_rtol1e-7.dat`

Manuscript output panels:

- **(a)** `richards_celia_haverkamp_mixed_rtol1e-7_half_solution_mesh.pdf`
- **(b)** `richards_celia_haverkamp_mixed_rtol1e-7_full_solution_mesh.pdf`
- **(c)** `richards_celia_new_mexico_mixed_rtol1e-7_half_solution_mesh.pdf`
- **(d)** `richards_celia_new_mexico_mixed_rtol1e-7_full_solution_mesh.pdf`

### Figure 4: Accepted time step sizes

Input tables from `results/richards_conservation/RUN_ID/data/`:

- `richards_celia_haverkamp_mixed_rtol1e-5.dat`
- `richards_celia_haverkamp_mixed_rtol1e-7.dat`
- `richards_celia_haverkamp_mixed_rtol1e-9.dat`
- `richards_celia_new_mexico_mixed_rtol1e-5.dat`
- `richards_celia_new_mexico_mixed_rtol1e-7.dat`
- `richards_celia_new_mexico_mixed_rtol1e-9.dat`

Manuscript output panels:

- **(a)** `richards_celia_haverkamp_mixed_time_steps.pdf`
- **(b)** `richards_celia_new_mexico_mixed_time_steps.pdf`

### Figure 5: Absolute mass bias

Input tables from `results/richards_conservation/RUN_ID/data/`:

- `richards_celia_haverkamp_mixed_rtol1e-5.dat`
- `richards_celia_haverkamp_mixed_rtol1e-7.dat`
- `richards_celia_haverkamp_mixed_rtol1e-9.dat`
- `richards_celia_haverkamp_mixed_pressure_head_transfer_rtol1e-5.dat`
- `richards_celia_haverkamp_mixed_pressure_head_transfer_rtol1e-7.dat`
- `richards_celia_haverkamp_mixed_pressure_head_transfer_rtol1e-9.dat`
- `richards_celia_haverkamp_pressure_head_pressure_head_transfer_rtol1e-5.dat`
- `richards_celia_haverkamp_pressure_head_pressure_head_transfer_rtol1e-7.dat`
- `richards_celia_haverkamp_pressure_head_pressure_head_transfer_rtol1e-9.dat`
- `richards_celia_haverkamp_pressure_head_water_content_transfer_rtol1e-5.dat`
- `richards_celia_haverkamp_pressure_head_water_content_transfer_rtol1e-7.dat`
- `richards_celia_haverkamp_pressure_head_water_content_transfer_rtol1e-9.dat`
- `richards_celia_new_mexico_mixed_rtol1e-5.dat`
- `richards_celia_new_mexico_mixed_rtol1e-7.dat`
- `richards_celia_new_mexico_mixed_rtol1e-9.dat`
- `richards_celia_new_mexico_mixed_pressure_head_transfer_rtol1e-5.dat`
- `richards_celia_new_mexico_mixed_pressure_head_transfer_rtol1e-7.dat`
- `richards_celia_new_mexico_mixed_pressure_head_transfer_rtol1e-9.dat`
- `richards_celia_new_mexico_pressure_head_pressure_head_transfer_rtol1e-5.dat`
- `richards_celia_new_mexico_pressure_head_pressure_head_transfer_rtol1e-7.dat`
- `richards_celia_new_mexico_pressure_head_pressure_head_transfer_rtol1e-9.dat`
- `richards_celia_new_mexico_pressure_head_water_content_transfer_rtol1e-5.dat`
- `richards_celia_new_mexico_pressure_head_water_content_transfer_rtol1e-7.dat`
- `richards_celia_new_mexico_pressure_head_water_content_transfer_rtol1e-9.dat`

Manuscript output panels:

- **(a)** `richards_celia_haverkamp_mixed_mass_bias.pdf`
- **(b)** `richards_celia_haverkamp_pressure_head_mass_bias.pdf`
- **(c)** `richards_celia_new_mexico_mixed_mass_bias.pdf`
- **(d)** `richards_celia_new_mexico_pressure_head_mass_bias.pdf`

## License

This repository is released under the [MIT License](LICENSE).
