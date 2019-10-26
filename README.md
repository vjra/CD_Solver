# CD_Solver

This will be a object-oriented PDE-solver, that allows one to solve a coupled nonlinear cross-diffusion system using the finite element software package ngsolve (https://ngsolve.org/). Using implicit Euler in time, finite element in space and Newton to overcoming nonlinearities.

## Installation

Install the ngsolve

```sh
$ git clone https://github.com/vjra/CD_Solver.git
$ cd CD_package
```
## Introduction
Generates a experiment object to solve parabolic-parabolic system for $u=(\rho,c) \in \R^2$

$$u_t = \operatorname{div}(A(u) \nabla u) + f(u)$$

for t > 0 on a bounded domain Omega,
with source term
$f(u) = (0,-u[1]+u[0]^alpha)$
and parameters
with $eps >= 0$ and $alpha >0$.

The experiment object builds on the ngsolve library and contains:
* the geometry,
* the mesh,
* the weak formulation of the problem,
* Test and Trial functions,
* simulation parameters,
* finite element space,
* initial data for the problem,
* solution of the computation as a ngsolve gridfunction.

Methods
-------
set_file_structure(args*): Changes file and folder structure of the experiment.

meshgeneration_n_spaces(): Generates Mesh, trial and test function and finite element spaces.

set_initial_data(args*): Sets initial data.

weak_formulation(args*): Sets weak formulation for the problem.

run_experiment(): Solves problem using solver.py.

plotseries_and_more(args*): Plots every time step and exports it to vtk. Creates L^\infty plot vs time.

plot_and_export(args*): Plots specific time steps and exports it to vtk. Creates surfaces plots and exports them to eps.

## Usage

See solve_example.py and plot_example.py for examples to solve and plot an experiment.
