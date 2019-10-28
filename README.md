# CD_Solver

This will be a object-oriented PDE-solver, that allows one to solve a coupled nonlinear cross-diffusion system using the finite element software package ngsolve (https://ngsolve.org/). Using implicit Euler in time, finite element in space and Newton to overcoming nonlinearities.


## Installation

* Install ngsolve by following the instruciton on https://ngsolve.org/.
* Install matplotlib (https://matplotlib.org/) and pandas (https://pandas.pydata.org/)
* Clone and change directory to CD_package:

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


## Usage

See solve_example.py and plot_example.py for examples to solve and plot an experiment.
