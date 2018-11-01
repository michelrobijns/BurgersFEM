# Solving Burgers' Equation in 1D using the Finite Element Method

This is a C++ implementation of a finite element method (FEM) solver for Burgers' equation in 1D. There also a [version of this program in Python 3](https://github.com/michelrobijns/pyBurgersFEM).

## Usage

Open a terminal and `cd` to the top-level directory, compile using `make`, run `main` and then run `plot.sh` to plot the solution.

Parameters like the mesh spacing and the time step are defined in `main.cpp`. The forcing function, initial condition, and boundary values are defined as functions inside `main.cpp`. In other words, to change the problem all you need to do is modify `main.cpp` and recompile.

## Verification

To verify the correctness of the implementation, two test problems (one with Dirichlet and one with periodic boundary conditions) are included for which the exact reference solutions are known:

* `test_dirichlet.cpp` (compile with `make test_dirichlet`)
* `test_periodic.cpp` (compile with `make test_periodic`)

These test problems were generated using the method of manufactured solutions. After running the executables, solutions can be plotted by running `test_plot.sh`.

## Features

* Solves Burgers' equation in 1D
* Uses piecewise-linear basis functions
* Uses the backward Euler method to advance the solution through time
* Uses Newton's method to solve the system of nonlinear equations
* Supports Dirichlet boundary conditions and periodic boundary conditions
* Supports uniformly and nonuniformly spaced meshes
* All linear system solves are of O(n)
* Assembly of the residual and the Jacobian is parallelized using OpenMP

## Dependencies

* Apple Accelerate framework (yes, the code will only compile on macOS for now...)
* OpenMP (comes with most compilers, so not really a dependency)
* Python 3 with matplotlib to plot the solution
