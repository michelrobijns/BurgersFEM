# Solving Burgers' Equation in 1D using the Finite Element Method

Work in progress!

This is a C++ implementation of a finite element method (FEM) solver for Burgers' equation in 1D. There also a [version of this program in Python 3](https://github.com/michelrobijns/pyBurgersFEM).

## Features

* Solves Burgers' equation in 1D
* Uses piecewise-linear basis functions
* Uses the backward Euler method (or implicit Euler method) to advance the solution through time
* Uses Newton's method to solve the system of nonlinear equations
* Supports Dirichlet boundary conditions and periodic boundary conditions
* Supports uniformly spaced meshes

## Dependencies

* Apple Accelerate framework (yes, the code only compiles on macOS for now...)
* OpenMP (comes with most compilers, so not really a dependency)

## Usage

Compile using `make`, run `main` and then `plot.sh` to plot the solution.
