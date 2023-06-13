# ILT.jl
A Julia implementation of a numerical Inverse Laplace Transform solver, intended for use in NMR relaxation analysis

## Dependencies
- LinearAlgebra (standard library)
- JuMP
- Ipopt

## Optional dependencies
- Plots, LaTeXStrings 
- [HSL](https://licences.stfc.ac.uk/products/Software/HSL)
  - Alternative, fast, linear solvers for Ipopt
  - Supports multi-threading, unlike the default solver MUMPS used by Ipopt.jl, providing significant performance boosts for ILT-FT computations

## Description
This package provides routines to perform a brute force inversion of the Laplace Transform on discrete real data.
