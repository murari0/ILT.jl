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

[Background on the Laplace Transform](docs/overview.pdf)

The package exports the following symbols:
```julia
ilt(t, y, logr_min, logr_max, N, solveropts::Pair...; α = 1.0, fn::ILT.KernelFunction = ILT.t1ir())
ilt(t, y, r, solveropts::Pair...; α = 1.0, fn::ILT.KernelFunction = ILT.t1ir())
lcurve(t, y, r, logα_min, logα_max, Nα, solveropts::Pair...; fn::ILT.KernelFunction = ILT.t1ir())
read_bruker_p2d(filepath, f1dim, f2dim, f1length; stride_len = 1)
read_2d!(data::Array{T,2}, filepath::String) where {T<:Number}
```
The two `ilt()` methods accept 
  - a vector of delays `t` 
  - the corresponding intensities `y` 
  - optional options `solveropts` that are passed through directly to the solver (see [Ipopt's documentation](https://www.coin-or.org/Bonmin/option_pages/options_list_ipopt.html) for a complete description of available options) 
  - the Tikhonov regularisation factor `α` 
  - a Kernel function `fn` that determines what functional form to use for the transform (inversion recovery - default, saturation recovery, $T_2$ relaxation) 
  - either a vector of relaxation rates `r`, or values for the minimum (`logr_min`), maximum (`logr_max`) and length (`N`) used to generate such a vector, over which the inverse is calculated

Both methods return a tuple `(r, F, residual)`, where `r` is the input relaxation rate vector with an extra rate of `0.0` pushed to the start to represent y-offset, `F` is the computed ILT at each point in `r` and `residual` is the residual rms error from the data values in `y`.

The `lcurve()` method accepts, in addition to the parameters passed through to `ilt()`, a minimum (`logα_min`), maximum (`logα_max`) and length (`Nα`) to generate a vector of `α` values at which to compute the ILT. The method returns a tuple of vectors `(residuals, norms, alphas)` containing the vector of `α` values, the residual rms error, and the norm of the solution obtained at each `α` value. These can be used to plot the L-Curve to determine the optimum value of `α` to use.
