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

<details>
  <summary> Mathematical background </summary>
  
  The Laplace transform of a function $f(t)$ is defined as $$\mathcal{L} \lbrace f \rbrace (s) = \int_{0}^{\infty} f(t)e^{-st}dt$$ where, in general, the transformed function denoted as $\bar{f}(s)$ is a complex-valued function of the complex variable $s$.
  In an NMR $T_2$ relaxation experiment, the signal intensity of a signal with time constant $T_2$ decays as $I_0(t)e^{-t/T_2}$. In a sample with multiple peaks, each with its own relaxation behaviour, the total integrated signal decays as $$I(t) = \sum_n I_{0n}e^{-t/T_{2n}}$$ Relabelling signal components by their relaxation rates $R = 1/T_2$, we have $$I(t) = \sum_R I_R^0e^{-Rt}$$ which functionally resembles a discrete version of the Laplace transform.
</details>

