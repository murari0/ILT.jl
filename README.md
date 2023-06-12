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
  In an NMR $T_2$ relaxation experiment, the signal intensity of a signal with time constant $T_2$ decays as $I_0(t)e^{-t/T_2}$. In a sample with multiple peaks, each with its own relaxation behaviour, the total integrated signal decays as $$I(t) = \sum_n I_{0n}e^{-t/T_{2n}}$$ sampled at discrete time points $\lbrace t_m \rbrace$. This resembles a discrete version of the Laplace transform with $T_2 = 1/s$, where both the original and transformed functions are real-valued functions of real variables.
  The amplitudes of the relaxation components that exist in the signal, $\lbrace I_{0n} \rbrace$, can be extracted by inverting the Laplace transform. Unlike the Fourier transform, however, the Laplace transform is not unitary and the inverse cannot be easily calculatied[^1]. Writing the transformation above in matrix form with the following redefinitions
  ```math
  \mathbf{b} = \left[ I(t_0) I_(t_1) \dots I(t_m) \right]\intercal
  \mathbf{x} = \left[ I_{01} I_{02} \dots I_{0n} \right]\intercal
  \mathbf{A} = \begin{bmatrix} e^{-t_1s_1} & \cdots & e^{-t_1s_n} \\ \vdots & \ddots & \vdots \\ e^{-t_ms_1} & \cdots & e^{-t_ms_n} \end{bmatrix}
  ```
  we have a system of linear equations $\mathbf{Ax} - \mathbf{b} = 0$ which can be solved as an optimisation problem with constraints ${x_i} \ge 0$. The same approach can be used for $T_1$ recovery experiments by changing the functional form of the elements of the transformation matrix $A$ to recovery curves. 
</details>

[1]: Numerical algorithms to compute the inverse exist (<https://link.springer.com/article/10.1007/s11075-012-9625-3>), but generally require the function to be sampled at a large number of arbitrary points, which is not compatible with experimentally acquired datasets.
