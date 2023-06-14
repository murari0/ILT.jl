module ILT

using LinearAlgebra
using JuMP, Ipopt

include("./kernels.jl")
include("./utils.jl")

export
    ilt,
    lcurve,
    read_bruker_p2d,
    read_2d!

"""
    ilt(t, y, logr_min, logr_max, N[, solveropts...]; kwargs...)
    ilt(t, y, r, N[, solveropts...]; kwargs...)

Computes the inverse Laplace transform of `y(t)`

Sets up an optimisation problem to invert the Laplace transform of the function defined by the arguments `t` and `y`. The vector `r` (as given by the user or generated as a logarithmically spaced range using `logr_min`, `logr_max`, `N` parameters) is the set of points in r-space at which the inverse is computed. The optimisation problem is solved using Ipopt.

`solveropts` is a list of optional arguments specified as `"key"=>"value"` pairs that are passed to Ipopt. For example, the key `"linear_solver"` tells Ipopt what linear solver to use - specifying `"linear_solver"=>"ma27"` will instruct Ipopt to use the STFC's HSL_MA27 solver, if available.

The output tuple `(r, F, residual)` contains the values `F` of the inverse at each point in `r`. Both vectors have length `N+1` - an additional point, `0.0` is pushed to the start of the set of relaxation rates to account for any residual y-offset in the data. The vectors `F[2:end]` and `r[2:end]` contain the solution at finite relaxation times/rates.

# Keyword arguments
- `α = 1.0`: Tikhonov regularisation parameter
- `fn = ILT.t1ir()`: Type of recovery/relaxation behaviour, defined in the `KernelFunction` datatype

See also [`ILT.KernelFunction`](@ref), [`lcurve`](@ref).
"""
function ilt(t, y, logr_min, logr_max, N, solveropts::Pair...; α = 1.0, fn::KernelFunction = t1ir())
    r = exp10.(range(logr_min,logr_max,length=N))
    ilt(t,y,r,solveropts...;α,fn)
end

function ilt(t, y, r, solveropts::Pair...; α = 1.0, fn::KernelFunction = t1ir())
    N = length(r)
    r = [0.0; r]

    A = fn.(t,r')
    Areg = vcat(A, α*I(N+1))
    yreg = vcat(y, zeros(N+1,1))

    model = Model(Ipopt.Optimizer)
    if !isempty(solveropts)
        set_attributes(model, solveropts...);
    end
    set_silent(model)
    @variable(model,x[1:N+1])
    @constraint(model, [i=2:N+1], x[i] >=0);
    @objective(model, Min, sum((Areg*x-yreg).^2))
    optimize!(model)

    F = value.(x)[1:N+1]
    residual = norm(A*F-y)
    return (r,F,residual)
end

"""
    lcurve(t, y, logr_min, logr_max, N, logα_min, logα_max, Nα[, solveropts...]; fn = ILT.t1ir())
    lcurve(t, y, r, logα_min, logα_max, Nα[, solveropts...]; fn = ILT.t1ir())

Computes the L-Curve for a given function `y(t)`, to choose the optimum value of `α` for Tikhonov regularisation.

Creates a vector of `α` values using the `logα_min`, `logα_max`, `Nα` parameters given, computes the inverse Laplace transform of `y` for each value of `α`, and returns the residual and solution norm for each solution, which define the L-Curve used to choose the optimum value of `α`. The vector `r` at which the inverse function is computed is either used as provided by the user or generated as a logarithmically spaced vector using the `logr_min`, `logr_max`, `N` values given.

`solveropts` is a list of optional arguments specified as `"key"=>"value"` pairs that are passed to Ipopt. For example, the key `"linear_solver"` tells Ipopt what linear solver to use - specifying `"linear_solver"=>"ma27"` will instruct Ipopt to use the STFC's HSL_MA27 solver, if available.

# Keyword arguments
- `fn = ILT.t1ir()`: Type of recovery/relaxation behaviour, defined in the `KernelFunction` datatype

See also: [`ILT.KernelFunction`](@ref), [`ilt`](@ref).
"""
function lcurve(t, y, logr_min, logr_max, N, logα_min, logα_max, Nα, solveropts::Pair...; fn::KernelFunction = t1ir())
    r = exp10.(range(logr_min,logr_max,length=N))
    lcurve(t,y,r,logα_min,logα_max,Nαsolveropts...;fn)
end

function lcurve(t, y, r, logα_min, logα_max, Nα, solveropts::Pair...; fn::KernelFunction = t1ir())
    alphas = exp10.(range(logα_min,logα_max,length=Nα))
    result = [ilt(t,y,r,solveropts...;α=a,fn) for a in alphas]
    residuals = [r[3] for r in result]
    norms = [norm(r[2]) for r in result]
    (residuals, norms, alphas)
end

end
