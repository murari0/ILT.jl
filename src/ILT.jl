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

function lcurve(t, y, r, logα_min, logα_max, Nα, solveropts::Pair...; fn::KernelFunction = t1ir())
    alphas = exp10.(range(logα_min,logα_max,length=Nα))
    result = [ilt(t,y,r,solveropts...;α=a,fn) for a in alphas]
    residuals = [r[3] for r in result]
    norms = [norm(r[2]) for r in result]
    (residuals, norms, alphas)
end

end