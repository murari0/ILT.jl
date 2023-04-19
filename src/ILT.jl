using LinearAlgebra
using JuMP, Ipopt

include("./kernels.jl")

function ilt(t, y, logr_min, logr_max, N; α = 1, fn::KernelFunction = t1ir())

    r = exp10.(range(logr_min, logr_max, length=N))
    # y-offset
    push!(r, 0)

    # Populate kernel matrix
    A = fn.(t,r')

    # Add rows for regularisation
    R = α*I(N+1)
    R[end] = 0;
    AR = vcat(A, R)
    yR = vcat(y, zeros(N+1,1))

    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model,x[1:N+1])
    @constraint(model, [i=1:N], x[i] >=0);
    @objective(model, Min, sum((AR*x-yR).^2))
    optimize!(model)
    F = value.(x)[1:N+1]
    residual = norm(A*F-y)
    return (r,F,residual)
end
