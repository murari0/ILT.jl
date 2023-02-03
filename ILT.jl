using LinearAlgebra
using JuMP, Ipopt 
using Plots

function ilt(t, y, logr_min, logr_max, N; α=1)

    r = exp10.(range(logr_min, logr_max, length=N))
    # zero rate constant for y-offset in input data
    push!(r, 0)

    # kernel function
    a(i, j) = 1-2*exp(-t[i]*r[j])

    # convert a(i,j) to  matrix
    A = [a(i, j) for i = eachindex(t), j = eachindex(r)]

    # add rows for regularisation, do not regularise offset
    R = α*I(N+1)
    R[end]=0;
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

t = [1E-4, 0.01, 0.05, 0.1, 0.25, 0.4, 1, 2, 4, 8, 16, 30, 45, 60]
y = [-0.9877, -0.9675, -0.8528, -0.7522, -0.4481, -0.217, 0.3236, 0.6201, 0.7275, 0.8054, 0.8883, 0.9491, 1, 0.993]