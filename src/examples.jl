using Plots

include("./ILT.jl")

t = [1E-4, 0.01, 0.05, 0.1, 0.25, 0.4, 1, 2, 4, 8, 16, 30, 45, 60];
y = [-0.9877, -0.9675, -0.8528, -0.7522, -0.4481, -0.217, 0.3236, 0.6201, 0.7275, 0.8054, 0.8883, 0.9491, 1, 0.993];

arange = exp10.(range(start=-5, stop=0, length=12));
Lcurve = [ilt(t, y, -2, 2, 200, α=a) for a in arange];
residuals = [L[3] for L in Lcurve];
norms = [norm(L[2]) for L in Lcurve];
scatter(residuals, norms, title="L-Curve", xaxis=("Residuals (||A*F(α) - y||^2)", :log), yaxis=("Norms (||F(α)||^2)", :log), legend=false)

(r, F, ) = ilt(t, y, -2, 2, 200, α=arange[10]);
plot((1 ./r[1:end-1], F[1:end-1], title="Inverse Laplace Transform", xaxis=("T1 (s)",:log), yaxis=("Intensity (arb)"), legend=false)
