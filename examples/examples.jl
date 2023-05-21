using ILT
using Plots
using LaTeXStrings

function demo_1dilt()
    t = [1E-4, 0.01, 0.05, 0.1, 0.25, 0.4, 1, 2, 4, 8, 16, 30, 45, 60]
    y = [-0.9877, -0.9675, -0.8528, -0.7522, -0.4481, -0.217, 0.3236, 0.6201, 0.7275, 0.8054, 0.8883, 0.9491, 1, 0.993]
    r = exp10.(range(-2,2,length=200))

    (residuals, norms, alphas) = lcurve(t,y,r,-5,1,12)
    ann = [text(string(i)*", "*string(round(alphas[i],digits=2)),:bottom,:left,pointsize=8) for i in 1:length(alphas)]
    p1 = scatter(residuals, norms, series_annotations=ann, title="L-Curve", xaxis=(L"\textrm{Residuals\quad}(||A*F(α) - y||^2)", :log), yaxis=(L"\textrm{Solution norms\quad}(||F(α)||^2)", :log), legend=false)

    (r_res, F, ) = ilt(t,y,r,α=alphas[9])
    t1s = reverse(1./r_res[2:end])
    a = round(alphas[9],digits=2)
    p2 = plot(t1s, reverse(F[2:end]), title="Inverse Laplace Transform, α = $a", xaxis=(L"T_1\ (\mathrm{s})",:log), yaxis=(L"\textrm{Intensity (arb)}"), legend=false)

    plot(p1, p2, layout=(1,2))
end

function demo_iltft(data::Array{Float64,2}, t::Vector{Float64}, f::Vector{Float64}, α::Float64; multithreaded = false, lsolver = "")
    if size(data) != (length(t), length(f))
        throw(DimensionMismatch("Length of delay list and/or frequency axis does not match data"))
    end

    r = exp10.(range(-2,2,length=200))
    vdata = eachcol(data)
    S = Array{Float64,2}(undef,200,size(data,2))
    vS = eachcol(S)

    if multithreaded
        if isempty(lsolver) || lsolver == "mumps"
            throw(ErrorException("MUMPS does not support multi-threading. Please specify an alternative solver"))
        end
        Threads.@threads for i in 1:length(vS)
            vS[i] = ilt(t,vdata[i],r,"linear_solver"=>lsolver;α)[2][2:end]
        end
    else
        vS .= [ilt(t,y,r,"linear_solver"=>lsolver;α)[2][2:end] for y in vdata]
    end

    t1 = reverse(1 ./r)
    contour(reverse(f),t1,reverse(S),xaxis=("Resonance offset δ (ppm)",:flip),yaxis=(L"T_1 (\mathrm{s})",:log),levels=20,guidefontfamily="Computer Modern")
    return S
end