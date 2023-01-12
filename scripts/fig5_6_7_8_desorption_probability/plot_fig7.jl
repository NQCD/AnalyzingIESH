using DrWatson
@quickactivate "AnalyzingIESH"
using CairoMakie
using DataFrames
using CSV
using JamesPlots
using ColorSchemes

dash = [0.0, 3.2, 4.0]

include(srcdir("data.jl"))
include(srcdir("krr.jl"))

function generate_headers(file)
    headers = readline(file)
    headers = split(headers, ",")[begin:2:end]
    collect(Iterators.flatten(zip(headers .* "x", headers .* "y")))
end

file = datadir("physrevB2018", "physrevB_fig4bottom.csv")
reference = CSV.read(file, DataFrame; skipto=3, header=generate_headers(file))
results = collect_results!(datadir("fig5_6_7_8_desorption_probability"))
colormap = JamesPlots.NICECOLORS

function fit_curve(x, y)
    i = findfirst(yi->isapprox(yi, 1.0, atol=0.01), y)
    if i === nothing
        i = length(x)
    end
    x_train = x[1:i]
    y_train = y[1:i]
    x_test = 0.0:0.01:x[i]
    lambda = 1e-4
    y_pred = kernel_ridge_regression(with_lengthscale(MaternKernel(), 1), x_train, y_train, x_test, lambda)
    clamp!(y_pred, 0, 1.0)
    return x_test[:], y_pred
end

function select_reference_data(mu)
    header = "mu=$(mu)eV"
    perm = sortperm(reference[!,header*"x"])
    x = reference[!,header*"x"][perm]
    y = reference[!,header*"y"][perm]

    i = findfirst(ismissing, x)
    if i === nothing
        x = Float64.(x)
        y = Float64.(y)
    else
        x = Float64.(x[1:i-1])
        y = Float64.(y[1:i-1])
    end
    return x, y
end

function plot_data()
    fig = Figure(resolution=(JamesPlots.RESOLUTION[1],1.2*JamesPlots.RESOLUTION[2]), figure_padding=(1, 5, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))
    ax1 = MyAxis(fig[1,1], xlabel="Γ /eV", ylabel="Desorption probability", limits=(0.0, 1.0, 0, nothing))

    mus = [0.0, 0.5, 1.0, 1.5, 2.0]
    for (color_index, mu) in enumerate(mus)
        if mu == 2.0
            x, y = select_reference_data(mu)
            lines!(ax1, x, y, color=colormap[color_index], linestyle=dash)
            scatter!(ax1, x, y, color=colormap[color_index], linestyle=dash, marker=:utriangle, markersize=5)
        else
            x, y = select_reference_data(mu)
            x_test, y_test = fit_curve(x, y)
            lines!(ax1, x_test, y_test, color=colormap[color_index], linestyle=dash)
            scatter!(ax1, x, y, color=colormap[color_index], linestyle=dash, marker=:utriangle, markersize=5)
        end
    end

    width = 64.0
    discretisation = :ShenviGaussLegendre
    nstates = 100
    dt = 1.0
    # gammas = [0.02, 0.1, 0.25, 0.4, 0.5, 0.8, 1.0]
    gammas = [0.02, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 0.8, 1.0]
    trajectories = 10000
    initial_state = :Adiabatic
    method = :AdiabaticIESH
    for (color_index, voltage) in enumerate(0.0:1.0:4.0)
        probabilities = Float64[]
        for gamma in gammas
            try
            dat = select_data(results; width, discretisation, nstates, dt, gamma, voltage, trajectories, initial_state, method)
            push!(probabilities, dat.probability[end])
            catch
            push!(probabilities, 0)
            end
        end

        if voltage == 4.0
            lines!(ax1, vcat(0.0, gammas), vcat(0.0, probabilities), color=colormap[color_index])
            scatter!(ax1, gammas, probabilities, color=colormap[color_index])
        else
            x_test, y_test = fit_curve(vcat(0.0, gammas), vcat(0.0, probabilities))
            lines!(ax1, x_test, y_test, color=colormap[color_index])
            scatter!(ax1, gammas, probabilities, color=colormap[color_index])
        end

        writedlm(projectdir("figure_data", "fig7", "chemicalpotential=$(voltage/2).txt"), [gammas probabilities])
    end

    Colorbar(fig[0,1]; limits=(-0.25, 2.25), colormap=cgrad(colormap[1:5], categorical=true), ticks=mus, nsteps=5, ticklabelpad=-14, ticksvisible=false, label="Chemical potential /eV", flip_vertical_label=true, vertical=false, labelpadding=6, minorticksvisible=true, minortickalign=1.0, minorticks=[0.25, 0.75, 1.25, 1.75], minorticksize=16)

    fig
end

save_figure(plotsdir("fig7_longtime_chemicalpotential.png"), plot_data, ncolors=6)
