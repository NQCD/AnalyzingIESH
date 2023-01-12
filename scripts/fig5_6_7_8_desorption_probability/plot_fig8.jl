using DrWatson
@quickactivate "ErpenbeckThossModel"
using CairoMakie
using DataFrames
using CSV
using JamesPlots
using ColorSchemes, Colors

include(srcdir("data.jl"))
include(srcdir("krr.jl"))

dash = [0.0, 3.2, 4.0]

function generate_headers(file)
    headers = readline(file)
    headers = split(headers, ",")[begin:2:end]
    collect(Iterators.flatten(zip(headers .* "x", headers .* "y")))
end

file = datadir("physrevB2018", "physrevB_fig4top.csv")
reference = CSV.read(file, DataFrame; skipto=3, header=generate_headers(file))
results = collect_results!(datadir("fig5_6_7_8_desorption_probability"))
colormap = ColorScheme(parse.(Colorant, ["#B7E6A5", "#7CCBA2", "#46AEA0", "#089099", "#00718B", "#045275", "#003147"]))
colors = [ColorSchemes.get(colormap, i/5) for i=0:4]

function fit_curve(x, y)
    i = findfirst(yi->isapprox(yi, 1.0, atol=0.01), y)
    if i === nothing
        i = length(x)
    end
    x_train = x
    y_train = y
    x_test = 0.0:0.01:x[end]
    lambda = 0.001
    y_test = kernel_ridge_regression(with_lengthscale(MaternKernel(), 1), x_train, y_train, x_test, lambda)
    clamp!(y_test, 0, 1.0)
    for i in eachindex(x_test, y_test)
        if x_test[i] > 2.0
            y_test[i] = 1.0
        end
    end
    return x_test, y_test
end

function select_reference_data(gamma)
    header = "gamma=$(gamma)eV"
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
    fig = Figure(resolution=(JamesPlots.RESOLUTION[1],1.2*JamesPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))
    ax1 = MyAxis(fig[1,1], xlabel="Chemical potential /eV", ylabel="Desorption probability", limits=(0.0, 3.0, 0, nothing))

    gammas = [0.02, 0.1, 0.25, 0.5, 1.0]
    for (color_index, gamma) in enumerate(gammas)
        x, y = select_reference_data(gamma)
        x_test, y_test = fit_curve(x, y)
        lines!(ax1, x_test, y_test, color=colors[color_index], linestyle=dash)
        scatter!(ax1, x, y, color=colors[color_index], linestyle=dash, marker=:utriangle, markersize=5)
    end

    width = 64.0
    discretisation = :ShenviGaussLegendre
    nstates = 100
    dt = 1.0
    voltages = [0.0, 0.5, 1.0, 1.5, 1.75, 2.0, 2.5, 3.0, 3.25, 3.5, 4.0, 5.0, 6.0]
    trajectories = 10000
    initial_state = :Adiabatic
    method = :AdiabaticIESH
    for (color_index, gamma) in enumerate(gammas)
        probabilities = Float64[]
        for voltage in voltages
            try
            dat = select_data(results; width, discretisation, nstates, dt, gamma, voltage, trajectories, initial_state, method)
            push!(probabilities, dat.probability[end])
            catch
            push!(probabilities, 0)
            end
        end
        x_test, y_test = fit_curve(voltages./2, probabilities)
        lines!(ax1, x_test, y_test, color=colors[color_index])
        scatter!(ax1, voltages./2, probabilities, color=colors[color_index])

        writedlm(projectdir("figure_data", "fig8", "gamma=$gamma.txt"), [voltages./2 probabilities])
    end

    vlines!(ax1, 2.33, linestyle=dash, color=:black)
    ticks = (1:5, string.(gammas))
    Colorbar(fig[0,1]; limits=(0.5, 5.5), colormap=cgrad(colors, categorical=true), ticks=ticks, nsteps=5, ticklabelpad=-14, ticksvisible=false, label="Γ /eV", flip_vertical_label=true, vertical=false, labelpadding=6, minorticksvisible=true, minortickalign=1.0, minorticks=[1.5, 2.5, 3.5, 4.5], minorticksize=16)

    fig
end

save_figure(plotsdir("fig8_longtime_chemicalpotential.png"), plot_data, ncolors=6)
