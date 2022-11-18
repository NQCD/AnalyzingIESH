using DrWatson
@quickactivate "AnalyzingIESH"
using CairoMakie
using DataFrames
using CSV
using JamesPlots
using ColorSchemes
using DelimitedFiles: writedlm

include(srcdir("data.jl"))
include(srcdir("krr.jl"))
dash = [0.0, 3.2, 4.0]

results = collect_results!(datadir("fig5_6_7_8_desorption_probability"))
reference = CSV.read(datadir("physrevB2018", "fig6_gamma=1.0.csv"), DataFrame; header=[:t, :probability])

function plot_data()
    width = [4, 8, 16, 32, 64, 128, 256]

    fig = Figure(figure_padding=(1, 1, 1, 2), font=projectdir("fonts", "MinionPro-Capt.otf"))
    ax1 = MyAxis(fig[1,1], xscale=log2, xlabel="Band width /eV", ylabel="Desorption probability",
        limits=(nothing, nothing, nothing, nothing), xticks=width
    )

    hlines!(ax1, reference.probability[end] / 100, color=:black, linestyle=dash)

    dt = 1.0
    gamma = 1.0
    voltage = 0.0
    tmax = 200
    distribution = :wigner
    trajectories = 1000

    function plot_longtime!(ax, nstates, discretisation, color, linestyle, marker; kwargs...)
        longtimepopulation = Float64[]
        for w in width
            try
            dat = select_data(results; width=w, discretisation, nstates, dt, gamma, voltage, tmax, distribution, trajectories)
            push!(longtimepopulation, dat.probability[end])
            catch
            push!(longtimepopulation, 0.0)
            end
        end

        lambda = 0.005
        x_test = range(2, 8; length=200)
        y_pred = kernel_ridge_regression(with_lengthscale(MaternKernel(), 10), log2.(width), longtimepopulation, x_test, lambda)
        lines!(ax, 2 .^ x_test, y_pred; color, linestyle)

        scatter!(ax, width, longtimepopulation; color, marker, kwargs...)
        writedlm(projectdir("figure_data", "fig5", "$(discretisation)_$nstates.txt"), [width longtimepopulation])
    end

    cmap = JamesPlots.NICECOLORS
    colors = [ColorSchemes.get(cmap, (i-1)/2) for i=1:3]
    marker = :rect
    colors = COLORS[1:2:5]
    markersize = 4
    plot_longtime!(ax1, 100, :ShenviGaussLegendre, colors[1], :dash, :circle; markersize)
    plot_longtime!(ax1, 200, :ShenviGaussLegendre, colors[1], dash, :utriangle; markersize)
    plot_longtime!(ax1, 400, :ShenviGaussLegendre, colors[1], :solid, :rect; markersize)

    marker = :circle
    plot_longtime!(ax1, 100, :TrapezoidalRule, colors[2], :dash, :xcross; markersize)
    plot_longtime!(ax1, 200, :TrapezoidalRule, colors[2], dash, :dtriangle; markersize)
    plot_longtime!(ax1, 400, :TrapezoidalRule, colors[2], :solid, :diamond; markersize)

    state_labels = ["100 states", "200 states", "400 states"]

    elements_gl = [
        [LineElement(color=colors[1], linestyle=:dash), MarkerElement(;color=colors[1], marker=:circle, markersize)],
        [LineElement(color=colors[1], linestyle=dash), MarkerElement(;color=colors[1], marker=:utriangle, markersize)],
        [LineElement(color=colors[1], linestyle=:solid), MarkerElement(;color=colors[1], marker=:rect, markersize)],
    ]
    elements_t = [
        [LineElement(color=colors[2], linestyle=:dash), MarkerElement(;color=colors[2], marker=:xcross, markersize)],
        [LineElement(color=colors[2], linestyle=dash), MarkerElement(;color=colors[2], marker=:dtriangle, markersize)],
        [LineElement(color=colors[2], linestyle=:solid), MarkerElement(;color=colors[2], marker=:diamond, markersize)],
    ]

    elements = [elements_gl, elements_t]
    labels = [state_labels, state_labels]

    Legend(fig[1,1], elements, labels, ["Gauss-Legendre", "Trapezoid"];
        tellwidth=false, tellheight=false,
        valign=:bottom, halign=0.35,
        orientation=:horizontal,
        nbanks=3,
        titlefont=projectdir("fonts", "MinionPro-Capt.otf"),
        margin=(5,5,5,5),
    )

    fig
end

save_figure(plotsdir("fig5_convergence.png"), plot_data)
