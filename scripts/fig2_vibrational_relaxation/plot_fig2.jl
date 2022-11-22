using DrWatson
@quickactivate "MiaoSubotnikIESH"
using JamesPlots
using CairoMakie
using DataFrames
using CSV
using Unitful, UnitfulAtomic
using StatsBase
using LinearAlgebra: normalize
using Peaks
using ColorSchemes, Colors

include(srcdir("data.jl"))

function generate_headers()
    headers = readline(datadir("ref", "jcp2009_fig5.csv"))
    headers = split(headers, ",")[begin:2:end]
    collect(Iterators.flatten(zip(headers .* "x", headers .* "y")))
end

results = collect_results!(datadir("fig2_vibrational_relaxation"))
dash = [0.0, 3.2, 4.0]
Ms = [20, 40, 60, 80, 90]
nsteps = length(Ms)

colormap = ColorScheme(parse.(Colorant, ["#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977"]))
colors = [ColorSchemes.get(colormap, i/(nsteps)) for i in range(0, nsteps, length=nsteps)]

function plot_data!(ax, name, color, marker)
    data = readdlm(projectdir("figure_data", "fig2", "$(name).txt"))
    scatter!(ax, data[:,1], data[:,2]; color, marker, markersize=5)
end

function plot_data()
    fig = Figure(figure_padding=(1, 6, 1, 2), font=projectdir("fonts", "MinionPro-Capt.otf"), resolution=(JamesPlots.RESOLUTION[1],1.2*JamesPlots.RESOLUTION[2]))

    limits = (0, 500, nothing, nothing)
    xlabel = "Time /fs"
    ylabel = "Vibrational energy /eV"
    xticks = 0:100:500
    ax = MyAxis(fig[1,1]; limits, xlabel, ylabel, xticks)

    plot_data!(ax, "fortran_40_adiabatic", :black, :utriangle)
    for (i,n) in enumerate([20, 40, 60, 80, 90])
        plot_data!(ax, "fortran_$n", colors[i], :utriangle)
    end

    plot_data!(ax, "nqcd_40_adiabatic", :black, :circle)
    for (i,n) in enumerate([20, 40, 60, 80, 90])
        plot_data!(ax, "nqcd_$n", colors[i], :circle)
    end

    Colorbar(fig[0,1]; limits=(0,5), colormap=cgrad(colors, categorical=true), ticklabelsvisible=false,
        nsteps, ticklabelpad=-14, ticksvisible=false, label="Number of metal states", vertical=false, flip_vertical_label=true,
        labelpadding=0, minorticksvisible=true, minortickalign=1.0, minorticksize=16, minorticks=1:4)
    Label(fig[0,1], "20"; tellwidth=false, tellheight=false, halign=0.085, valign=0.40)
    Label(fig[0,1], "40"; tellwidth=false, tellheight=false, halign=0.29, valign=0.40)
    Label(fig[0,1], "60"; tellwidth=false, tellheight=false, halign=0.5, valign=0.40)
    Label(fig[0,1], "80"; tellwidth=false, tellheight=false, halign=0.71, valign=0.40)
    Label(fig[0,1], "90"; tellwidth=false, tellheight=false, halign=0.915, valign=0.40)

    return fig
end

save_figure(plotsdir("fig2_vibrational_relaxation.png"), plot_data, ncolors=5)
