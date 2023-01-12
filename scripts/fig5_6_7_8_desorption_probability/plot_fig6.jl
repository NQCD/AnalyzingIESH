using DrWatson
@quickactivate "AnalyzingIESH"
using CairoMakie
using DataFrames
using CSV
using JamesPlots
using ColorSchemes

include(srcdir("data.jl"))
dash = [0.0, 3.2, 4.0]

function generate_headers()
    headers = readline(datadir("physrevB2018", "physrevB_fig6.csv"))
    headers = split(headers, ",")[begin:2:end]
    collect(Iterators.flatten(zip(headers .* "x", headers .* "y")))
end

results = collect_results!(datadir("fig5_6_7_8_desorption_probability"))
reference = CSV.read(datadir("physrevB2018", "physrevB_fig6.csv"), DataFrame; skipto=3, header=generate_headers())
# cmeresults = collect_results!(datadir("cme_nonadaptive"))

missingcolumns = ["0.02eV_1.0V", "0.02eV_0.0V", "0.1eV_0.0V", "0.25eV_0.0V"]
colormap = JamesPlots.NICECOLORS

function plot_data()
    fig = Figure(resolution=(JamesPlots.RESOLUTION[1],2.0*JamesPlots.RESOLUTION[2]), figure_padding=(1, 6, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))
    axes = [MyAxis(fig[i,1], xlabel="Time /fs", limits=(0, 200, 0, 1.2)) for i=1:5]

    hidexdecorations!.(axes[1:end-1], ticks=false, minorticks=false)

    gammas = [0.02, 0.1, 0.25, 0.5, 1.0]
    volts = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    i = 0
    for g in gammas
        i += 1
        for (color_index, v) in enumerate(volts)
            header = "$(g)eV_$(v)V"
            header in missingcolumns && continue
            lines!(axes[i], reference[!,header*"x"], reference[!,header*"y"], color=ColorSchemes.get(colormap, color_index/length(volts)), linestyle=dash)
        end
    end


    width = 64
    dt = 1.0
    nstates = 100
    discretisation = :ShenviGaussLegendre
    trajectories = 10000
    shift_band = true
    initial_state = :Adiabatic
    method = :AdiabaticIESH

    i = 0
    for gamma in gammas
        i += 1
        for (color_index, voltage) in enumerate(volts)
            try
            dat = select_data(results; gamma, voltage, width, dt, nstates, discretisation, trajectories, shift_band, initial_state, method)
            t = range(0, dat.tmax, length=length(dat.probability))
            lines!(axes[i], t, dat.probability; color=ColorSchemes.get(colormap, color_index/length(volts)), label="$voltage V")
            writedlm(projectdir("figure_data", "fig6", "iesh_gamma=$(gamma)_chemicalpotential=$(voltage/2).txt"), [t dat.probability])
            catch;end
        end
    end

    dt = 1.0
    trajectories = 10000
    method = :CME

    i = 0
    for gamma in gammas[1:2]
        i += 1
        for (color_index, voltage) in enumerate(volts)
            try
            dat = select_data(results; gamma, voltage, dt, trajectories, method)
            t = range(0, dat.tmax, length=length(dat.probability))
            lines!(axes[i], t, dat.probability; color=ColorSchemes.get(colormap, color_index/length(volts)), label="$voltage V", linestyle=:dot)
            writedlm(projectdir("figure_data", "fig6", "cme_gamma=$(gamma)_chemicalpotential=$(voltage/2).txt"), [t dat.probability])
            catch;end
        end
    end

    for (i, gamma) in enumerate(gammas)
        Label(fig[i,1], "Γ = $gamma eV", tellwidth=false, tellheight=false, halign=:left, valign=:top, padding=(8, 8, 8, 8))
    end

    Label(fig[:,0], "Desorption probability", rotation=π/2, padding=(0, 5, 0, 0))

    Colorbar(fig[0,1]; limits=(-0.25, 3.25), colormap=cgrad(colormap[1:7], categorical=true), ticks=0:0.5:3,
        nsteps=7, ticklabelpad=-14, ticksvisible=false, label="Chemical potential /eV", vertical=false,
        flip_vertical_label=true, labelpadding=6, minorticksvisible=true, minortickalign=1.0, minorticks=0.25:0.5:2.75, minorticksize=16)

    fig
end

save_figure(plotsdir("fig6_desorption_probabilities.png"), plot_data, ncolors=6)
