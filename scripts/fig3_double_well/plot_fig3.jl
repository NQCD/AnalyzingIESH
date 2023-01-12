using DrWatson
@quickactivate "AnalyzingIESH"
using JamesPlots
using CairoMakie
using DataFrames
using CSV

results = collect_results!(datadir("fig3_double_well"))

dash = [0.0, 3.2, 4.0]
include(srcdir("data.jl"))

function plot_data()
    fig = Figure(figure_padding=(1, 10, 1, 1), resolution=(JamesPlots.COLUMN_WIDTH, 2*JamesPlots.RESOLUTION[2]), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))

    ax1 = MyAxis(fig[2,1], width=Relative(0.80), height=Relative(0.8), halign=0.80, valign=0.5, limits=(0, 75, 0, 1))
    ax2 = MyAxis(fig[4,1], width=Relative(0.80), height=Relative(0.8), halign=0.80, valign=0.5, limits=(0, 75, 2.0, 4.0), xlabelpadding=-10)

    ax3 = MyAxis(fig[2:3,1], xlabel="ωt", ylabel="Impurity population", limits=(0, 10000, 0.0, 1.0), xticks=0:2500:10000, ylabelpadding=2, yticks=0.20:0.20:1)
    ax4 = MyAxis(fig[4:5,1], xlabel="ωt", ylabel="Kinetic energy /kT", limits=(0, 10000, 0, 5), xticks=0:2500:10000, ylabelpadding=2)


    hidexdecorations!(ax1; ticks=false, minorticks=false, ticklabels=false)
    hidexdecorations!(ax2; ticks=false, minorticks=false, ticklabels=false)
    hidexdecorations!(ax3; ticks=false, minorticks=false)

    hideydecorations!(ax1; ticks=false, minorticks=false, ticklabels=false, label=true)
    hideydecorations!(ax2; ticks=false, minorticks=false, ticklabels=false, label=true)
    hideydecorations!(ax3; ticks=false, minorticks=false, ticklabels=false, label=false)
    hideydecorations!(ax4; ticks=false, minorticks=false, ticklabels=false, label=false)

    function plot_data!(dat, label, ax1, ax2, color, filename)
        t = range(0, dat.tmax, length=size(dat.output,1))
        population = dat.output[:,1]
        lines!(ax1, t, population; color)
        writedlm(projectdir("figure_data", "fig3", "population_$filename.txt"), [t population])

        kinetic = dat.output[:,2]
        kT = 9.5e-4
        lines!(ax2, t, kinetic ./ kT; label, color)
        writedlm(projectdir("figure_data", "fig3", "kinetic_$filename.txt"), [t kinetic ./ kT])
    end

    reference = DataFrame(CSV.File(datadir("ref", "reference.csv"), header=["t", "population"]))
    lines!(ax1, reference.t, reference.population, color=:black)
    reference = DataFrame(CSV.File(datadir("ref", "kinetic.csv"), header=["t", "kinetic"]))
    lines!(ax2, reference.t, reference.kinetic, label="Reference", color=:black)

    interval = 1
    reference = DataFrame(CSV.File(datadir("ref", "longtimepopulation.csv"), header=["t", "population"]))
    lines!(ax3, reference.t[begin:interval:end], reference.population[begin:interval:end], color=:black)
    reference = DataFrame(CSV.File(datadir("ref", "longtimekinetic.csv"), header=["t", "kinetic"]))
    lines!(ax4, reference.t[begin:interval:end], reference.kinetic[begin:interval:end], label="Reference", color=:black)

    data = select_data(results; dt=100.0, tmax=75.0, trajectories=1e3, bath_type=:ShenviGaussLegendre, initial_state=:Adiabatic)
    plot_data!(data, "Gauss-Legendre", ax1, ax2, COLORS[1], "short_gauss")
    data = select_data(results; dt=100.0, tmax=75.0, trajectories=1e3, bath_type=:TrapezoidalRule, initial_state=:Adiabatic)
    plot_data!(data, "Trapezoid", ax1, ax2, COLORS[3], "short_trapezoid")

    data = select_data(results; dt=100.0, tmax=10000.0, trajectories=1e3, bath_type=:ShenviGaussLegendre, initial_state=:Adiabatic)
    plot_data!(data, "NQCDynamics.jl", ax3, ax4, COLORS[1], "long_gauss")
    data = select_data(results; dt=100.0, tmax=10000.0, trajectories=1e3, bath_type=:TrapezoidalRule, initial_state=:Adiabatic)
    plot_data!(data, "NQCDynamics.jl", ax3, ax4, COLORS[3], "long_trapezoid")

    hlines!(ax3, 0.1, linestyle=dash, color=(:black, 0.6))
    hlines!(ax4, 0.5, linestyle=dash, color=(:black, 0.6))

    Legend(fig[1,1], ax2; tellwidth=false, tellheight=true, orientation=:horizontal, valign=:center, halign=:center, margin=(0, 0, 0, 0), nbanks=1)

    return fig
end

save_figure(plotsdir("fig3_double_well.png"), plot_data, ncolors=4)
