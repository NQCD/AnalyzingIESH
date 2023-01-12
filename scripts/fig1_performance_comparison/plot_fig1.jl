using DrWatson
@quickactivate "MiaoSubotnikIESH"
using JamesPlots
using CairoMakie
using DataFrames
using CSV
using Statistics: mean, median
using StatsBase: sem
using ColorTypes: ColorTypes
using DelimitedFiles: writedlm

results = collect_results!(datadir("fig1_performance_comparison"))
gdf = groupby(results, [:estimate_probability, :gamma, :blas, :threads])
estimate = combine(gdf, :runtime => mean, :runtime => sem, :runtime => minimum, :runtime => median, :runtime=>maximum, nrow)

reference = Dict(
    1e-4 => 600,
    4e-4 => 27,
    1.6e-3 => 15,
    6.4e-3 => 10,
)
palette = JamesPlots.choose_color(3)

dash = [0.0, 3.2, 4.0]
include(srcdir("data.jl"))

function plot_chart!(ax, gamma, i)

    heights = Float64[]
    otherheights = Float64[]
    estimates = (true, false)
    for (container, estimate_probability) in zip((heights, otherheights), estimates)
        for blas in [:MKL, :OpenBLAS]
            for threads in [1, 4, 8]
                y = 0.0
                try
                    y = select_data(estimate; estimate_probability, gamma, blas, threads).runtime_minimum / 60
                catch
                end
                push!(container, y)
            end
        end
    end
    push!(heights, reference[gamma])
    push!(otherheights, reference[gamma])
    writedlm(projectdir("figure_data", "fig1", "with_estimate_$i.txt"), heights)
    writedlm(projectdir("figure_data", "fig1", "without_estimate_$i.txt"), otherheights)

    xs = repeat(1:3; inner=2)
    dodge = repeat(1:2; outer=3)

    push!(xs, 4)
    push!(dodge, 1)
    
    color = repeat(palette[1:2]; inner=3)
    push!(color, palette[3])

    errorbars!(ax, xs[begin:end-1] .+ vcat(repeat([-0.25, 0.25], outer=3)), heights[begin:end-1], zero(heights)[begin:end-1], otherheights[begin:end-1] .- heights[begin:end-1]; whiskerwidth=3, linewidth=1)
    barplot!(ax, xs, heights; dodge, color, strokewidth=1, gap=0.00, dodge_gap=0.0)

    ax.xticks = (xs[begin:end-1] .+ vcat(repeat([-0.25, 0.25], outer=3)), string.([1, 4, 8, 1, 4, 8]))
    ax.xticklabelpad = 0.0
end

function plot_data()
    fig = Figure(figure_padding=(1, 1, 1, 1),
        resolution=(JamesPlots.COLUMN_WIDTH, 1.0*JamesPlots.RESOLUTION[2]),
        fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")),
    )

    limits = (0.3, nothing, 0.0, nothing)
    yticklabelpad = 1.0
    axes = [Axis(fig[1,i]; limits, yticklabelpad) for i=1:4]
    hidespines!.(axes, :r, :t)
    axes[1].ylabel = "Elapsed time per trajectory /min"
    hidexdecorations!.(axes; label=false, ticklabels=false)

    axes[4].yticks=0:150:600

    gammas = reverse([1e-4, 4e-4, 1.6e-3, 6.4e-3])

    for (i, gamma) in enumerate(gammas)
        plot_chart!(axes[i], gamma, i)
    end

    Ms = [
        "M = 40"
        "M = 40"
        "M = 80"
        "M = 200"
    ]
    Ts = [
        "T = 1.2 × 10⁵"
        "T = 1.8 × 10⁵"
        "T = 3.0 × 10⁵"
        "T = 5.0 × 10⁵"
    ]
    for (i, (M, T)) in enumerate(zip(Ms, Ts))
        Label(fig[-1,i], M; tellwidth=false, tellheight=true, halign=:left, valign=:top, padding=(4, 0, 0, 0), justification=:left)
        Label(fig[0,i], T; tellwidth=false, tellheight=true, halign=:left, valign=:top, padding=(4, 0, 0, 0), justification=:left)
    end

    legend_elements = [PolyElement(color=palette[i], strokecolor=:black, strokewidth=1) for i=1:3]
    legend_labels = [
        "MKL"
        "OpenBLAS"
        "Ref. 11"
    ]

    Legend(fig[end+1,:], legend_elements, legend_labels,
        tellwidth=true, tellheight=true, orientation=:horizontal, margin=(0, 0, 0, 2),
        nbanks=1, patchsize=(6, 6)
    )
    colgap!(fig.layout, 1.0)

    return fig
end

save_figure(plotsdir("fig1_performance_comparison.png"), plot_data, ncolors=2)
