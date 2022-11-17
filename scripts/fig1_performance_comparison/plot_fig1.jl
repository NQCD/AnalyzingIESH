using DrWatson
@quickactivate "AnalyzingIESH"
using JamesPlots
using CairoMakie
using DataFrames
using CSV
using Statistics: mean, median
using StatsBase: sem
using ColorTypes: ColorTypes
using DelimitedFiles: writedlm

results = collect_results!(datadir("fig1_performance_comparison"))
gdf = groupby(results, [:gamma, :blas, :threads])
averaged = combine(gdf, :runtime => mean, :runtime => sem, :runtime => minimum, :runtime => median, nrow)

reference = Dict(
    1e-4 => 600,
    4e-4 => 27,
    1.6e-3 => 15,
    6.4e-3 => 10,
)
palette = JamesPlots.choose_color(3)

dash = [0.0, 3.2, 4.0]
include(srcdir("data.jl"))

function plot_data()
    fig = Figure(figure_padding=(1, 1, 1, 1),
        resolution=(JamesPlots.COLUMN_WIDTH, 1.0*JamesPlots.RESOLUTION[2]),
        font=projectdir("fonts", "MinionPro-Capt.otf"),
    )

    limits = (0.3, nothing, 0.0, nothing)
    yticklabelpad = 1.0
    axes = [Axis(fig[1,i]; limits, yticklabelpad) for i=1:4]
    hidespines!.(axes, :r, :t)
    axes[1].ylabel = "Elapsed time per trajectory /min"
    hidexdecorations!.(axes; label=false)

    axes[4].yticks=0:150:600

    gammas = reverse([1e-4, 4e-4, 1.6e-3, 6.4e-3])

    for i=1:4
        gamma = gammas[i]

        heights = Float64[]
        for blas in [:MKL, :OpenBLAS]
            for threads in [1, 4, 8]
                y = 0.0
                try
                    y = select_data(averaged; gamma, blas, threads).runtime_minimum / 60
                catch
                end
                push!(heights, y)
            end
        end
        push!(heights, reference[gamma])
        @show heights
        writedlm(projectdir("figure_data", "fig1", "$i.txt"), heights)

        xs = repeat(1:2; inner=3)
        dodge = repeat(1:3; outer=2)
        push!(xs, 3)
        push!(dodge, 1)
        
        color = repeat(palette[1:2]; inner=3)
        push!(color, palette[3])

        barplot!(axes[i], xs, heights; dodge, color, strokewidth=1, gap=0.00, dodge_gap=0.0)
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
        tellwidth=true, tellheight=true, orientation=:horizontal, margin=(0, 0, 0, 4),
        nbanks=1, patchsize=(6, 6)
    )
    colgap!(fig.layout, 1.5)

    offset = 0.335
    for i=1:4
        text!(axes[i], 1-offset, 0; text="1", textsize=8, align=(:center, :bottom))
        text!(axes[i], 1, 0; text="4", textsize=8, align=(:center, :bottom))
        text!(axes[i], 1+offset, 0; text="8", textsize=8, align=(:center, :bottom))

        text!(axes[i], 2-offset, 0; text="1", textsize=8, align=(:center, :bottom))
        text!(axes[i], 2, 0; text="4", textsize=8, align=(:center, :bottom))
        text!(axes[i], 2+offset, 0; text="8", textsize=8, align=(:center, :bottom))
    end

    return fig
end

save_figure(plotsdir("fig1_performance_comparison.png"), plot_data, ncolors=2)
