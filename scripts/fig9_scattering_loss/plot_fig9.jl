using DrWatson
@quickactivate "ErpenbeckThossModel"
using CairoMakie
using DataFrames
using CSV
using JamesPlots
using ColorSchemes, Colors
using Unitful, UnitfulAtomic
using HDF5
using DelimitedFiles

include(srcdir("data.jl"))
include(srcdir("krr.jl"))
dash = [0.0, 3.2, 4.0]

colormap = ColorScheme(parse.(Colorant, ["#B7E6A5", "#7CCBA2", "#46AEA0", "#089099", "#00718B", "#045275", "#003147"]))
colors = [ColorSchemes.get(colormap, i/5) for i=0:4]

function plot_lines!(ax, label, lambda)
    gammas = [0.02, 0.1, 0.25, 0.5, 1.0]
    for i=1:length(gammas)
        data = readdlm(projectdir("figure_data", "fig9", "$(label)_gamma=$(gammas[i]).txt"))
        x = data[:,1]
        y = data[:,2]

        if size(data,2) == 3
            err = data[:,3]
            errorbars!(ax, x, y, err, color=colors[i], whiskerwidth=6)
        end

        x_test = range(0.25, 5; length=200)

        y_pred = kernel_ridge_regression(with_lengthscale(MaternKernel(), 1.5), x, y, x_test, lambda)
        lines!(ax, x_test, y_pred; color=colors[i])

        scatter!(ax, x, y, color=colors[i])
    end
end

function plot_data()
    fig = Figure(resolution=(JamesPlots.RESOLUTION[1],1.8*JamesPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), font=projectdir("fonts", "MinionPro-Capt.otf"))
    axes = [MyAxis(fig[i,j], xlabel="Initial kinetic energy /eV", ylabel="Kinetic energy loss /%", limits=(0.25, 5.0, nothing, nothing), yticks=0:20:60) for i=1:2, j=1:1]

    linkyaxes!(axes...)

    plot_reference_line!(axes)

    plot_lines!(axes[1], "IESH", 1e-1)
    plot_lines!(axes[2], "MDEF", 1e-2)

    ticks = (1:5, string.(all_params["gamma"]))
    Colorbar(fig[0,1]; limits=(0.5, 5.5), colormap=cgrad(colors, categorical=true), ticks=ticks, nsteps=5, ticklabelpad=-14, ticksvisible=false, label="Γ /eV", flip_vertical_label=true, vertical=false, labelpadding=6, minorticksvisible=true, minortickalign=1.0, minorticks=[1.5, 2.5, 3.5, 4.5], minorticksize=16)

    hidexdecorations!.(axes[1,:]; ticks=false, minorticks=false)
    Label(fig[1,1], "IESH"; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10))
    Label(fig[2,1], "MDEF"; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10))

    fig
end

function plot_reference_line!(axes)
    energies = range(0.25, 5; length=200)
    loss = (energies .- 1.8314118436225917) ./ energies .* 100
    clamp!(loss, 0, Inf)
    lines!(axes[1], energies, loss, color=:black, linestyle=dash)
    lines!(axes[2], energies, loss, color=:black, linestyle=dash)
end

save_figure(plotsdir("fig9_scattering_loss.png"), plot_data, ncolors=6)
