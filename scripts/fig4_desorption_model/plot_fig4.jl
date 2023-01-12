using DrWatson
@quickactivate "AnalyzingIESH"
using NQCDynamics
using NQCDynamics: Calculators
using CairoMakie
using Unitful, UnitfulAtomic
using LinearAlgebra
using JamesPlots
using ColorSchemes
using Colors
using Combinatorics

colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]))
colormap = JamesPlots.NICECOLORS

atoms = Atoms(10.54u"u")
x_ang = range(1, 5, length=200)
x = austrip.(x_ang .* u"Å")

function build_sim(gamma)
    Γ = austrip(gamma * u"eV")
    thossmodel = ErpenbeckThoss(;Γ)
    nstates = 100
    bandmin = -austrip(16u"eV")
    bandmax = austrip(16u"eV")
    bath = ShenviGaussLegendre(nstates, bandmin, bandmax)
    model = AndersonHolstein(thossmodel, bath)
    return Simulation(atoms, model)
end

function plot_system_hamiltonian!(ax, gamma, r)
    sim = build_sim(gamma)

    U0 = zeros(length(r))
    U1 = zeros(length(r))
    coupling = zeros(length(r))
    for (i, x) in enumerate(r)
        x = austrip(x * u"Å")
        V = NQCModels.potential(sim.calculator.model.model, hcat(x))
        U0[i] = ustrip.(auconvert.(u"eV", V[1,1]))
        U1[i] = ustrip.(auconvert.(u"eV", V[2,2]))
        coupling[i] = ustrip(auconvert(u"eV", V[2,1]))^2
    end
    lines!(ax, r, U0, color=COLORS[3], label=L"U_0(x)")
    lines!(ax, r, U1, color=COLORS[5], label=L"U_1(x)")
    lines!(ax, r, coupling, color=COLORS[1], label=L"V_k^2(x)")
end

function plot_iesh_surfaces!(ax, gamma, r)
    sim = build_sim(gamma)
    model = sim.calculator.model
    nelectrons = NQCModels.nelectrons(model)
    nstates = NQCModels.nstates(model)-1

    x = austrip.(r * u"Å")
    matrix_x = hcat.(x)
    extra = ustrip.(auconvert.(u"eV", NQCModels.state_independent_potential.(model, matrix_x)))
    eigs = eigvals.(potential.(model, matrix_x))
    config = vcat(ones(nelectrons), zeros(nelectrons+1))

    perms = multiset_permutations(config, nstates+1)
    energies = []
    sortpoints = []
    for config in Iterators.take(perms, 20000)
        e = [sum(e[config .!= 0]) for e in eigs]
        push!(energies, ustrip.(auconvert.(u"eV", e)))
        push!(sortpoints, e[end])
    end
    perm = sortperm(sortpoints)#[1:1600]
    energies = energies[perm]
    energies = energies[1:50:2000]

    color_range = []
    max_e = energies[end][end]
    min_e = energies[1][end]
    for e in energies
        c = (e[end] - min_e) / (max_e - min_e)
        push!(color_range, c)
    end
    colors = [ColorSchemes.get(colormap, 1-c) for c in color_range]

    for (i, e) in Iterators.reverse(enumerate(energies))
        y = e .+ extra
        lines!(ax, x_ang, y, color=colors[i])
    end
end

function plot_model_diabats()
    fig = Figure(figure_padding=(1, 2, 1, 4), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), resolution=(JamesPlots.RESOLUTION[1], JamesPlots.RESOLUTION[2]*1.5))

    ax1 = MyAxis(fig[1,1]; ylabel="Energy /eV", limits=(1, 5, -2, 6), xlabel="x /Å")
    ax2 = MyAxis(fig[2,1]; ylabel="Energy /eV", limits=(1, 5, -402, -394), xlabel="x /Å", yticks=-402:2:-396)
    hidexdecorations!(ax1, ticks=false, minorticks=false)

    gamma = 1
    plot_system_hamiltonian!(ax1, gamma, x_ang)
    plot_iesh_surfaces!(ax2, 1, x_ang)

    Legend(fig[1,1], ax1, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:horizontal)
    return fig
end

save_figure(plotsdir("fig4_desorption_model.png"), plot_model_diabats)
