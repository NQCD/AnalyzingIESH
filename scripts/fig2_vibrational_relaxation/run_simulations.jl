
using DrWatson
@quickactivate "MiaoSubotnikIESH"
using MKL
using Libdl: DL_LOAD_PATH
using NQCDynamics
using FortranNOAu111
using Unitful, UnitfulAtomic
using LinearAlgebra: diagm, norm
using LinearAlgebra
BLAS.set_num_threads(1)

using DelimitedFiles

dt = 0.1u"fs"
filename = srcdir("surface_Au111.dat")
freeze_layers = 1
Evib = austrip(0.58u"eV")

function read_cell(filename)
    open(filename) do io
        _ = parse(Int, readline(io)) # discard number of gold atoms
        cellsize = [parse(Float64, readline(io)) for _=1:3] .* u"Å"
        order = [2, 3, 1]
        cell = PeriodicCell(diagm(austrip.(cellsize[order])))
        return cell
    end
end

function load_positions(filename, site)
    # All quantities in Angstrom, convert to atomic units on return
    xy_coords = Dict(:fcc => (1.47609, 0.852218), :hcp => (2.95217, 1.70444))
    heights = Dict(:fcc => 1.808320359272221, :hcp => 1.7358814769105728)

    no_positions = zeros(3,2)
    r₀ = 1.15077

    x, y = xy_coords[site]
    height = heights[site]
    no_positions[1,:] .= x
    no_positions[2,:] .= y
    no_positions[3,1] = height
    no_positions[3,2] = height + r₀

    au_positions = permutedims(readdlm(filename; skipstart=4))
    return austrip.(hcat(no_positions, au_positions) .* u"Å")
end

function output(sol, i)

    sim = sol.prob.p

    function calculate_vibrational_energy(sim, u)
        x = DynamicsUtils.get_positions(u)
        v = DynamicsUtils.get_velocities(u)
        r = norm(x[:,1] .- x[:,2])

        ṙ = 0
        for j=1:3
            ṙ += (v[j,1]-v[j,2])*(x[j,1]-x[j,2]) / r
        end

        μ = InitialConditions.QuantisedDiatomic.reduced_mass(sim.atoms.masses[1:2])
        return μ/2*ṙ^2
    end

    vibrational_energy = zeros(length(sol.u))
    for (i, u) in enumerate(sol.u)
        vibrational_energy[i] = calculate_vibrational_energy(sim, u)
    end

    return vibrational_energy
end

function run_simulation(params)

    @unpack steps, trajectories, bath_states, site, method = params

    push!(DL_LOAD_PATH, projectdir("lib"))

    atoms = Atoms(vcat([:N, :O], fill(:Au, 528)))
    r = load_positions(filename, site)

    v = zero(r)
    vibrational_velocities = zeros(2)
    vibrational_velocities[1] = +sqrt(Evib / sum(atoms.masses[1]))
    vibrational_velocities[2] = -sqrt(Evib / sum(atoms.masses[2]))
    v[3,1:2] .= vibrational_velocities

    bandwidth = austrip(7u"eV")
    bath = ReferenceGaussLegendre(bath_states, -bandwidth/2, bandwidth/2)
    @info "Ignore the warning, it is necessary to additionally scale the coupling parameter in the Fortran model."
    model = AndersonHolstein(FortranNOAu111Model(r; freeze_layers), bath)

    sim = Simulation{eval(method)}(atoms, model)
    distribution = DynamicalDistribution(v, r, (3, length(atoms)))

    return run_dynamics(sim, (0.0, steps*dt), distribution;
        trajectories, dt, output,
        reduction=MeanReduction()
    )
end

all_params = Dict{String,Any}(
    "batch"=>1,
    "trajectories"=>[1000],
    "steps"=>[5000],
    "bath_states"=>[40, 60, 80, 90],
    "method"=>[:AdiabaticIESH],
    "site"=>[:hcp],
)
params = dict_list(all_params)[1]

@time result = run_simulation(params)
params["vibrational"] = result[:output]

@tagsave(datadir("fig2_vibrational_relaxation"), savename(params, "jld2"))
