using DrWatson
@quickactivate "AnalyzingIESH"
using NQCDynamics
using Unitful, UnitfulAtomic
using Distributions
using LinearAlgebra
using InteractiveUtils

all_params = Dict{String,Any}(
    "repeat"=>collect(1:5),
    "trajectories"=>[1],
    "dt"=>[10.0],
    "gamma"=>[1e-4, 4e-4, 1.6e-3, 6.4e-3],
    "blas"=>[:OpenBLAS],
    "bath_type"=>[:TrapezoidalRule],
    "threads"=>[1, 4, 8],
    "estimate_probability"=>[true, false]
)
params = dict_list(all_params)[end] # Change index to change parameters

@unpack threads = params 
BLAS.set_num_threads(threads)

display(versioninfo())
println()

display(BLAS.get_config())
println()

@info "Threading info" BLAS.get_num_threads() Threads.nthreads()
println()

function run_simulation(params)
    @unpack dt, trajectories, bath_type, gamma, estimate_probability = params

    kT = 9.5e-4

    if gamma == 1e-4
        W = 100gamma
        tmax = 5e5 * dt
        M = 200
    elseif gamma == 4e-4
        W = 40gamma
        tmax = 3e5 * dt
        M = 80
    elseif gamma == 1.6e-3
        W = 5gamma
        tmax = 1.8e5 * dt
        M = 40
    elseif gamma == 6.4e-3
        W = 5gamma
        tmax = 1.2e5 * dt
        M = 40
    end

    model = MiaoSubotnik(;Γ=gamma)
    bath = eval(bath_type)(M, -W, W)
    metal_model = AndersonHolstein(model, bath)
    atoms = Atoms(2000)
    sim = Simulation{AdiabaticIESH}(atoms, metal_model; estimate_probability)

    velocity = VelocityBoltzmann(5kT, atoms.masses, (1,1))
    mω² = atoms.masses[1] * model.ω^2
    σ = sqrt(5kT / mω²)
    position = Normal(0.0, σ)

    distribution = DynamicalDistribution(velocity, position, (1,1))

    tspan = (0.0, tmax)
    saveat = range(tspan[1], tspan[2], length=10)

    @info "Compile run"
    @time run_dynamics(sim, (0.0, dt), distribution; output=OutputPosition, trajectories, dt, saveat)
    @info "Compiled"

    return @elapsed run_dynamics(sim, tspan, distribution;
        output=OutputPosition,
        trajectories, dt, saveat,
    )

end

@time output = run_simulation(params)
params["runtime"] = output

@tagsave(datadir("fig1_performance_comparison", savename(params, "jld2")), params)
