
using DrWatson
@quickactivate "MiaoSubotnikIESH"
using MKL
using NQCDynamics
using Unitful, UnitfulAtomic
using Distributions
using LinearAlgebra

function impurity_population(sol, i)
    sim = sol.prob.p
    out = zeros(length(sol), 2)
    for i in eachindex(sol.u)
        out[i,1] = 1 - Estimators.diabatic_population(sim, sol.u[i])[1]
        out[i,2] = Estimators.kinetic_energy(sim, sol.u[i])
    end
    return out
end

function run_simulation(params)
    @unpack dt, trajectories, tmax, bath_type, initial_state = params

    kT = 9.5e-4
    M = 40 # number of bath states
    Γ = 6.4e-3
    W = 10Γ / 2 # bandwidth  parameter
    model = MiaoSubotnik(;Γ)
    bath = eval(bath_type)(M, -W, W)
    metal_model = AndersonHolstein(model, bath)
    atoms = Atoms(2000)
    sim = Simulation{AdiabaticIESH}(atoms, metal_model)

    velocity = VelocityBoltzmann(5kT, atoms.masses, (1,1))
    mω² = atoms.masses[1] * model.ω^2
    σ = sqrt(5kT / mω²)
    position = Normal(0.0, σ)

    if initial_state == :Diabatic
        electronic_dist = FermiDiracState(0.0, 0.0, statetype=Diabatic(), available_states=2:M+1)
    elseif initial_state == :Adiabatic
        electronic_dist = FermiDiracState(0.0, 0.0, statetype=Adiabatic())
    end

    distribution = DynamicalDistribution(velocity, position, (1,1)) * electronic_dist

    tmax = tmax / model.ω
    tspan = (0.0, tmax)
    saveat = range(tspan[1], tspan[2], length=200)

    return run_dynamics(sim, tspan, distribution;
        output=impurity_population,
        trajectories, dt, saveat,
        reduction=MeanReduction()
    )

end

all_params = Dict{String,Any}(
    "trajectories"=>[1000],
    "dt"=>[500.0],
    "tmax"=>[75.0, 10000.0],
    "bath_type"=>[:ShenviGaussLegendre, :TrapezoidalRule],
    "initial_state"=>[:Adiabatic, :Diabatic]
)
params = dict_list(all_params)[1] # Change index to change parameters

@time output = run_simulation(params)
params["output"] = output[:impurity_population]

@tagsave(datadir("fig3_double_well", savename(params, "jld2")), params)
