
using DrWatson
@quickactivate "ErpenbeckThossModel"
using MKL
using NQCDynamics
using Unitful, UnitfulAtomic
using LinearAlgebra: BLAS
BLAS.set_num_threads(1)

function termination_condition(u, t, integrator)::Bool
    return (t > austrip(10u"fs")) && (DynamicsUtils.get_positions(u)[1] > austrip(5u"Å"))
end
terminate = DynamicsUtils.TerminatingCallback(termination_condition)

function run_simulation(params)
    @unpack dt, nstates, gamma, trajectories, width, mass, discretisation, voltage, method, incident_energy, tmax, temperature = params

    atoms = Atoms(mass*u"u")
    Γ = austrip(gamma * u"eV")
    thossmodel = ErpenbeckThoss(;Γ, m=atoms.masses[1])

    bandmin = - austrip(width * u"eV") / 2 + austrip(voltage/2 * u"eV")
    bandmax = + austrip(width * u"eV") / 2 + austrip(voltage/2 * u"eV")
    @info "Band parameters" bandmin bandmax
    bath = eval(discretisation)(nstates, bandmin, bandmax)

    model = AndersonHolstein(thossmodel, bath; fermi_level=voltage/2*u"eV")
    temperature = austrip(temperature * u"K")

    m = atoms.masses[1]
    position = austrip(5u"Å")
    ke = austrip(incident_energy * u"eV")
    velocity = - sqrt(2ke / m)
    tmax = austrip(tmax * u"fs")

    sim = Simulation{eval(method)}(atoms, model; temperature,
        friction_method=ClassicalMethods.WideBandExact((nstates+1)/(bandmax-bandmin), 1/austrip(100u"K"))
    )
    dist = DynamicalDistribution(velocity, position, (1,1))

    return run_dynamics(sim, (0.0, tmax), dist;
        output=(OutputKineticEnergy, OutputPosition),
        dt=dt*u"fs",
        callback=terminate,
        trajectories,
        reduction=FileReduction(datadir("fig9_scattering_loss", savename(params, "h5"))),
    )
end

all_params = Dict{String,Any}(
    "trajectories"=>[1],
    "nstates"=>[300],
    "dt"=>[0.05],
    "width"=>[256],
    "mass"=>[10.54],
    "gamma"=>[0.02, 0.05, 0.1, 0.25, 0.5, 1.0],
    "temperature"=>[0.0],
    "tmax"=>[300],
    "voltage"=>[0.0],
    "discretisation"=>[:TrapezoidalRule],
    "method"=>[:DiabaticMDEF],
    "incident_energy"=>collect(0.25:0.25:5)
)
params = dict_list(all_params)[1] # Change index to change parameters

@time run_simulation(params)
