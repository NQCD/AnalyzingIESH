
using DrWatson
@quickactivate "ErpenbeckThossModel"
using MKL
using NQCDynamics
using Unitful, UnitfulAtomic
using LinearAlgebra: BLAS
BLAS.set_num_threads(1)
using Distributions: Normal

function termination_condition(u, t, integrator)::Bool
    return DynamicsUtils.get_positions(u)[1] > austrip(5u"Å")
end
terminate = DynamicsUtils.TerminatingCallback(termination_condition)

function desorption(sol, i)
    tfinish = ustrip(auconvert(u"fs", sol.t[end]))
    return convert_to_probability(tfinish)
end

function convert_to_probability(t)
    tvals = range(0, 200; length=100)
    out = zero(tvals)
    out[tvals .> t] .= 1
    return out
end

function run_simulation(params)
    @unpack dt, nstates, gamma, trajectories, width, distribution, mass, tmax, discretisation, voltage, shift_band, method = params
    if method === :AdiabaticIESH
        @unpack initial_state = params
    end

    atoms = Atoms(mass*u"u")
    Γ = austrip(gamma * u"eV")
    thossmodel = ErpenbeckThoss(;Γ, m=atoms.masses[1])

    if shift_band
        bandmin = - austrip(width * u"eV") / 2 + austrip(voltage/2 * u"eV")
        bandmax = + austrip(width * u"eV") / 2 + austrip(voltage/2 * u"eV")
    else
        bandmin = - austrip(width * u"eV") / 2
        bandmax = + austrip(width * u"eV") / 2
    end
    bath = eval(discretisation)(nstates, bandmin, bandmax)

    model = AndersonHolstein(thossmodel, bath; fermi_level=voltage/2*u"eV")
    T = 300u"K"
    β = 1/austrip(T)

    ω = NQCModels.AdiabaticModels.getω₀(thossmodel.morse)
    m = atoms.masses[1]
    if distribution === :wigner
        position = PositionHarmonicWigner(ω, β, m; centre=thossmodel.morse.x₀)
        velocity = VelocityHarmonicWigner(ω, β, m)
    elseif distribution === :classical
        βmω² = m * ω^2 * β
        σ = sqrt(1 / βmω²)
        position = Normal(thossmodel.morse.x₀, σ)
        velocity = Normal(0, 1/sqrt(m*β))
    else
        throw(error("Distribution: $distribution not recognised."))
    end

    sim = Simulation{eval(method)}(atoms, model)

    if method === :AdiabaticIESH
        if initial_state == :Diabatic
            electronic_dist = FermiDiracState(voltage/2 * u"eV", T, statetype=Diabatic(), available_states=2:nstates+1)
        elseif initial_state == :Adiabatic
            electronic_dist = FermiDiracState(voltage/2 * u"eV", T, statetype=Adiabatic())
        end
        dist = DynamicalDistribution(velocity, position, (1,1)) * electronic_dist
    elseif method == :Classical
        dist = DynamicalDistribution(velocity, position, (1,1))
    end

    return run_dynamics(sim, (0.0, tmax*u"fs"), dist;
        output=desorption,
        dt=dt*u"fs",
        callback=terminate,
        trajectories,
        reduction=MeanReduction(),
        savetime=false
    )
end

all_params = Dict{String,Any}(
    "distribution"=>[:wigner],
    "trajectories"=>[10000],
    "nstates"=>[100, 200, 400],
    "dt"=>[1.0],
    "gamma"=>[1.0],
    "width"=>[4, 8, 16, 32, 54, 128, 256],
    "mass"=>[10.54],
    "tmax"=>[200],
    "voltage"=>[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    "shift_band"=>[true],
    "discretisation"=>[:ShenviGaussLegendre, :TrapezoidalRule],
    "method"=>[:AdiabaticIESH],
    "initial_state"=>@onlyif("method" == :AdiabaticIESH, [:Adiabatic]),
)
params = dict_list(all_params)[1] # Change index for other parameters

@time output = run_simulation(params)
params["probability"] = output[:desorption]

@tagsave(datadir("fig5_6_7_8_desorption_probability", savename(params, "jld2")), params)
