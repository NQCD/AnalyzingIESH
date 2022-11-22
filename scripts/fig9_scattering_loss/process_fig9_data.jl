using DrWatson
@quickactivate "ErpenbeckThossModel"
using DataFrames
using CSV
using Unitful, UnitfulAtomic
using Statistics: mean, std
using HDF5
using DelimitedFiles

include(srcdir("data.jl"))


function calculate_final_kinetic_energy(dat)
    samples = Float64[]
    for traj in eachindex(dat)
        if dat[traj]["OutputPosition"][end][1] > austrip(4u"Å")
            Ei = dat[traj]["OutputKineticEnergy"][1]
            Ef = dat[traj]["OutputKineticEnergy"][end]
            push!(samples, (Ei - Ef)/Ei * 100)
        else
            Ei = dat[traj]["OutputKineticEnergy"][1]
            Ef = austrip(1.8314118436225917u"eV")
            loss = (Ei - Ef)/Ei * 100
            loss = clamp(loss, 0, Inf)
        end
    end
    m = mean(samples)
    return (m, std(samples; mean=m) / sqrt(length(samples)))
end

function process_data(label, results, all_params)

    energies = all_params["incident_energy"]
    gammas = all_params["gamma"]
    energy_mean = zeros(length(energies), length(gammas))
    energy_std = zeros(length(energies), length(gammas))
    for (i, incident_energy) in enumerate(energies)
        for (j, gamma) in enumerate(gammas)
            try
            dat = select_data_entry(results; gamma, incident_energy)
            energy_mean[i,j], energy_std[i,j] = calculate_final_kinetic_energy(dat)
            catch; end;
            @show energy_mean[i,j], energy_std[i,j]
        end
    end

    for i=1:length(gammas)
        slice = (!isnan).(energy_mean[:,i])
        energies_slice = energies[slice]
        mean_slice = energy_mean[slice,i]
        if label == "IESH"
            writedlm(projectdir("figure_data", "fig9", "$(label)_gamma=$(gammas[i]).txt"), [energies_slice mean_slice energy_std[slice,i]])
        else
            writedlm(projectdir("figure_data", "fig9", "$(label)_gamma=$(gammas[i]).txt"), [energies_slice mean_slice])
        end
    end
end

all_params = Dict{String,Any}(
    "trajectories"=>[200],
    "nstates"=>[200],
    "dt"=>[0.01],
    "width"=>[64],
    "mass"=>[10.54],
    "gamma"=>[0.02, 0.1, 0.25, 0.5, 1.0],
    "temperature"=>[0.0],
    "tmax"=>[300],
    "voltage"=>[0.0],
    "discretisation"=>[:ShenviGaussLegendre],
    "method"=>[:AdiabaticIESH],
    "incident_energy"=>collect(0.25:0.25:5)
)

results = read_data("fig9_scattering_loss", all_params; accesses=["gamma", "incident_energy"])
@time process_data("IESH", results, all_params)

all_params = Dict{String,Any}(
    "trajectories"=>[1],
    "nstates"=>[300],
    "dt"=>[0.05],
    "width"=>[256],
    "mass"=>[10.54],
    "gamma"=>[0.02, 0.1, 0.25, 0.5, 1.0],
    "temperature"=>[0.0],
    "tmax"=>[300],
    "voltage"=>[0.0],
    "discretisation"=>[:TrapezoidalRule],
    "method"=>[:DiabaticMDEF],
    "incident_energy"=>collect(0.25:0.25:5)
)

results = read_data("fig9_scattering_loss", all_params; accesses=["gamma", "incident_energy"])
@time process_data("MDEF", results, all_params)
