using DrWatson
@quickactivate "MiaoSubotnikIESH"
using DataFrames
using CSV
using Unitful, UnitfulAtomic
using StatsBase
using LinearAlgebra: normalize
using Peaks

include(srcdir("data.jl"))

results = collect_results!(datadir("fig2_vibrational_relaxation"))
Ms = [20, 40, 60, 80, 90]

function process_relaxation(batches; kwargs...)
    if batches > 1
        avg = zeros(5002)
        for b in 1:batches
            data = select_data(results; batch=b, kwargs...)
            avg += data.vibrational
        end
        relaxation = ustrip(auconvert.(u"eV", avg)) ./ batches
    else
        data = select_data(results; kwargs...)
        relaxation = ustrip(auconvert.(u"eV", data.vibrational))
    end

    t = range(0, 500, length=5002)
    pks, vals = findmaxima(relaxation, 10)
    pushfirst!(pks, 1)
    pushfirst!(vals, relaxation[1])
    if !ismissing(kwargs[:method])
        writedlm(projectdir("figure_data", "fig2", "nqcd_$(kwargs[:bath_states])_adiabatic.txt"), [t[pks] vals])
    else
        writedlm(projectdir("figure_data", "fig2", "nqcd_$(kwargs[:bath_states]).txt"), [t[pks] vals])
    end
end

process_relaxation(1; bath_states=40, site=:hcp, method=:Classical)
process_relaxation(20; bath_states=20, site=:hcp, method=missing)
process_relaxation(20; bath_states=40, site=:hcp, method=missing)
process_relaxation(20; bath_states=60, site=:hcp, method=missing)
process_relaxation(20; bath_states=80, site=:hcp, method=missing)
process_relaxation(20; bath_states=90, site=:hcp, method=missing)
