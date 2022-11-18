
using DrWatson
using HDF5
using DataFrames

get_check(::Missing) = x->ismissing.(x)
get_check(v::Any) = x->x .== v

function select_data(results; kwargs...)
    df = copy(results)
    for (k, v) in kwargs
        subset!(df, k => get_check(v); skipmissing=true)
    end
    if nrow(df) == 1
        return df[1,:]
    else
        display(kwargs)
        throw(error("Data incorrect? Found $(nrow(df)) rows instead of 1."))
    end
end

function select_data_multiple(results; kwargs...)
    df = copy(results)
    for (k, v) in kwargs
        subset!(df, k => get_check(v); skipmissing=true)
    end
    return df
end

function read_data(directory, all_params; accesses)
    parameter_files = dict_list(all_params)
    output = Dict{String,Any}()
    for p in parameter_files
        filename = savename(p, "h5")
        keyname = savename(p; accesses)
        try
        h5open(datadir(directory, filename), "r") do fid
            haskey(output, keyname) && throw(error("Duplicate data found. "))
            output[keyname] = []
            for traj in keys(fid)
                group = fid[traj]
                push!(output[keyname], Dict(quantity => Array(group[quantity]) for quantity in keys(group)))
            end
        end
        catch
            @show filename
        end
    end
    return output
end

select_data_entry(results; kwargs...) = results[savename(kwargs)]

