
using DrWatson
using DataFrames
using FileIO
using Statistics
using StatsBase: skewness, kurtosis
include("StrainSurveys.jl")
include("SpatialCorrelations.jl")

abstract type AbstractStrainStatistic end
name(stat::Type{<: AbstractStrainStatistic}, variable) = Symbol( String(stat |> Symbol) * String(variable) )

# Cumulant expressions
struct Mean <: AbstractStrainStatistic end
struct Variance <: AbstractStrainStatistic end
struct Skewness <: AbstractStrainStatistic end
struct Kurtosis <: AbstractStrainStatistic end

Mean(x) = mean(x)
Variance(x) = var(x)
Skewness(x) = skewness(x)
Kurtosis(x) = kurtosis(x)

function cycleStatistics(x, varname)
    colnames = Symbol[]
    colvals  = Vector{Real}[]
    for stat ∈ subtypes(AbstractStrainStatistic)
        push!(colnames, name(stat, varname))
        push!(colvals, [stat(x)])
    end
    return colnames, colvals
end

function load_dataset_names(dir)
    names = String[]
    for fl ∈ readdir(dir; join = true)
        fl |> isfile ? push!(names, fl) : nothing
    end
    return names
end

function read_results_parameters(filename, param_type::Type{<: SimulationParameters} = DistributionParameters)
    results = FileIO.load(filename)
    params = parse_savename(param_type, filename)
    colnames = [:L, :ndis, :rtol, :cratio, :nsamples]
    colvals  = [ [params.Lx], [params.ndislocations], [params.rtol], [params.cratio], [params.nsamples] ]
    return results, params, colnames, colvals
end

function extract_single_DataFrame(filename, param_type::Type{<: SimulationParameters} = DistributionParameters )
    results, params, colnames, colvals = read_results_parameters(filename, param_type)

    # Get ndata
    nvals = length.(results["dislocations"])
    n_names, n_vals = cycleStatistics(nvals, :NumDislocations)
    append!(colnames, n_names)
    append!(colvals, n_vals)

    # Get fields data
    fields_df = send_to_DataFrame(results["strains"], params)
    for field ∈ names(fields_df)
        field_names, field_stats = cycleStatistics(fields_df[!, field], field)
        append!(colnames, field_names)
        append!(colvals, field_stats)
    end

    scalar_df =  DataFrame(colvals, colnames)
    corr_df = extract_single_Correlations(filename, param_type)
    return hcat(scalar_df, corr_df)
end

function find_DataFrames(all_files, args...)
    single_dfs = Vector{DataFrame}(undef, length(all_files))

    @inbounds Threads.@threads for idx ∈ eachindex(all_files)
        single_dfs[idx] = extract_single_DataFrame(all_files[idx], args...)
    end

    output_df = single_dfs[begin]
    @inbounds for df_idx ∈ eachindex(single_dfs)
        df_idx == oneunit(df_idx) ? continue : nothing
        output_df = vcat(output_df, single_dfs[df_idx])
    end
    return output_df
end


# begin
#     plt = plot( dfs[!, :n], dfs[!, :ΔB1g]; label = "ΔB1g", markershape = :cross)
#     plot!(plt, dfs[!, :n], dfs[!, :ΔB2g]; label = "ΔB2g", markershape = :xcross)
#     plot!(plt, dfs[!, :n], dfs[!, :ΔΔ]; label = "ΔΔ", markershape = :utriangle)
#     plot!(plt; xlabel = "Expected Number of Dislocations", ylabel = "Field Standard Deviation", title = "\$L = $(dfs[1, :L]) \$" )
# end