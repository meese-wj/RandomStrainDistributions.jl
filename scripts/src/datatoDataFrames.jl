
using DrWatson
using DataFrames
using FileIO
using Statistics
using StatsBase: skewness, kurtosis
using Distributions
include("StrainSurveys.jl")
include("ChiSqGoF.jl")
include("SpatialCorrelations.jl")

abstract type AbstractStrainStatistic end
abstract type AbstractCyclableStatistic <: AbstractStrainStatistic end
abstract type AbstractSplittingStatistic <: AbstractStrainStatistic end
name(stat::Type{<: AbstractStrainStatistic}, variable) = Symbol( String(stat |> Symbol) * String(variable) )

# Cumulant expressions
struct Mean <: AbstractCyclableStatistic end
struct Variance <: AbstractCyclableStatistic end
struct Skewness <: AbstractCyclableStatistic end
struct Kurtosis <: AbstractCyclableStatistic end

# Δ distribution fit (these aren't cyclable because they only apply for Δ)
struct HistogramFit <: AbstractSplittingStatistic end
struct GammaDistFit <: AbstractSplittingStatistic end
struct χ2GammaDistFit <: AbstractSplittingStatistic end

const HISTBINS::Int = 16

Mean(x) = mean(x)
Variance(x) = var(x)
Skewness(x) = skewness(x)
Kurtosis(x) = kurtosis(x)
HistogramFit(x, nbins = HISTBINS) = histogram_fit(x, nbins)
GammaDistFit(x) = fit(Gamma, x)
χ2GammaDistFit(x, nbins = HISTBINS) = χsqGoFTest(x, nbins, GammaDistFit(x))

function computeStatistics(x, varname, stattype::Type{<: AbstractStrainStatistic})
    (varname !== :Δ && stattype === AbstractSplittingStatistic) && throw(MethodError(computeStatistics, (x, varname, stattype)))

    colnames = Symbol[]
    colvals  = []
    for stat ∈ subtypes(stattype)
        push!(colnames, name(stat, varname))
        push!(colvals, [stat(x)])
    end
    return colnames, colvals
end

cycleStatistics(x, varname) = computeStatistics(x, varname, AbstractCyclableStatistic)
splittingStatistics(x) = computeStatistics(x, :Δ, AbstractSplittingStatistic)

function load_dataset_names(dir)
    names = String[]
    for fl ∈ readdir(dir; join = true)
        fl |> isfile ? push!(names, fl) : nothing
    end
    @info "$(length(names)) files found in $dir"
    return names
end

function read_results_parameters(filename, param_type::Type{<: SimulationParameters} = DistributionParameters)
    results = FileIO.load(filename)
    params = parse_savename(param_type, filename)
    colnames = [:L, :ndis, :con, :rtol, :cratio, :nsamples]
    colvals  = Vector{Any}[ [params.Lx], [params.ndislocations], [concentration(params)], [params.rtol], [params.cratio], [params.nsamples] ]
    return results, params, colnames, colvals
end

function extract_single_DataFrame(filename, param_type::Type{<: SimulationParameters} = DistributionParameters; correlations::Bool = true )
    results, params, colnames, colvals = read_results_parameters(filename, param_type)

    # Get ndata
    nvals = length.(results["dislocations"])
    n_names, n_vals = cycleStatistics(nvals, :NumDislocations)
    append!(colnames, n_names)
    append!(colvals, n_vals)

    # Get cycleable fields data
    fields_df = send_to_DataFrame(results["strains"], params)
    for field ∈ names(fields_df)
        field_names, field_stats = cycleStatistics(fields_df[!, field], field)
        append!(colnames, field_names)
        append!(colvals, field_stats)
    end

    # Get splitting data
    Δnames, Δstats = splittingStatistics(fields_df[!, :Δ])
    append!(colnames, Δnames)
    append!(colvals, Δstats)

    # Collect DataFrame and include correlations (if correlations == true)
    scalar_df =  DataFrame(colvals, colnames)
    if correlations
        corr_df = extract_single_Correlations(filename, param_type)
        scalar_df = hcat(scalar_df, corr_df)
    end

    return scalar_df
end

function find_DataFrames(all_files, args...; kwargs...)
    single_dfs = Vector{DataFrame}(undef, length(all_files))

    @info "Calculating single DataFrames..."
    @inbounds Threads.@threads for idx ∈ eachindex(all_files)
        single_dfs[idx] = extract_single_DataFrame(all_files[idx], args...; kwargs...)
    end

    @info "Concatenating DataFrames..."
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