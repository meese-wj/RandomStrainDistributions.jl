
using DrWatson
using DataFrames
using FileIO
using Statistics
using StatsBase: skewness, kurtosis
using Distributions
include("recursive_subtypes.jl") # for fancy subtyping
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
abstract type AbstractHistogramFit <: AbstractSplittingStatistic end
abstract type AbstractDenseHistogramFit <: AbstractHistogramFit end
struct HistogramFit <: AbstractHistogramFit end
struct DenseHistogramFit <: AbstractDenseHistogramFit end
struct NormHistogramFit <: AbstractHistogramFit end
struct DenseNormHistogramFit <: AbstractDenseHistogramFit end

const HISTBINS::Int = 64
default_num_bins(::AbstractHistogramFit) = HISTBINS
const DENSEHISTBINS::Int = 8 * HISTBINS
default_num_bins(::AbstractDenseHistogramFit) = DENSEHISTBINS

hist_normalized(::HistogramFit) = false
hist_normalized(::DenseHistogramFit) = false
hist_normalized(::NormHistogramFit) = true
hist_normalized(::DenseNormHistogramFit) = true

hist_density(::AbstractHistogramFit) = false
hist_density(::AbstractDenseHistogramFit) = true

struct GammaDistFit <: AbstractSplittingStatistic end
struct PseudoRayleighFit <: AbstractSplittingStatistic end
struct χ2GammaDistFit <: AbstractSplittingStatistic end
struct χ2PseudoRayleighFit <: AbstractSplittingStatistic end
struct GammaDistcV <: AbstractSplittingStatistic end
struct PseudoRayleighcV <: AbstractSplittingStatistic end
struct HistogramcV <: AbstractSplittingStatistic end

Mean(x) = mean(x)
Variance(x) = var(x)
Skewness(x) = skewness(x)
Kurtosis(x) = kurtosis(x)

# Histograms -- generate their own constructors based on traits
# This will in principle make it less likely for me to mess stuff up...
for hist_type ∈ recursive_subtypes(AbstractHistogramFit)
    if isconcretetype(hist_type)
        hist_symb = Symbol(hist_type)
        quote 
            $(hist_symb)(x) = histogram_fit(x, default_num_bins( $(hist_symb)() ); 
                                            normalize = hist_normalized( $(hist_symb)() ), 
                                            density = hist_density( $(hist_symb)() ) )
        end |> eval
    end
end

# DenseNormHistogramFit(x) = histogram_fit(x, default_num_bins(DenseNormHistogramFit()); normalize = true, density = true)
GammaDistFit(x) = fit(Gamma, x)
PseudoRayleighFit(x) = fit_pseudo_rayleigh(NormHistogramFit(x))
χ2GammaDistFit(x, nbins = HISTBINS) = χsqGoFTest(x, nbins, GammaDistFit(x))
χ2PseudoRayleighFit(x, nbins = HISTBINS) = χsqGoFTest(x, nbins, pseudo_rayleigh_distribution(PseudoRayleighFit(x).param))
GammaDistcV(x) = fitcV(x, GammaDistFit(x))
PseudoRayleighcV(x) = fitcV(x, pseudo_rayleigh_distribution(PseudoRayleighFit(x).param))
HistogramcV(x, nbins = HISTBINS) = histcV(x, nbins)

# SCROLL DOWN FOR similarity_coeff ETC.

function computeStatistics(x, varname, stattype::Type{<: AbstractStrainStatistic})
    (varname !== :Δ && stattype === AbstractSplittingStatistic) && throw(MethodError(computeStatistics, (x, varname, stattype)))

    colnames = Symbol[]
    colvals  = []
    for stat ∈ recursive_subtypes(stattype)
        if isconcretetype(stat)
            push!(colnames, name(stat, varname))
            push!(colvals, [stat(x)])
        end
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
    # cV comparisons
    gammacV = scalar_df[!, :GammaDistcVΔ][1]
    rayleighcV = scalar_df[!, :PseudoRayleighcVΔ][1]
    histcV = scalar_df[!, :HistogramcVΔ][1]
    cVdf = DataFrame( [ [rmse(gammacV..., histcV[2])], [similarity_coeff(gammacV..., histcV[2])],
                        [rmse(rayleighcV..., histcV[2])], [similarity_coeff(rayleighcV..., histcV[2])] ], 
                      [:GammacVRMSEΔ, :GammacVSimCoeffΔ,
                       :PseudoRayleighcVRMSEΔ, :PseudoRayleighcVSimCoeffΔ] )
    scalar_df = hcat(scalar_df, cVdf )

    if correlations
        crosscorr_df = extract_cross_Correlations(filename, param_type)
        scalar_df = hcat(scalar_df, crosscorr_df)
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
