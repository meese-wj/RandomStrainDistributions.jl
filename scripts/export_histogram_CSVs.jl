using Pkg
Pkg.activate(@__DIR__)

# *********************
# These are for JLD2 reconstruction
# to get rid of warnings. 
using Distributions
using HypothesisTests
using LsqFit
# *********************

# *********************
# Useful packages are below
using DataFrames
using CSV
using DrWatson
using JLD2
using Pipe
using StatsBase
# *********************

simdatadir(args...) = datadir("Agate", "DataFrames", args...)
outputdatadir(args...) = simdatadir("Histograms", args...)

data = @pipe load_object( simdatadir("unit_cratio_sweep.JLD2") ) |> sort(_, :L; rev = true)

Lsize = 256
histogram_names = [:DenseNormHistogramFitB1g, :DenseNormHistogramFitB2g, :DenseNormHistogramFitΔ]
collected_names = [:L, :con, :ndis, :nsamples, histogram_names...]

collected_data = @pipe filter(:L => ==(Lsize), data) |> _[!, collected_names] |> sort(_, :con)

function bin_centers(hist)
    edgez = hist.edges[begin]
    return edgez[1:end-1] + 0.5 * diff(edgez)
end

function extract_histogram_columns(hist)
    bins = bin_centers(hist)
    vals = hist.weights
    return bins, vals
end

function extract_new_dataframes(culled_data)
    dfs = []
    for (condx, conval) ∈ enumerate(culled_data[!, :con])
        # B1g
        bins, vals = extract_histogram_columns(culled_data[condx, :DenseNormHistogramFitB1g])
        data_df = DataFrame([fill(conval, length(bins)), bins, vals], [:c, :B1g, :pB1g])
        # B2g
        bins, vals = extract_histogram_columns(culled_data[condx, :DenseNormHistogramFitB2g])
        data_df = hcat(data_df, DataFrame([bins, vals], [:B2g, :pB2g]))
        # Delta
        bins, vals = extract_histogram_columns(culled_data[condx, :DenseNormHistogramFitΔ])
        data_df = hcat(data_df, DataFrame([bins, vals], [:Delta, :pDelta]))
        push!(dfs, data_df)
    end
    return dfs
end

all_hists = extract_new_dataframes(collected_data)
for df ∈ all_hists
    conval = round(df[1, :c]; sigdigits = 3)
    filename = outputdatadir("con=$conval.csv")
    CSV.write(filename, df)
end
