using Random
using RandomStrainDistributions
using DrWatson
using FileIO
using JLD2
using DataFrames
using CSV
include("DrWatsonHelpers.jl")
include("DistributionParameters.jl")

function survey_trial(rsd, params)
    dis_vectors = collect_dislocations(rsd)
    strains = generate_disorder(params.Lx, params.Ly, dis_vectors; include_Δ = true, tolerance = params.rtol, coupling_ratio = params.cratio)
    return dis_vectors, strains
end

function survey!(dislocations, strain_fields, rsd, params)
    @inbounds Threads.@threads for trial ∈ 1:params.nsamples
        dis_vectors, strains = survey_trial(rsd, params)
        append!(dislocations, dis_vectors)
        strain_fields[:, :, :, trial] .= strains
    end
    return nothing
end

function send_to_DataFrame(strain_fields, params)
    fields_per_strain = nsize(params) * params.nsamples
    values = [
        (@view strain_fields[:, :, 1, :]),
        (@view strain_fields[:, :, 2, :]),
        (@view strain_fields[:, :, 3, :])
    ]
    return DataFrame( map( x -> reshape(x, (fields_per_strain, )), values ), [:B1g, :B2g, :Δ] )
end

function main(params; seed = 42, prefix = "strain-survey", save_output = true)
    Random.seed!(seed)
    
    rsd = RandomDislocationDistribution(; concentration = concentration(params),
                                          Lx = params.Lx, Ly = params.Ly, 
                                          burgers_vectors = tetragonal_burgers_vectors )

    dislocations = Dislocation2D{typeof(params.cratio)}[]                                    
    strain_fields = zeros(typeof(params.cratio), params.Lx, params.Ly, 3, params.nsamples)

    survey!(dislocations, strain_fields, rsd, params)

    results = Dict("dislocations" => dislocations, "strains" => strain_fields)

    jld2name = savename(prefix, params, "jld2") |> datadir
    
    if save_output
        save(jld2name, results)
    end
    return nothing 
end