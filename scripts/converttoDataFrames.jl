
using DrWatson
using CSV
include("src/datatoDataFrames.jl")

function convertToDataFrames(dir, filesavename)
    filenames = load_dataset_names(dir)
    @info "Loading data as a single DataFrame"
    df = find_DataFrames(filenames)
    write_file = joinpath(dir, "DataFrames", filesavename)
    @info "Writing single DataFrame as $(write_file)"
    CSV.write(write_file, df)
    return df
end

@time convertToDataFrames(datadir("Agate"), "unit_cratio_sweep.csv")
