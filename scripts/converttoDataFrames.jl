#!/usr/bin/bash -l
#SBATCH --time=23:30:00
#SBATCH --partition=msibigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=101
#SBATCH --mem-per-cpu=15g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --job-name=ToData
#SBATCH -o %x-%A.out
#=
    pwd
    module swap julia julia/1.8.0
    module list
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK converttoDataFrames.jl
    exit
=#

using Pkg
Pkg.activate(@__DIR__)

using DrWatson
using JLD2
include("src/datatoDataFrames.jl")

function convertToDataFrames(dir, filesavename)
    filenames = load_dataset_names(dir)
    @info "Loading data as a single DataFrame"
    df = find_DataFrames(filenames)
    write_file = joinpath(dir, "DataFrames", filesavename)
    @info "Writing single DataFrame as $(write_file)"
    JLD2.save_object(write_file, df)
    return df
end

@time convertToDataFrames(datadir(), "unit_cratio_sweep.jld2")
