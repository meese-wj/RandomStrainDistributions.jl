#!/usr/bin/bash -l
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --array=1-10
#SBATCH --job-name=L-112
#SBATCH -o %x-%A_%a.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK distribution_generation.jl
    exit
=#

using Pkg
Pkg.activate(@__DIR__)

include("src/StrainSurveys.jl")
include("src/jobname_parser.jl")
using Random

Random.seed!(42)

ndis_vals = 2 .^ UnitRange(1, parse(Int, ENV["SLURM_ARRAY_TASK_COUNT"]))
@show jobname = ENV["SLURM_JOB_NAME"]
@show jobid = ENV["SLURM_ARRAY_TASK_ID"]

nsamples = 2^13
myLx = jobname_parser( jobname, "L", Int; connector = "-" )
mynumber = ndis_vals[parse(Int, jobid)]

params = DistributionParameters(; Lx = myLx, ndislocations = mynumber, rtol = 0.001, nsamples = nsamples |> Int)

@time main( params; save_output = true )




