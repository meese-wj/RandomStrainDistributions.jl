#!/usr/bin/bash -l
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10g
#SBATCH --mail-type=all
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --array=1-10
#SBATCH --job-name=L-128
#SBATCH -o %x-%A_%a.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK distribution_generation.jl
    exit
=#

include("src/StrainSurveys.jl")
include("src/jobname_parser.jl")
using Random

Random.seed!(42)

ndis_vals = 2 .^ UnitRange(1, ENV["SLURM_ARRAY_TASK_COUNT"])

myLx = jobname_parser( ENV["SLURM_JOB_NAME"], "L", Int )
mynumber = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])

params = DistributionParameters(; Lx = myLx, ndislocations = mynumber, nsamples = 1e3 |> Int)

@time main( params; save_output = true )




