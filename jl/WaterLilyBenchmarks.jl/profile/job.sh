#!/bin/bash

#SBATCH --job-name=waterlily_profiling
#SBATCH --account=bsc21
#SBATCH --qos=acc_debug
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00

ml load julia
export JULIA_DEPOT_PATH="/home/bsc/bsc021850/.julia-acc"
export WATERLILY_DIR="/gpfs/projects/bsc21/bsc021850/WaterLily/WaterLily.jl"

sh profile.sh -c "tgv sphere cylinder" -p "8 5 6"
