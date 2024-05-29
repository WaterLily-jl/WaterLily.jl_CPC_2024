#!/bin/bash
#SBATCH --job-name=waterlily_sphere
#SBATCH --mail-user=bernat.font@bsc.es
#SBATCH --mail-type=all
#SBATCH --account=bsc21
#SBATCH --qos=acc_bsccase
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00

ml load julia
export JULIA_DEPOT_PATH="/home/bsc/bsc021850/.julia-acc"

julia --project sphere.jl
