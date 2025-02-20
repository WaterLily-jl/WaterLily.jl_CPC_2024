#!/bin/bash
#SBATCH --job-name=waterlily_tgv
#SBATCH --account=ehpc26
#SBATCH --qos=acc_ehpc
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

export JULIA_DEPOT_PATH="$HOME/.julia-acc"

julia --project tgv.jl
