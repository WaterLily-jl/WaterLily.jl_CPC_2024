#!/bin/bash
#SBATCH --job-name=waterlily_sphere
#SBATCH --account=ehpc26
#SBATCH --qos=acc_ehpc
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00

export JULIA_DEPOT_PATH="$HOME/.julia-acc"

julia --proj sphere.jl --D=88
julia --proj sphere.jl --D=88 --Lc=6
julia --proj sphere.jl --D=128
julia --proj sphere.jl --D=168
julia --proj sphere.jl --run=false --plot=true
