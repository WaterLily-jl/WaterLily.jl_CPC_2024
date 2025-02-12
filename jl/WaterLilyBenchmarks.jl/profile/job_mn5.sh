#!/bin/bash

#SBATCH --job-name=waterlily_profiling
#SBATCH --account=ehpc26
#SBATCH --qos=acc_ehpc
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00

export JULIA_DEPOT_PATH="$HOME/.julia-acc"

sh profile.sh \
	-b "CuArray" \
	-c "tgv sphere cylinder" \
	-p "8 5 6" \
	-r "1"

