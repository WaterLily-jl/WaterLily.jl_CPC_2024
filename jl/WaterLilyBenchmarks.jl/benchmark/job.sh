#!/bin/bash

#SBATCH --job-name=waterlily_benchmarks
#SBATCH --mail-user=bernat.font@bsc.es
#SBATCH --mail-type=all
#SBATCH --account=bsc21
#SBATCH --qos=acc_bsccase
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00

ml load julia
export JULIA_DEPOT_PATH="/home/bsc/bsc021850/.julia-acc"
export WATERLILY_DIR="/gpfs/projects/bsc21/bsc021850/WaterLily/WaterLily.jl"

sh benchmark.sh -b "Array CuArray" \
	-t "1 2 4 8 16" \
	-c "tgv sphere cylinder" \
	-p "6,7,8,9 3,4,5,6 4,5,6,7" \
	-s "100 100 100" \
	-ft "Float32 Float32 Float32"
