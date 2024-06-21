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
#SBATCH --time=15:00:00

ml load julia

export JULIA_DEPOT_PATH="/home/bsc/bsc021850/.julia-acc"

sh benchmark.sh -b "CuArray" \
	-c "tgv sphere cylinder donut" \
	-p "6,7,8,9 3,4,5,6 4,5,6,7 5,6,7,8" \
	-s "100 100 100 100" \
	-ft "Float32 Float32 Float32 Float32"
