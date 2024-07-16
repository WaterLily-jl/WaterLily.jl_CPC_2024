#!/bin/bash

#SBATCH --job-name=waterlily_benchmarks
#SBATCH --mail-user=bernat.font@bsc.es
#SBATCH --mail-type=all
#SBATCH --account=bsc21
###SBATCH --qos=debug
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00


ml load rocm julia
export JULIA_DEPOT_PATH="/home/bsc/bsc021850/.julia-amd"
export WATERLILY_DIR="/gpfs/projects/bsc21/bsc021850/WaterLily/WaterLily.jl"

sh benchmark.sh -b "ROCArray" \
	-c "tgv sphere cylinder donut" \
	-p "6,7,8 3,4,5 4,5,6 5,6,7" \
	-s "100 100 100 100" \
	-ft "Float32 Float32 Float32 Float32"
