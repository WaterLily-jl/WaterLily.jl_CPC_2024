#!/bin/bash

#SBATCH --job-name=waterlily_benchmarks
#SBATCH --account=project_465001133
#SBATCH --partition=dev-g
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00

module use /appl/local/csc/modulefiles
module load julia
export WATERLILY_DIR="/users/fontbern/WaterLily.jl"

sh benchmark.sh -b "ROCArray" \
	-c "tgv sphere cylinder" \
	-p "6,7,8,9 3,4,5,6 4,5,6,7" \
	-s "100 100 100" \
	-ft "Float32 Float32 Float32"

