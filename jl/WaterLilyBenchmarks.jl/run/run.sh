#!/bin/bash

#SBATCH --job-name=waterlily_benchmarks
#SBATCH --mail-user=bernat.font@bsc.es
#SBATCH --mail-type=all
#SBATCH --account=bsc21
#SBATCH --queue=acc_bsccase
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

ml load julia

sh $WATERLILY_DIR/benchmark/benchmark.sh \
	-b "Array CuArray" -t "1 2 4 8 16 32 64" -c "tgv" -p "7,8" -s "100" -ft "Float32"
