#!/bin/bash
## Results from WaterLily-Benchmarks/benchmark.sh, which is run on MareNostrum5 using job_benchmark_mn5.sh, and on LUMI using job_benchmark_lumi.sh
## WaterLily-Benchmarks/compare.jl is used to post-process these results.

THIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WB_DIR=$THIS_DIR"/../WaterLily-Benchmarks"

DATA_DIR=$THIS_DIR"/data"
PLOT_DIR=$THIS_DIR"/plots"
TEX_IMG_DIR=$THIS_DIR"/../../tex/img"

# Post-process profiling
julia --project=$WB_DIR $WB_DIR/compare.jl --data_dir=$DATA_DIR --plot_dir=$PLOT_DIR --patterns="7e07b6b" \
    --speedup_base="CPUx01" --backend_color="lightrainbow" --sort=9
# Copy plots to TeX folder
cd $PLOT_DIR
for f in *.pdf; do cp "$f" "$TEX_IMG_DIR/$(echo "$f" | cut -f1,2 -d'_' | awk '{print $1".pdf"}')"; done
cd $THIS_DIR
