#!/bin/bash
## Runs ../WaterLily-Benchmarks/profile.sh, post-process the data, and generates the plot for the TeX file.

THIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WB_DIR=$THIS_DIR"/../WaterLily-Benchmarks"

DATA_DIR=$THIS_DIR"/data"
PLOT_DIR=$THIS_DIR"/plots"
TEX_IMG_DIR=$THIS_DIR"/../../tex/img"

# Run and post-process profiling of kernel timings distribution
sh $WB_DIR/profile.sh -c "tgv sphere cylinder" -p "8 5 6" -s 1000 -r 2 -dd $DATA_DIR -pd $PLOT_DIR
# Copy plots to TeX folder
cp $PLOT_DIR/* $TEX_IMG_DIR

# Run and post-process profiling of most expensive kernels (451, 395) for 3 different grid levels of the cylinder
sh $WB_DIR/profile.sh -ns "ncu" -k "327 355" -dd $DATA_DIR -c "cylinder" -p "4,5,6" -s 20
# Check <arithmetic intensity [FLOP/byte], performance [FLOP/s]> values from NCU UI roofline model and hardcode them in roofline.jl
# Run roofline.jl
julia --proj roofline.jl