#!/bin/bash
## Runs ../WaterLily-Benchmarks/profile.sh, post-process the data, and generates the plot for the TeX file.

THIS_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WB_DIR=$THIS_DIR"/../WaterLily-Benchmarks"

DATA_DIR=$THIS_DIR"/data"
PLOT_DIR=$THIS_DIR"/plots"
TEX_IMG_DIR=$THIS_DIR"/../../tex/img"

# Run profiling
sh $WB_DIR/profile.sh -c "tgv sphere cylinder" -p "8 5 6" -s 1000 -r 2 -dd $DATA_DIR -pd $PLOT_DIR

# Copy plots to TeX folder
cp $PLOT_DIR/* $TEX_IMG_DIR
