# This script should be run using run_benchmarks.sh, eq:
# $ sh $WATERLILY_DIR/benchmark/benchmark.sh -b "Array CuArray" -t "32 64" -c "tgv" -p "7,8" -s "100" -ft "Float32"

include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using CUDA
using AMDGPU
using GPUArrays: allowscalar

allowscalar(false)

backend_str = Dict(Array => "CPUx$(Threads.nthreads())", CuArray => "GPU-NVIDIA", ROCArray => "GPU-AMD")
cases, log2p, max_steps, ftype, backend = parse_cla(ARGS;
    cases=["tgv", "jelly"], log2p=[(6,7), (5,6)], max_steps=[100, 100], ftype=[Float32, Float32], backend=Array
)

# Generate benchmark data
datadir = "./data/MN5_"*git_hash
mkpath(datadir)
run_benchmarks(cases, log2p, max_steps, ftype, backend, backend_str[backend]; datadir)