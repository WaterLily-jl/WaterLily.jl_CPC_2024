# This script should be run using run_benchmarks.sh, eq:
# $ sh $WATERLILY_DIR/benchmark/benchmark.sh -b "Array CuArray" -t "32 64" -c "tgv" -p "7,8" -s "100" -ft "Float32"

using WaterLilyBenchmarks

# using CUDA: CuArray, allowscalar
# allowscalar(false)

backend_str = Dict(Array => "CPUx$(Threads.nthreads())")#, CuArray => "GPU")
cases, log2p, max_steps, ftype, backend = parse_cla(ARGS;
    cases=["tgv", "jelly"], log2p=[(6,7), (5,6)], max_steps=[100, 100], ftype=[Float32, Float32], backend=Array
)

# Generate benchmark data
run_benchmarks(cases, log2p, max_steps, ftype, backend, backend_str[backend])

# Plot

