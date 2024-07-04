using Revise
include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using WaterLily
import WaterLily: CFL
using CUDA

function main(p, backend; Re=1600, T=Float32, t_max=20.0, verbose=true)
    sim = tgv(p, backend; Re=Re, T=T)
    sim_step!(sim; remeasure=false)
    return sim
end

Re = 1600
t_max = 0.001
T = Float64
backend = Array

sim = main(9, backend; Re=Re, T=T, t_max=t_max)

