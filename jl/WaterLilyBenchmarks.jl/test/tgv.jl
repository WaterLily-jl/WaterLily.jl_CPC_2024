using Revise
include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using WaterLily
using AMDGPU

function main(p, backend; Re=1600, T=Float32, t_max=1.0, verbose=true)
    sim = tgv(p, backend; Re=Re, T=T)
    while WaterLily.time(sim.flow)*(sim.U/sim.L) < t_max
        sim_step!(sim; remeasure=false)
        verbose && println("tU/L=", round(WaterLily.time(sim.flow)*(sim.U/sim.L), digits=4), ", Δt=",round(sim.flow.Δt[end], digits=3))
    end
end

main(5, ROCArray)
