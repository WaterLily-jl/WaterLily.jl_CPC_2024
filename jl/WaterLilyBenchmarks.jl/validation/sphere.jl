using Revise
include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using WaterLily
import WaterLily: CFL
using CUDA
using DelimitedFiles
using Plots, StatsPlots, LaTeXStrings, CategoricalArrays, Printf, ColorSchemes
using JLD2

Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 1,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 5),
    right_margin = Plots.Measures.Length(:mm, 5),
    bottom_margin = Plots.Measures.Length(:mm, 5),
    top_margin = Plots.Measures.Length(:mm, 5),
    titlefontsize = 15,
    legendfontsize = 14,
    tickfontsize = 14,
    labelfontsize = 14,
)
linewidth = 2

function main(p, backend; Re=3600, T=Float32, t_max=20.0, verbose=true)
    sim = tgv(p, backend; Re=Re, T=T)
    E0, Z0 = EZ!(sim.flow.σ, sim.flow.u, sim.flow.ν, sim.L)
    E, Z, t = T[E0], T[Z0], T[0.0]
    while WaterLily.time(sim.flow)*(sim.U/sim.L) < t_max
        sim_step!(sim; remeasure=false)
        verbose && println("tU/L=", round(WaterLily.time(sim.flow)*(sim.U/sim.L), digits=4), ", Δt=",round(sim.flow.Δt[end], digits=3))
        Ei, Zi = EZ!(sim.flow.σ, sim.flow.u, sim.flow.ν, sim.L)
        push!(E, Ei); push!(Z, Zi); push!(t, WaterLily.time(sim.flow)*(sim.U/sim.L))
    end
    return E, Z, t
end

ps = [5,7,8]
Re = 3600
t_max = 20.0
T = Float32
backend = CuArray
run = true

data = readdlm("data/tgv/TGV_Re1600.dat", skipstart=43)
t_dns, E_dns, Z_dns = data[:,1], data[:,2], data[:,3]

E, Z, t = main(ps[1], backend; Re=Re, T=T, t_max=t_max)