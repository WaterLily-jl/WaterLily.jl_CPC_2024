using WaterLily
using CUDA
using DelimitedFiles
using Plots, StatsPlots, LaTeXStrings, CategoricalArrays, Printf, ColorSchemes
using JLD2

Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 3,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 5),
    right_margin = Plots.Measures.Length(:mm, 5),
    bottom_margin = Plots.Measures.Length(:mm, 5),
    top_margin = Plots.Measures.Length(:mm, 5),
    titlefontsize = 20,
    legendfontsize = 20,
    tickfontsize = 20,
    labelfontsize = 20,
)
linewidth = 3

function tgv(p, backend; Re=1600, T=Float32)
    L = 2^p; U = 1; κ=π/L; ν = 1/(κ*Re)
    function uλ(i,xyz)
        x,y,z = @. (xyz-0.0)/L*π
        i==1 && return -U*sin(x)*cos(y)*cos(z)
        i==2 && return  U*cos(x)*sin(y)*cos(z)
        return 0.
    end
    Simulation((L, L, L), (0, 0, 0), 1/κ; U=U, uλ=uλ, ν=ν, T=T, mem=backend)
end

function EZ!(σ, u, ν, L)
    Ω = prod(size(inside(σ)))
    @inside σ[I] = WaterLily.ke(I, u)
    KE = mapreduce(identity,+,@inbounds(σ[inside(σ)]))
    @inside σ[I] = WaterLily.ω_mag(I, u)
    Z = mapreduce(abs2,+,@inbounds(σ[inside(σ)]))
    return KE/Ω, Z*ν*L/Ω
end

function main(p, backend; Re=1600, T=Float32, t_max=20.0, verbose=true)
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

ps = [6,7,8,9]
Re = 1600
t_max = 20.0
T = Float64
backend = CuArray
run = false
datadir = "data/tgv"
pdf_file = "../../tex/img/tgv_validation.pdf"

data_dns = readdlm("data/tgv/TGV_Re1600.dat", skipstart=43)
t_dns, E_dns, Z_dns = data_dns[:,1], data_dns[:,2], data_dns[:,3]

function main()
    if run
        mkpath(datadir)
        for p in ps
            E, Z, t = main(p, backend; Re=Re, T=T, t_max=t_max)
            jldsave(joinpath(datadir,"p$p.jld2"); E=E, Z=Z, t=t)
        end
    end

    p1 = plot(t_dns, E_dns, label=" DNS ", color=:black, linewidth=linewidth)
    p2 = plot(t_dns, Z_dns, label=" DNS ", color=:black, linewidth=linewidth)
    for p in ps
        E, Z, t = jldopen(joinpath(datadir,"p$p.jld2"))["E"], jldopen(joinpath(datadir,"p$p.jld2"))["Z"], jldopen(joinpath(datadir,"p$p.jld2"))["t"]
        plot!(p1, t, E, label=L"~%$(2^p)^3", linewidth=linewidth, linestyle=:dash)
        plot!(p1, xlabel=L"$t\pi/L$", ylabel="Kinetic energy", framestyle=:box, grid=true, size=(1200, 600), ylims=(0, 0.15), xlims=(0, 20))
        plot!(p2, t, Z, label=L"~%$(2^p)^3\,", linewidth=linewidth, linestyle=:dash)
        plot!(p2, xlabel=L"$t\pi/L$", ylabel="Dissipation", framestyle=:box, grid=true, size=(1200, 600), ylims=(0, 0.015), xlims=(0, 20), legend=false)
    end
    plot(p1, p2, layout=(1, 2))
    fig_path = joinpath(string(@__DIR__), pdf_file)
    println("Figure stored in $(fig_path)")
    savefig(fig_path)
end

main()