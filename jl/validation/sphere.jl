using StaticArrays
using CUDA
using WaterLily
using ReadVTK, WriteVTK
using FFTW, Interpolations, JLD2, Plots, LaTeXStrings, Printf, DSP
import Base: time

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

function sphere(D, backend; L=(8,2,2), center=SA[2,1,1], Re=3700, T=Float32)
    U = 1; ν = U*D/Re
    body = AutoBody((x,t)-> √sum(abs2, x .- (center .* D)) - D/2)
    Simulation(L.*D, (U, 0, 0), D; U=U, ν=ν, body=body, T=T, mem=backend, exitBC=true)
end

struct MeanFlow{T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Mf<:AbstractArray{T}}
    P :: Sf # pressure scalar field
    U :: Vf # velocity vector field
    UU :: Mf # squared velocity tensor
    τ :: Mf # Reynolds stress tensor
    t :: Vector{T} # time
    function MeanFlow(flow::Flow{D,T}; t_init=0.0) where {D,T}
        f = typeof(flow.u).name.wrapper
        P = zeros(T, size(flow.p)) |> f
        U = zeros(T, size(flow.u)) |> f
        UU = zeros(T, size(flow.p)...,D,D) |> f
        τ = zeros(T, size(UU)) |> f
        new{T,typeof(P),typeof(U),typeof(UU)}(P,U,UU,τ,T[t_init])
    end
end
time(meanflow::MeanFlow) = meanflow.t[end]-meanflow.t[1]
function reset!(meanflow::MeanFlow; t_init=0.0)
    fill!(meanflow.P, 0); fill!(meanflow.U, 0); fill!(meanflow.UU, 0); fill!(meanflow.τ, 0)
    deleteat!(meanflow.t, collect(1:size(meanflow.t)[1]))
    push!(meanflow.t, t_init)
end
function load!(meanflow::MeanFlow, fname::String; dir="data/")
    obj = jldopen(dir*fname)
    f = typeof(meanflow.U).name.wrapper
    meanflow.P .= obj["P"] |> f
    meanflow.U .= obj["U"] |> f
    meanflow.UU .= obj["UU"] |> f
    meanflow.τ .= obj["τ"] |> f
    meanflow.t .= obj["t"]
end
function write!(fname, meanflow::MeanFlow; dir="data/", vtk=false, sim=nothing)
    jldsave(dir*fname*".jld2";
        P=Array(meanflow.P),
        U=Array(meanflow.U),
        UU=Array(meanflow.UU),
        τ=Array(meanflow.τ),
        t=meanflow.t)
    if vtk && sim isa Simulation
        copy!(sim.flow, meanflow)
        wr = vtkWriter(fname; dir=dir)
        WaterLily.write!(wr, sim)
        close(wr)
    end
end
function update!(meanflow::MeanFlow, flow::Flow; stats_turb=true)
    dt = WaterLily.time(flow) - meanflow.t[end]
    ε = dt / (dt + (meanflow.t[end] - meanflow.t[1]) + eps(eltype(flow.p)))
    WaterLily.@loop meanflow.P[I] = ε*flow.p[I] + (1.0 - ε)*meanflow.P[I] over I in CartesianIndices(flow.p)
    WaterLily.@loop meanflow.U[Ii] = ε*flow.u[Ii] + (1.0 - ε)*meanflow.U[Ii] over Ii in CartesianIndices(flow.u)
    if stats_turb
        for i in 1:ndims(flow.p), j in 1:ndims(flow.p)
            WaterLily.@loop meanflow.UU[I,i,j] = ε*(flow.u[I,i].*flow.u[I,j]) + (1.0 - ε)*meanflow.UU[I,i,j] over I in CartesianIndices(flow.p)
            WaterLily.@loop meanflow.τ[I,i,j] = meanflow.UU[I,i,j] - meanflow.U[I,i,j]*meanflow.U[I,i,j] over I in CartesianIndices(flow.p)
        end
    end
    push!(meanflow.t, meanflow.t[end] + dt)
end
function copy!(a::Flow, b::MeanFlow)
    a.u .= b.U
    a.p .= b.P
end
function read_forces(fname::String; dir="data/")
    obj = jldopen(dir*fname)
    return obj["force"], obj["u_probe"], obj["time"]
end

function run_sim(D, backend; L=(7,3,3), center=SA[1.5,1.5,1.5], u_probe_loc=(4.5,2.1,1.5), u_probe_component=2, Re=3700, T=Float32, restart=false)
    sim = sphere(D, backend; L, center, Re, T)
    meanflow = MeanFlow(sim.flow)
    force,u_probe,time = Vector{T}[],T[] ,T[] # force coefficients, u probe location, time
    u_probe_loc_n = @. (u_probe_loc * sim.L) |> floor |> Int
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=false)
        # Force stats
        push!(force, WaterLily.total_force(sim)/(0.5*sim.U^2*sim.L^2))
        push!(u_probe, view(sim.flow.u,u_probe_loc_n...,u_probe_component) |> Array |> x->x[]) # WaterLily.interp(SA[7D,5D,4D], sim.flow.u[:,:,:,1]))
        push!(time, sim_time(sim))
        cd = round(force[end][1],digits=4)
        verbose && println("tU/D = $(time[end]); Cd = $cd")
        if WaterLily.sim_time(sim)%dump_interval < sim.flow.Δt[end]*sim.U/sim.L
            verbose && println("Writing force and probe values")
            jldsave(datadir*"force_D$D.jld2"; force=force, time=time, u_probe=u_probe)
        end
        # Mean flow stats
        if stats && sim_time(sim) > stats_init
            length(meanflow.t) == 1 && reset!(meanflow; t_init=WaterLily.time(sim))
            verbose && println("Computing stats")
            update!(meanflow, sim.flow; stats_turb=stats_turb)
            if WaterLily.sim_time(sim)%dump_interval < sim.flow.Δt[end]*sim.U/sim.L
                verbose && println("Writing stats")
                write!(fname_output*"_D$D", meanflow; dir=datadir)
            end
        end
    end
    verbose && println("Writing force and probe values")
    jldsave(datadir*"force_D$D.jld2"; force=force, u_probe=u_probe, time=time)
    verbose && println("Writing stats")
    write!(fname_output*"_D$D", meanflow; dir=datadir)
    wr = vtkWriter("flow_D$D"; dir=datadir)
    WaterLily.write!(wr, sim)
    close(wr)
    println("Done!")
    return sim, meanflow, force
end

backend = CuArray
T = Float32
Re = 3700
stats = true
stats_turb = false
time_max = 400.0 # in CTU
stats_init = 100.0 # in CTU
stats_interval = 0.1 # in CTU
dump_interval = 5000 # in CTU
Ds = [88,128,168] # diameter resolution
L = (7,3,3) # domain size in D # (7,3,3)
center = SA[1.5,1.5,1.5]
u_probe_loc = (4.5,2.1,1.5) # in D
u_probe_component = 2
datadir = "data/sphere/"
fname_output = "meanflow"
verbose = true
run = 0 # 0: postproc, 1: run
_plot = true

function main()
    run == 1 && mkpath(datadir)
    p_cd = plot()
    for D in Ds
        println("Running D = $D")
        if run == 1
            _, _, force = run_sim(D, backend; L, center, u_probe_loc, u_probe_component, Re, T)
        end
        # postproc forces
        force, u_probe, t = read_forces("force_D$D.jld2"; dir=datadir)
        force = mapreduce(permutedims, vcat, force)
        fx, fy, fz = force[:,1], force[:,2], force[:,3]
        idx = t .> stats_init
        fx, fy, fz, u, t = fx[idx], fy[idx], fz[idx], u_probe[idx], t[idx]
        CD_mean = sum(fx[2:end].*diff(t))/sum(diff(t))
        println("▷ ΔT [CTU] = "*@sprintf("%.4f", t[end]-t[1]))
        println("▷ CD_mean = "*@sprintf("%.4f", CD_mean))
        if _plot
            scatter!(p_cd, [D], [-CD_mean], grid=true, ms=10, ma=1, ylims=(0.25,0.4501), xlims=(80,176), xticks=Ds,
                xlabel=L"$D$", lw=0, framestyle=:box, size=(600, 600), legend=:bottomright, color=:black,
                legendfontsize=14, tickfontsize=18, labelfontsize=18, left_margin=Plots.Measures.Length(:mm, 5),
                ylabel=L"$\overline{C_D}$", label=D==168 ? "Present" : ""
            )
            # cd_plot = plot(t, -fx, linewidth=2, label=@sprintf("%.1f", prod(L.*D)/1e6)*" M")
            # plot!(cd_plot, xlabel=L"$tU/D$", ylabel=L"$C_D$", framestyle=:box, grid=true, size=(600, 600), ylims=(0.20, 0.40), xlims=(t[1], t[end]))
            # savefig(cd_plot, string(@__DIR__) * "../../../tex/img/sphere_D$(D)_CD.pdf")
        end
    end
    hline!(p_cd, [0.394], linestyle=:dash, color=:blue, label=L"\mathrm{Rodriguez}\,\,et\,\,al\mathrm{.\,\,(DNS)}")
    hline!(p_cd, [0.355], linestyle=:dashdot, color=:green, label=L"\mathrm{Yun}\,\,et\,\,al\mathrm{.\,\,(LES)}")
    savefig(p_cd, string(@__DIR__) * "../../../tex/img/sphere_validation.pdf")
end

main()
