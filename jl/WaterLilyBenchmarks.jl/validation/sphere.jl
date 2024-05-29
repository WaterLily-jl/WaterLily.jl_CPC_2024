using StaticArrays
using CUDA
using WaterLily
using ReadVTK, WriteVTK
using JLD2
import Base: time

function sphere(p, backend; Re=3700, T=Float32)
    D = 2^p; U = 1; ν = U*D/Re
    L = (14D, 8D, 8D)
    center = SA[4D, 4D, 4D]
    body = AutoBody((x,t)-> √sum(abs2, x .- center) - D/2)
    Simulation(L, (U, 0, 0), D; U=U, ν=ν, body=body, T=T, mem=backend, exitBC=true)
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
function load!(meanflow::MeanFlow, fname::String)
    obj = jldopen(fname)
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

backend = CuArray
T = Float32
Re = 3700
ps = [4,5,6]

stats = true
stats_turb = false
time_max = 300.0 # in CTU
stats_init = 100.0 # in CTU
stats_interval = 0.1 # in CTU
dump_interval = 5000 # in CTU
datadir = "data/sphere/"
fname_output = "meanflow"
verbose = true

function run(p, backend; Re=3700, T=Float32)
    sim = sphere(p, backend; Re=Re, T=T)
    sim_out = sphere(p, backend; Re=Re, T=T)
    meanflow = MeanFlow(sim.flow)
    force,time = Vector{T}[],T[] # force coefficients (and time)
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=verbose)
        # Force stats
        push!(force, WaterLily.∮nds(sim.flow.p,sim.flow.f,sim.body,sim_time(sim))/(0.5*sim.U^2*sim.L^2))
        push!(time, sim_time(sim))
        verbose && println("Cd=",round(force[end][1],digits=4))
        if WaterLily.sim_time(sim)%dump_interval < sim.flow.Δt[end]*sim.U/sim.L
            verbose && println("Writing forces")
            jldsave(datadir*"force_p$p.jld2"; force=force, time=time)
        end
        # Mean flow stats
        if stats && sim_time(sim) > stats_init
            length(meanflow.t) == 1 && reset!(meanflow; t_init=WaterLily.time(sim))
            verbose && println("Computing stats")
            update!(meanflow, sim.flow; stats_turb=stats_turb)
            if WaterLily.sim_time(sim)%dump_interval < sim.flow.Δt[end]*sim.U/sim.L
                verbose && println("Writing stats")
                write!(fname_output*"_p$p", meanflow; dir=datadir, vtk=true, sim=sim_out)
            end
        end
    end
    verbose && println("Writing forces")
    jldsave(datadir*"force_p$p.jld2"; force=force, time=time)
    verbose && println("Writing stats")
    write!(fname_output*"_p$p", meanflow; dir=datadir, vtk=true, sim=sim_out)
    wr = vtkWriter("flow_p$p"; dir=datadir)
    WaterLily.write!(wr, sim)
    close(wr)
    println("Done!")
    return sim, meanflow, force
end

mkpath(datadir)
for p in ps
    sim, meanflow, force = run(p, backend; Re=Re, T=T)
end