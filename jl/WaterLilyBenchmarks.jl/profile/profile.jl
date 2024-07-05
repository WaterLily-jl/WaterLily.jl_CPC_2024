# Run profiling with (eg): nsys profile -o "tgv" julia --project profile.jl --case="tgv" --log2p=8
# Analyse as (eg): nsys stats -r nvtx_startend_sum --force-export=true "tgv.nsys-rep"

include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using WaterLily
using CairoMakie
using Printf
using CUDA
CUDA.allowscalar(false)

iarg(arg) = occursin.(arg, ARGS) |> findfirst
arg_value(arg) = split(ARGS[iarg(arg)], "=")[end]
metaparse(x) = eval(Meta.parse(x))
getf(str) = eval(Symbol(str))

function main(sim, max_steps)
    for _ in 1:10 sim_step!(sim; remeasure=false) end # compilation and warmup
    for i in 1:max_steps sim_step_profile!(sim; remeasure=false) end
end

const max_steps = 300
const T = Float32
const B = CuArray

# if "--case" in ARGS, run profiling
if !isnothing(iarg("case"))
    case = arg_value("case")
    log2p = !isnothing(iarg("log2p")) ? arg_value("log2p") |> metaparse : log2p
    max_steps = !isnothing(iarg("max_steps")) ? arg_value("max_steps") |> metaparse : max_steps
    ftype = !isnothing(iarg("ftype")) ? arg_value("ftype") |> metaparse : T
    backend = !isnothing(iarg("backend")) ? arg_value("backend") |> x -> eval(Symbol(x)) : B

    sim = getf(case)(log2p, backend; T=ftype)
    main(sim, max_steps)
# else parse and plot Nsys profile results
else
    nsys_fields = ["time_pc", "time_total", "instances", "time_avg", "time_med", "time_min", "time_max", "time_std", "range"]
    kernels = ["project!", "CFL!", "BDIM!", "BC!", "conv_diff!", "scale_u!", "BCTuple", "accelerate!", "exitBC!", "measure!"]
    kernel_instances_per_dt = Float64[2, 1, 2, 4, 2, 2, 2, 2, 1, 1]
    kernels_dict = Dict(k=>Dict{String,Any}("ipdt"=>kernel_instances_per_dt[i]) for (i,k) in enumerate(kernels))
    data = readlines(`nsys stats -r nvtx_startend_sum --force-export=true "data/tgv/tgv.nsys-rep"`)[8:end-1] .|> x->split(x,' '; keepempty=false)
    for kernel in data
        kernel_name = split(kernel[end],':')[end]
        for (i,k) in enumerate(nsys_fields[1:end-1])
            kernels_dict[kernel_name][k] = replace(kernel[i],","=>"") |> x->parse(Float64,x)
        end
    end
end

nsys_fields = ["time_pc", "time_total", "instances", "time_avg", "time_med", "time_min", "time_max", "time_std", "range"]
kernels = ["project!", "CFL!", "BDIM!", "BC!", "conv_diff!", "scale_u!", "BCTuple", "accelerate!", "exitBC!", "measure!"]
kernel_instances_per_dt = Float64[2, 1, 2, 4, 2, 2, 2, 2, 1, 1]
kernels_dict = Dict(k=>Dict{String,Any}("ipdt"=>kernel_instances_per_dt[i]) for (i,k) in enumerate(kernels))
println(read(`nsys stats -r nvtx_startend_sum --force-export=true "data/tgv/tgv.nsys-rep"`,String))
data = readlines(`nsys stats -r nvtx_startend_sum --force-export=true "data/tgv/tgv.nsys-rep"`)[8:end-1] .|> x->split(x,' '; keepempty=false)
for kernel in data
    kernel_name = split(kernel[end],':')[end]
    for (i,k) in enumerate(nsys_fields[1:end-1])
        kernels_dict[kernel_name][k] = replace(kernel[i],","=>"") |> x->parse(Float64,x)
    end
    kernels_dict[kernel_name]["time_weighted"] = kernels_dict[kernel_name]["time_med"] * kernels_dict[kernel_name]["ipdt"]
end


labels, kernel_weighted_time, i = [],Float64[], 0
bc_labels = ["BDIM!", "BCTuple", "exitBC!", "BC!"]
for (k,v) in kernels_dict
    if k in bc_labels && i==0
        push!(labels, "BC!")
        push!(kernel_weighted_time, v["time_weighted"])
        global i = length(labels)
    elseif k in bc_labels
        kernel_weighted_time[i] += v["time_weighted"]
    else
        push!(labels, k)
        push!(kernel_weighted_time, v["time_weighted"])
    end
end
sortidx = sortperm(lowercase.(labels))
labels = labels[sortidx]
kernel_weighted_time = kernel_weighted_time[sortidx]
colors = [:yellow, :orange, :red, :blue, :purple, :green, :cyan]
with_theme(theme_latexfonts()) do
fig, ax, plt = pie(
    kernel_weighted_time,
    color = colors,
    radius = 1,
    inner_radius = 0.6 ,
    strokecolor = :white,
    strokewidth = 1,
    axis = (aspect = AxisAspect(1), autolimitaspect = 1),
)
hidedecorations!(ax); hidespines!(ax)
nc = length(kernel_weighted_time)
cbar = Colorbar(fig[1,2], colormap = cgrad(colors, categorical = true))
cbar.ticks = (range(0+1/2nc, 1-1/2nc, nc), string.(labels))
data = mod.(kernel_weighted_time/sum(kernel_weighted_time),2π)
for (i, c) in enumerate(colors)
    θ = (sum(data[1:i-1]) + data[i]/2)*2π
    x = 0.5*cos(θ)
    y = 0.5*sin(θ)
    pc = kernel_weighted_time[i]/sum(kernel_weighted_time)*100
    pc > 2 && Makie.text!(x/0.6, y/0.6, text=@sprintf("%.0f", pc), color=:white, fontsize=20)
end
save("figure.pdf", fig,)
end

# savefig(p, "test.pdf")