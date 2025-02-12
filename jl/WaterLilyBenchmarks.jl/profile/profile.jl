# Run profiling with (eg): nsys profile -o "./data/tgv/tgv.nsys-rep" --force-overwrite=true --export=sqlite julia --project profile.jl --case="tgv" --log2p=8 --run=1
# Analyse as (eg): nsys stats -r nvtx_gpu_proj_sum "./data/tgv/tgv.sqlite"
# Profiling results is stored in data/. Data not included in the repository because it weights 1GB approximately, but it is available upon request.
# WaterLily#kernel_profiling2 branch must be used for traces.

# Analyse stats of the main kernels for the main ranges (project! and conv_diff!)
# ncu --set full --kernel-name gpu___kern__446 --launch-skip 10586 --launch-count 1 -o myreport julia --project=/home/b-fg/Documents/tudelft/documents/papers/journals/WaterLily.jl_CPC_2024/jl/WaterLilyBenchmarks.jl/profile --startup-file=no /home/b-fg/Documents/tudelft/documents/papers/journals/WaterLily.jl_CPC_2024/jl/WaterLilyBenchmarks.jl/profile/profile.jl --case=tgv --log2p=8 --backend=CuArray --max_steps=1000 --ftype=Float32 --run=1

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

function main(sim, max_steps; remeasure=false)
    for i in 1:max_steps sim_step!(sim; remeasure=remeasure) end # profiling
end

const max_steps_const = 1000
const T_const = Float32
const B_const = CuArray

# if "--run=1" in ARGS, run profiling
isnothing(iarg("case")) && @error "No case specified."
case = arg_value("case")
if metaparse(arg_value("run")) == 1 # run profiling
    log2p = !isnothing(iarg("log2p")) ? arg_value("log2p") |> metaparse : log2p
    max_steps = !isnothing(iarg("max_steps")) ? arg_value("max_steps") |> metaparse : max_steps_const
    ftype = !isnothing(iarg("ftype")) ? arg_value("ftype") |> metaparse : T_const
    backend = !isnothing(iarg("backend")) ? arg_value("backend") |> x -> eval(Symbol(x)) : B_const

    sim = getf(case)(log2p, backend; T=ftype)
    main(sim, max_steps; remeasure=case=="cylinder")
else # postprocess profiling
    nsys_fields = ["range", "style", "total_proj_time", "total_range_time", "instances", "proj_avg", "proj_median", "proj_min", "proj_max", "proj_std", "total_gpu_ops", "avg_gpu_ops", "avg_range_level", "avg_num_child"]
    kernels = ["project!", "CFL!", "BDIM!", "BC!", "conv_diff!", "scale_u!", "copy_u0!", "exitBC!", "measure!"] # "BCTuple", "accelerate!"
    kernel_instances_per_dt = Float64[2, 1, 2, 4, 2, 2, 2, 2, 1, 1]
    kernels_dict = Dict(k=>Dict{String,Any}("ipdt"=>kernel_instances_per_dt[i], "time_weighted"=>0.0) for (i,k) in enumerate(kernels))
    data = readlines(`nsys stats -r nvtx_gpu_proj_sum "data/$case/$case.sqlite"`)[7:end-1] .|> x->split(x,' '; keepempty=false)
    # println(read(`nsys stats -r nvtx_gpu_proj_sum "data/$case/$case.sqlite"`, String))
    length(data) < 5 && @error "Profiling was not successful."
    for kernel in data
        kernel_name = split(kernel[1],':')[end]
        for (i,k) in enumerate(nsys_fields[3:end])
            kernels_dict[kernel_name][k] = replace(kernel[i+2],","=>"") |> x->parse(Float64,x)
        end
        kernels_dict[kernel_name]["time_weighted"] = kernels_dict[kernel_name]["proj_median"] * kernels_dict[kernel_name]["ipdt"]
    end
    total_weighted_time = sum(v["time_weighted"] for (_,v) in kernels_dict)
    for (k,v) in kernels_dict
        kernels_dict[k]["time_weighted_pc"] = kernels_dict[k]["time_weighted"] / total_weighted_time
    end
    bc_labels = ["BDIM!", "exitBC!", "BC!"]
    kernels_dict["BCs!"] = Dict{String,Any}("time_weighted"=>0.0)
    for (k,v) in kernels_dict
        !(k in bc_labels) && continue
        kernels_dict["BCs!"]["time_weighted"] += kernels_dict[k]["time_weighted"]
    end
    kernels_dict["BCs!"]["time_weighted_pc"] = kernels_dict["BCs!"]["time_weighted"] / total_weighted_time

    labels, kernel_weighted_time, i = [], Float64[], 0
    for (k,v) in kernels_dict
        k in bc_labels && continue
        push!(labels, k)
        push!(kernel_weighted_time, v["time_weighted"])
    end
    sortidx = sortperm(lowercase.(labels), rev=true)
    labels = labels[sortidx]
    kernel_weighted_time = kernel_weighted_time[sortidx]
    colors = [:pink, :grey, :orange, :red, :blue, :purple, :green]

    with_theme(theme_latexfonts(), fontsize=25, figure_padding=1) do
        fig, ax, plt = pie(
            kernel_weighted_time,
            color = colors,
            radius = 1,
            inner_radius = 0.6 ,
            strokecolor = :white,
            strokewidth = 1,
            axis = (aspect = AxisAspect(1), autolimitaspect = 1),
        )
        fig.scene.theme[:figure_padding] = 1
        hidedecorations!(ax); hidespines!(ax)
        # colorbar
        if case == "cylinder"
            nc = length(kernel_weighted_time)
            cbar = Colorbar(fig[1,2], colormap=cgrad(colors, categorical=true))
            cbar.ticks = (range(0+1/2nc, 1-1/2nc, nc), string.(labels))
        end
        global data = mod.(kernel_weighted_time/sum(kernel_weighted_time),2π)
        for (i, c) in enumerate(colors)
            θ = (sum(data[1:i-1]) + 0.5*data[i])*2π
            x = 0.8*cos(θ)
            y = 0.8*sin(θ)
            pc = kernel_weighted_time[i]/sum(kernel_weighted_time)*100
            pc > 1.9 && Makie.text!(x, y, text=@sprintf("%.0f", pc), color=:white, align=(:center, :center))
        end
        fig_path = string(@__DIR__) * "../../../../tex/img/$(case)_profile.pdf"
        save(fig_path, fig)
        println("Figure stored in $(fig_path)")
    end
end