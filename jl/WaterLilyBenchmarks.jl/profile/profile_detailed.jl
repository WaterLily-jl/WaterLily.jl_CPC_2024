# Run profiling with (eg): nsys profile -o "tgv" julia --project profile.jl --case="tgv" --log2p=8 --run
# Analyse as (eg): nsys stats -r nvtx_startend_sum --force-export=true "tgv.nsys-rep"

include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks
using WaterLily
using Plots: distinguishable_colors
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

# if "--run" in ARGS, run profiling
isnothing(iarg("case")) && @error "No case specified."
case = arg_value("case")
if metaparse(arg_value("run")) == 1
    log2p = !isnothing(iarg("log2p")) ? arg_value("log2p") |> metaparse : log2p
    max_steps = !isnothing(iarg("max_steps")) ? arg_value("max_steps") |> metaparse : max_steps_const
    ftype = !isnothing(iarg("ftype")) ? arg_value("ftype") |> metaparse : T_const
    backend = !isnothing(iarg("backend")) ? arg_value("backend") |> x -> eval(Symbol(x)) : B_const

    sim = getf(case)(log2p, backend; T=ftype)
    main(sim, max_steps; remeasure=case=="cylinder")
else
    # Postproc results
    nsys_fields = ["time_pc", "time_total", "instances", "time_avg", "time_med", "time_min", "time_max", "time_std", "range"]
    kernels = ["project!", "CFL!", "BDIM!", "BC!", "conv_diff!", "scale_u!", "BCTuple", "accelerate!", "exitBC!", "measure!",
        "residual!", "Jacobi!", "restrict!", "smooth!", "smooth2!", "Vcycle!", "norms"]
    kernel_instances_per_dt = Float64[2, 1, 2, 4, 2, 2, 2, 2, 1, 1]
    kernels_dict = Dict(k=>Dict{String,Any}() for k in kernels)
    println(read(`nsys stats -r nvtx_startend_sum --force-export=true "data/$case/$case.nsys-rep"`,String))
    data = readlines(`nsys stats -r nvtx_startend_sum --force-export=true "data/$case/$case.nsys-rep"`)[8:end-1] .|> x->split(x,' '; keepempty=false)
    length(data) < 5 && @error "Profiling was not successful."
    for kernel in data
        kernel_name = split(kernel[end],':')[end]
        for (i,k) in enumerate(nsys_fields[1:end-1])
            kernels_dict[kernel_name][k] = replace(kernel[i],","=>"") |> x->parse(Float64,x)
        end
        kernels_dict[kernel_name]["time_weighted"] = kernels_dict[kernel_name]["time_med"] # * kernels_dict[kernel_name]["ipdt"]
    end

    labels, kernel_weighted_time, i = [],Float64[], 0
    for (k,v) in kernels_dict
        push!(labels, k)
        push!(kernel_weighted_time, v["time_weighted"])
    end
    sortidx = sortperm(lowercase.(labels))
    labels = labels[sortidx]
    kernel_weighted_time = kernel_weighted_time[sortidx]
    colors = distinguishable_colors(length(kernels))

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
        fig.scene.theme[:figure_padding]=1
        hidedecorations!(ax); hidespines!(ax)
        # colorbar
        nc = length(kernel_weighted_time)
        cbar = Colorbar(fig[1,2], colormap = cgrad(colors, categorical = true))
        cbar.ticks = (range(0+1/2nc, 1-1/2nc, nc), string.(labels))
        data = mod.(kernel_weighted_time/sum(kernel_weighted_time),2π)
        for (i, c) in enumerate(colors)
            θ = (sum(data[1:i-1]) + data[i]/2)*2π
            x = 0.5*cos(θ)
            y = 0.5*sin(θ)
            pc = kernel_weighted_time[i]/sum(kernel_weighted_time)*100
            pc > 1.5 && Makie.text!(x/0.6, y/0.6, text=@sprintf("%.0f", pc), color=:white)
        end
        save(string(@__DIR__) * "../../../../tex/img/$(case)_profile_detailed.pdf", fig)
    end
end