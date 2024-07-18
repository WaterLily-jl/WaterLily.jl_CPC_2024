# Run with
# julia --project compare.jl $(find data/MN5_36dc237 \( -name "tgv*json" -o -name "sphere*json" -o -name "cylinder*json" \) -printf "%T@ %Tc %p\n" | sort -n | awk '{print $7}') --sort=0

using BenchmarkTools, PrettyTables, Plots, StatsPlots, LaTeXStrings, CategoricalArrays, Printf, ColorSchemes
include("../src/WaterLilyBenchmarks.jl")
using .WaterLilyBenchmarks

# Utils
fontsize = 20
speedup_fontsize = 14
Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 1,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 24),
    right_margin = Plots.Measures.Length(:mm, 0),
    bottom_margin = Plots.Measures.Length(:mm, 5),
    top_margin = Plots.Measures.Length(:mm, 5),
    legendfontsize = fontsize,
    tickfontsize = fontsize,
    labelfontsize = fontsize,
)

function Base.unique(ctg::CategoricalArray)
    l = levels(ctg)
    newctg = CategoricalArray(l)
    levels!(newctg, l)
end

function annotated_groupedbar(xx, yy, group; series_annotations="", bar_width=1.0, plot_kwargs...)
    gp = groupedbar(xx, yy, group=group, series_annotations="", bar_width=bar_width; plot_kwargs...)
    m = length(unique(group))       # number of items per group
    n = length(unique(xx))          # number of groups
    xt = (1:n) .- 0.5               # plot x-coordinate of groups' centers
    dx = bar_width/m                # each group occupies bar_width units along x
    # dy = diff([extrema(yy)...])[1]
    x2 = [xt[i] + (j - m/2 - 0.3)*dx for j in 1:m, i in 1:n][:]
    k = 1
    for i in 1:n, j in 1:m
        y0 = gp[1][2j][:y][i]*1.3# + 0.04*dy
        if isfinite(y0)
            tex = series_annotations[k]
            annotate!(x2[(i-1)*m + j]*1.01, y0, text(series_annotations[k], :center, :black, speedup_fontsize))
            k += 1
        end
    end
    gp
end

tests_dets = Dict(
    "tgv" => Dict("size" => (1, 1, 1), "title" => "TGV"),
    "sphere" => Dict("size" => (16, 6, 6), "title" => "Sphere"),
    "cylinder" => Dict("size" => (12, 6, 2), "title" => "Moving cylinder"),
    "donut" => Dict("size" => (2, 1, 1), "title" => "Donut"),
)
tests_ordered = String["tgv", "sphere", "cylinder"]

# Load benchmarks
benchmarks_all = [BenchmarkTools.load(fname)[1] for fname in ARGS if !occursin("--sort", fname)]
sort_cla = findfirst(occursin.("--sort", ARGS))
sort_idx = !isnothing(sort_cla) ? ARGS[sort_cla] |> x -> split(x, "=")[end] |> x -> parse(Int, x) : 0

# Separate benchmarks by test case
cases_str = [b.tags[1] for b in benchmarks_all] |> unique
benchmarks_all_dict = Dict(Pair{String, Vector{BenchmarkGroup}}(k, []) for k in cases_str)
for b in benchmarks_all
    push!(benchmarks_all_dict[b.tags[1]], b)
end

# Table and plots
p_cost = plot()
for (kk, case) in enumerate(tests_ordered)
    benchmarks = benchmarks_all_dict[case]
    # Get backends string vector and assert same case sizes for the different backends
    backends_str = [String.(k)[1] for k in keys.(benchmarks)]
    log2p_str = [String.(keys(benchmarks[i][backend_str])) for (i, backend_str) in enumerate(backends_str)]
    @assert length(unique(log2p_str)) == 1
    log2p_str = sort(log2p_str[1])
    f_test = benchmarks[1].tags[2]
    # Get data for PrettyTables
    header = ["Backend", "WaterLily", "Julia", "Precision", "Allocations", "GC [%]", "Time [s]", "Cost [ns/DOF/dt]", "Speed-up"]
    data, base_speedup = Matrix{Any}(undef, length(benchmarks), length(header)), 1.0
    global data_plot = Array{Float64}(undef, length(log2p_str), length(backends_str), 3) # times, cost, speedups
    printstyled("Benchmark environment: $case $f_test (max_steps=$(benchmarks[1].tags[4]))\n", bold=true)
    for (k,n) in enumerate(log2p_str)
        printstyled("▶ log2p = $n\n", bold=true)
        for (i, benchmark) in enumerate(benchmarks)
            datap = benchmark[backends_str[i]][n][f_test]
            speedup = i == 1 ? 1.0 : benchmarks[1][backends_str[1]][n][f_test].times[1] / datap.times[1]
            N = prod(tests_dets[case]["size"]) .* 2 .^ (3 .* eval(Meta.parse.(n)))
            cost = datap.times[1] / N / benchmarks[1].tags[4]
            data[i, :] .= [backends_str[i], benchmark.tags[end-1], benchmark.tags[end], benchmark.tags[end-3],
                datap.allocs, (datap.gctimes[1] / datap.times[1]) * 100.0, datap.times[1] / 1e9, cost, speedup]
        end
        sorted_cond, sorted_idx = 0 < sort_idx <= length(header), nothing
        if sorted_cond
            sorted_idx = sortperm(data[:, sort_idx])
            baseline_idx = findfirst(x->x==1, sorted_idx)
            data .= data[sorted_idx, :]
        end
        hl_base = Highlighter(f=(data, i, j) -> sorted_cond ? i == findfirst(x->x==1, sorted_idx) : i==1,
            crayon=Crayon(foreground=:blue))
        hl_fast = Highlighter(f=(data, i, j) -> i == argmin(data[:, end-1]), crayon=Crayon(foreground=(32,125,56)))
        pretty_table(data; header=header, header_alignment=:c, highlighters=(hl_base, hl_fast), formatters=ft_printf("%.2f", [6,7,8,9]))
        data_plot[k, :, :] .= data[:, end-2:end]
    end
    # Cost plot in GPU
    N = prod(tests_dets[case]["size"]) .* 2 .^ (3 .* eval(Meta.parse.(log2p_str)))
    N_str = (N./1e6) .|> x -> @sprintf("%.2f", x)
    scatter!(p_cost, N./1e6, data_plot[:,end,2], yaxis=:log10, xaxis=:log10, yminorgrid=true, xminorgrid=true,
        ms=10, ma=1, ylims=(1,30), xlims=(0.1,400),
        label=tests_dets[case]["title"], xlabel="DOF [M]", lw=0, framestyle=:box, grid=:xy, size=(600, 600), legend=true,
        legendfontsize=15, tickfontsize=18, labelfontsize=18, left_margin=Plots.Measures.Length(:mm, 5),
        ylabel="Cost [ns/DOF/dt]"
    )

    # Speedup plot
    N = prod(tests_dets[case]["size"]) .* 2 .^ (3 .* eval(Meta.parse.(log2p_str)))
    N_str = (N./1e6) .|> x -> @sprintf("%.2f", x)
    groups = repeat(N_str, inner=length(backends_str)) |> CategoricalArray
    levels!(groups, N_str)
    ctg = repeat(backends_str, outer=length(log2p_str)) |> CategoricalArray
    levels!(ctg, backends_str)
    p = annotated_groupedbar(groups, transpose(data_plot[:,:,1]), ctg;
        series_annotations=vec(transpose(data_plot[:,:,3])) .|> x -> @sprintf("%d", x) .|> latexstring, bar_width=0.92,
        Dict(:xlabel=>"DOF [M]", # :title=>tests_dets[case]["title"],
            :ylims=>(1e-1, 1e5), :lw=>0, :framestyle=>:box, :yaxis=>:log10, :grid=>true,
            :color=>reshape(palette([:cyan, :green], length(backends_str))[1:length(backends_str)], (1, length(backends_str))),
            :size=>(600, 600)
        )...
    )
    if kk == 1
        plot!(p, ylabel="Time [s]", legend=:topleft, left_margin=Plots.Measures.Length(:mm, 0))
    else
        plot!(p, ylabel="", legend=false, yformatter=Returns(""))
    end
    savefig(p, string(@__DIR__)*"../../../../tex/img/$(case)_benchmark.pdf")
end
fancylogscale!(p_cost)
savefig(p_cost, string(@__DIR__)*"../../../../tex/img/cost.pdf")

