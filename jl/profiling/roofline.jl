# Post-processing Nsight Compute results to plot roofline model for different kernels

using Plots


fontsize = 16
speedup_fontsize = 16
Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 1,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 5),
    right_margin = Plots.Measures.Length(:mm, 10),
    # bottom_margin = Plots.Measures.Length(:mm, 5),
    # top_margin = Plots.Measures.Length(:mm, 5),
    legendfontsize = fontsize-4,
    tickfontsize = fontsize,
    labelfontsize = fontsize,
)

const mem_bandwidth = 256e9
const max_perf_fp32 = 7004e9
const min_AI_fp32 = 27.38
const max_perf_fp64 = 109e9
const min_AI_fp64 = 0.43
const x_lims = (1e-2, 1e4)
const y_lims = (1e9, 1e14)

y_mem(x) = mem_bandwidth*x/1e12

kern_395_fp32 = (0.54, 128052e6/1e12)
kern_451_fp32 = (0.61, 81288e6/1e12)
kern_451_fp64 = (0.15, 20322e6/1e12)

p = plot(ylabel="Performance [TFLOP/s]", xlabel="Arithmetic Intensity [FLOP/byte]",
    ylims=y_lims./1e12, xlims=x_lims, yminorgrid=true, xminorgrid=true,
    yaxis=:log10, xaxis=:log10, legend=:bottomright, background_color_legend=RGBA{Float64}(1, 1, 1, 0.7)
)
# fancylogscale!(p)
plot!(p, y_mem, linewidth=2, color=:darkblue, primary=false)
plot!(p, [min_AI_fp32, x_lims[2]], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:black, primary=false)
plot!(p, [x_lims[1], min_AI_fp32], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:black, linestyle=:dash, primary=false)
plot!(p, [min_AI_fp64, x_lims[2]], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, primary=false)
plot!(p, [x_lims[1], min_AI_fp64], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, linestyle=:dash, primary=false)

scatter!(p, kern_395_fp32, ms=5, ma=1, color=:lightgreen, label="kernel_395 FP32", markershape=:diamond)
scatter!(p, kern_451_fp32, ms=5, ma=1, color=:lightblue, label="kernel_451 FP32", markershape=:diamond)
scatter!(p, kern_451_fp64, ms=5, ma=1, color=:lightblue, label="kernel_451 FP64")
savefig(p, joinpath(string(@__FILE__), "../../../tex/img/roofline.pdf"))