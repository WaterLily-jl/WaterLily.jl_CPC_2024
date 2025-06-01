# Post-processing Nsight Compute results to plot roofline model for different kernels
# Hardcoded <arithmetic intensity [FLOP/byte], performance [FLOP/s]> values for cylinder case
# Values extracted from NCU (UI) roofline model analysis.
kern_project_fp32_p4 = (0.62, 93018e6)
kern_convdiff_fp32_p4 = (1.21, 98592e6)
kern_project_fp32_p5 = (0.55, 113463e6)
kern_convdiff_fp32_p5 = (0.92, 124483e6)
kern_project_fp32_p6 = (0.54, 128107e6)
kern_convdiff_fp32_p6 = (0.76, 132356e6)

using Plots, Printf
fontsize = 11
Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 1,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 5),
    right_margin = Plots.Measures.Length(:mm, 10),
    bottom_margin = Plots.Measures.Length(:mm, 5),
    top_margin = Plots.Measures.Length(:mm, 5),
    legendfontsize = fontsize-1,
    tickfontsize = fontsize,
    labelfontsize = fontsize,
)
markersize = 5
fc = RGBA(0.5, 0.5, 0.5, 0.5)

const mem_bandwidth = 256e9
const max_perf_fp32 = 7004e9
const min_AI_fp32 = 27.38
const max_perf_fp64 = 109e9
const min_AI_fp64 = 0.43
const x_lims = (1e-1, 1e2)
const y_lims = (1e9, 1e14)

y_mem_normed(x) = mem_bandwidth*x/1e12
bandwidth_pc(x) = x[2]/x[1]/mem_bandwidth*100
norm_point(x) = (x[1], x[2]/1e12)

p = plot(ylabel="Performance [TFLOP/s]", xlabel="Arithmetic Intensity [FLOP/byte]",
    ylims=y_lims./1e12, xlims=x_lims, yminorgrid=true, xminorgrid=true,
    yaxis=:log10, xaxis=:log10, legend=:bottomright, background_color_legend=RGBA{Float64}(1, 1, 1, 1)
)
# fancylogscale!(p)
# plot!(p, y_mem_normed, linewidth=2, color=:darkblue, primary=false)
plot!(p, [x_lims[1], min_AI_fp32], [y_mem_normed(x_lims[1]),y_mem_normed(min_AI_fp32)], linewidth=2, color=:black, primary=false)
plot!(p, [min_AI_fp32, x_lims[2]], [y_mem_normed(min_AI_fp32),y_mem_normed(x_lims[2])], linewidth=2, color=:black, linestyle=:dash, primary=false)
plot!(p, [x_lims[1], min_AI_fp32], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:grey, linestyle=:dash, primary=false)
plot!(p, [min_AI_fp32, x_lims[2]], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:grey, primary=false)
# plot!(p, [min_AI_fp64, x_lims[2]], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, primary=false)
# plot!(p, [x_lims[1], min_AI_fp64], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, linestyle=:dash, primary=false)

scatter!(p, norm_point(kern_project_fp32_p4), label=@sprintf("project!   0.44M FP32 %d", bandwidth_pc(kern_project_fp32_p4))*"%",
    ms=markersize, ma=1, color=:berlin10, markershape=:diamond, alpha=fc)
scatter!(p, norm_point(kern_project_fp32_p5), label=@sprintf("project!   3.54M FP32 %d", bandwidth_pc(kern_project_fp32_p5))*"%",
    ms=markersize, ma=1, color=:cyan, markershape=:diamond, alpha=fc)
scatter!(p, norm_point(kern_project_fp32_p6), label=@sprintf("project!  28.31M FP32 %d", bandwidth_pc(kern_project_fp32_p6))*"%",
    ms=markersize, ma=1, color=:blue, markershape=:diamond, alpha=fc)

scatter!(p, norm_point(kern_convdiff_fp32_p4), label=@sprintf("convdiff!  0.44M FP32 %d", bandwidth_pc(kern_convdiff_fp32_p4))*"%",
    ms=markersize, ma=1, color=:berlin10, markershape=:circle)
scatter!(p, norm_point(kern_convdiff_fp32_p5), label=@sprintf("convdiff!  3.54M FP32 %d", bandwidth_pc(kern_convdiff_fp32_p5))*"%",
    ms=markersize, ma=1, color=:cyan, markershape=:circle)
scatter!(p, norm_point(kern_convdiff_fp32_p6), label=@sprintf("convdiff! 28.31M FP32 %d", bandwidth_pc(kern_convdiff_fp32_p6))*"%",
    ms=markersize, ma=1, color=:blue, markershape=:circle)

savefig(p, joinpath(string(@__FILE__), "../../../tex/img/roofline_cylinder.pdf"))