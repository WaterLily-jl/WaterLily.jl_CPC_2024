# Post-processing Nsight Compute results to plot roofline model for different kernels

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
bandwidth_pc(x::Tuple{T,T} where T<:Number) = x[2]/x[1]/mem_bandwidth*100
norm_point(x::Tuple{T,T} where T<:Number) = (x[1], x[2]/1e12)

kern_395_fp32_G1 = (0.62, 93064e6)
kern_451_fp32_G1 = (0.97, 68071e6)
kern_451_fp64_G1 = (0.24, 17018e6)
kern_395_fp32_G2 = (0.55, 116743e6)
kern_451_fp32_G2 = (0.73, 81190e6)
kern_451_fp64_G2 = (0.18, 20047e6)
kern_395_fp32_G3 = (0.54, 128052e6)
kern_451_fp32_G3 = (0.61, 81288e6)
kern_451_fp64_G3 = (0.15, 20322e6)


p = plot(ylabel="Performance [TFLOP/s]", xlabel="Arithmetic Intensity [FLOP/byte]",
    ylims=y_lims./1e12, xlims=x_lims, yminorgrid=true, xminorgrid=true,
    yaxis=:log10, xaxis=:log10, legend=:bottomright, background_color_legend=RGBA{Float64}(1, 1, 1, 1)
)
# fancylogscale!(p)
plot!(p, y_mem_normed, linewidth=2, color=:darkblue, primary=false)
plot!(p, [min_AI_fp32, x_lims[2]], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:black, primary=false)
plot!(p, [x_lims[1], min_AI_fp32], [max_perf_fp32/1e12, max_perf_fp32/1e12], linewidth=2, color=:black, linestyle=:dash, primary=false)
plot!(p, [min_AI_fp64, x_lims[2]], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, primary=false)
plot!(p, [x_lims[1], min_AI_fp64], [max_perf_fp64/1e12, max_perf_fp64/1e12], linewidth=2, color=:grey, linestyle=:dash, primary=false)

scatter!(p, norm_point(kern_395_fp32_G1), label=@sprintf("project!   0.44M FP32 %d", bandwidth_pc(kern_395_fp32_G1))*"%",
    ms=markersize, ma=1, color=:yellow, markershape=:diamond, alpha=fc)
scatter!(p, norm_point(kern_395_fp32_G2), label=@sprintf("project!   3.54M FP32 %d", bandwidth_pc(kern_395_fp32_G2))*"%",
    ms=markersize, ma=1, color=:lightgreen, markershape=:diamond, alpha=fc)
scatter!(p, norm_point(kern_395_fp32_G3), label=@sprintf("project!  28.31M FP32 %d", bandwidth_pc(kern_395_fp32_G3))*"%",
    ms=markersize, ma=1, color=:green, markershape=:diamond, alpha=fc)

scatter!(p, norm_point(kern_451_fp32_G1), label=@sprintf("convdiff!  0.44M FP32 %d", bandwidth_pc(kern_451_fp32_G1))*"%",
    ms=markersize, ma=1, color=:berlin10, markershape=:diamond)
scatter!(p, norm_point(kern_451_fp32_G2), label=@sprintf("convdiff!  3.54M FP32 %d", bandwidth_pc(kern_451_fp32_G2))*"%",
    ms=markersize, ma=1, color=:cyan, markershape=:diamond)
scatter!(p, norm_point(kern_451_fp32_G3), label=@sprintf("convdiff! 28.31M FP32 %d", bandwidth_pc(kern_451_fp32_G3))*"%",
    ms=markersize, ma=1, color=:blue, markershape=:diamond)

scatter!(p, norm_point(kern_451_fp64_G1), label=@sprintf("convdiff!  0.44M FP64 %d", bandwidth_pc(kern_451_fp64_G1))*"%",
    ms=markersize, ma=1, color=:berlin10)
scatter!(p, norm_point(kern_451_fp64_G2), label=@sprintf("convdiff!  3.54M FP64 %d", bandwidth_pc(kern_451_fp64_G2))*"%",
    ms=markersize, ma=1, color=:cyan)
scatter!(p, norm_point(kern_451_fp64_G3), label=@sprintf("convdiff! 28.31M FP64 %d", bandwidth_pc(kern_451_fp64_G3))*"%",
    ms=markersize, ma=1, color=:blue)



savefig(p, joinpath(string(@__FILE__), "../../../tex/img/roofline_cylinder.pdf"))