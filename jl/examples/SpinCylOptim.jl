using WaterLily,StaticArrays

# Spinning cylinder control simulation
using WaterLily: norm2                   # vector length
rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)] # rotation matrix
function drag_control_sim(ξ;D=96,Re=500,d_D=0.15f0,g_D=0.05f0)
    # set up big cylinder
    C,R,U = SA[2D,0],D÷2,1
    big = AutoBody((x,t)->norm2(x-C)-R)
    # set up small control cylinder
    r = d_D*R; c = C+(R+r+g_D*D)*SA{Float32}[1/2,√3/2]
    small = AutoBody((x,t)->norm2(x)-r,(x,t)->rot(ξ[]*U*t/r)*(x-c))
    # set up simulation
    Simulation((6D,2D),(U,0),D;ν=U*D/Re,body=big+small,T=typeof(ξ[]))
end

# Measure propulsive power (Thrust*U) normalized by the propulsor powering scale dξ³U³
scaled_power(sim::Simulation,dξ³) = WaterLily.total_force(sim)[1]*sim.U/(dξ³*sim.U^3*sim.L)

# Optimize the scaled power
include("Davidson.jl")
function cost(ξ,t_end=2;d_D=0.15,remeasure=false)
    println("ξ=",ξ)
    sim = drag_control_sim(ξ;d_D)
    sim_step!(sim,t_end;remeasure)
    -scaled_power(sim,d_D*ξ^3)
end
points = davidson(cost,3.,8.,tol=5e-2,verbose=true)

# Plot optimization points
using Plots,LaTeXStrings

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
    colorbar_titlefontsize = 14,
    tickfontsize = 14,
    labelfontsize = 14,
)

points = filter(a->a.f<0,points) # remove drag cases\
mine = palette([:cyan, :blue], 5)
scatter(points.x, points.f/opt.f, xlabel=L"$\xi$",ylabel=L"$c_p/\max(c_p)$", ylims=(0,1.1), label=nothing,
        zcolor = eachindex(points), markercolor=mine, colorbar_title="Iteration", size=(600,600))
savefig("SpinOptim.pdf")

# Plot the scaled power history
opt = (x = 6.252768481467935, f = -0.0007056811885623282, ∂ = 0.0956334509946055)
function history(ξ,duration=3;d_D=0.15,remeasure=false)
    sim = drag_control_sim(ξ;d_D)
    data = map(range(0,duration,300)) do tᵢ
        sim_step!(sim,tᵢ;remeasure)
        (t=sim_time(sim),cₚ=scaled_power(sim,d_D*ξ^2))
    end |> Table
    data,sim
end
histOpt,sim = history(opt.x);
hist3,_ = history(3.0);
plot(hist3,label=L"$\xi=3$")
scatter!(filter(a->abs(a.t-2)<3e-3,hist3),c=1,label=nothing)
plot!(histOpt,label=L"$\xi=6.253$",ylims=(-0.4,0.1),c=2)
scatter!(filter(a->abs(a.t-2)<3e-3,histOpt),c=2,label=nothing)
plot!(xlabel=L"$tU/D$",ylabel=L"$c_p$",size=(600,600))
savefig("SpinCylHist.pdf")

# Plot the flow & circles
include("TwoD_plots.jl")
circle(R,c,τ=2π) = reim([R*exp(im*θ) for θ ∈ range(0,τ,length=33)] .+ c)
ω = sim.flow.σ; @inside ω[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(ω[inside(ω)],shift=(-1,-1),clims=(-10,10),xlims=(100,500),border=:none,legend=false)
R,C = sim.L÷2, 2sim.L;
r = 0.15R; c = C+(1.1R+r)*exp(im*π/3);
addbody(circle(R,C,π)...,c=:lightgreen)
addbody(circle(r,c)...,c=:mediumorchid)
savefig("SpinCylFlood.pdf")