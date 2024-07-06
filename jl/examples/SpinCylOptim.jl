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

# Davidson minimizer (should be in Optim, but isn't yet)
using ForwardDiff: Dual,Tag,value,partials
using TypedTables
function davidson(f::F,ax::R,bx::R;tol=1e-2,verbose=true) where {F,R}
    function eval(f,x)
        T = typeof(Tag(f,R))
        fx = f(Dual{T}(x,one(R)))
        (x=x,f=value(fx),∂=partials(T,fx,1))
    end
    a,b,Δ = eval(f,ax),eval(f,bx),bx-ax
    points = [a,b]
    verbose && println("1: ",a,"\n2: ",b)
    while abs(Δ) > 2tol
        v = a.∂+b.∂-3*(b.f-a.f)/Δ; w = copysign(√(v^2-a.∂*b.∂),Δ)
        x = b.x-Δ*(b.∂+w-v)/(b.∂-a.∂+2w)
        (x+tol>max(a.x,b.x) || x-tol<min(a.x,b.x)) && (x = 0.5*(a.x+b.x))
        c = eval(f,x); push!(points,c)
        verbose && println(length(points),": ",c)
        c.f > max(a.f,b.f) && return ErrorException("Process failed")
        c.∂*a.∂ > 0 ? (a=c) : (b=c)
        Δ = b.x-a.x
    end
    (a.f < b.f ? a : b), points |> Table
end
opt,_=davidson(x->-log(x)/x,1.,10.,tol=1e-6); @assert opt.x ≈ exp(1) # test minimizer

# Optimize the scaled power
function cost(ξ,t_end=2;d_D=0.15,remeasure=false)
    println("ξ=",ξ)
    sim = drag_control_sim(ξ;d_D)
    sim_step!(sim,t_end;remeasure)
    -scaled_power(sim,d_D*ξ^3)
end
opt,points = davidson(cost,3.,8.)

# Plot optimization points
using Plots
points = filter(a->a.f<0,points) # remove drag cases
scatter(points.x, points.f/opt.f, xlabel="ξ",ylabel="ηₚ/max(ηₚ)", ylims=(0,1.1), label=nothing,
        zcolor = eachindex(points), markercolor=:Blues, colorbar_title="Iteration")
savefig("SpinOptim.png")

# Plot the scaled power history
# opt = (x = 6.252768481467935, f = -0.0007056811885623282, ∂ = 0.0956334509946055)
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
plot(hist3,label="ξ=3",xlabel="scaled time",ylabel="scaled power")
scatter!(filter(a->abs(a.t-2)<3e-3,hist3),c=1,label=nothing)
plot!(histOpt,label="ξ=6.253",ylims=(-0.4,0.1),c=2)
scatter!(filter(a->abs(a.t-2)<3e-3,histOpt),c=2,label=nothing)
savefig("SpinCylHist.png")

# Plot the flow & circles
include("TwoD_plots.jl")
circle(R,c,τ=2π) = reim([R*exp(im*θ) for θ ∈ range(0,τ,length=33)] .+ c)
ω = sim.flow.σ; @inside ω[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(ω[inside(ω)],shift=(-1,-1),clims=(-10,10),xlims=(100,500),border=:none,legend=false)
R,C = sim.L÷2, 2sim.L;
r = 0.15R; c = C+(1.1R+r)*exp(im*π/3);
addbody(circle(R,C,π)...,c=:lightgreen)
addbody(circle(r,c)...,c=:mediumorchid)
savefig("SpinCylFlood.png")