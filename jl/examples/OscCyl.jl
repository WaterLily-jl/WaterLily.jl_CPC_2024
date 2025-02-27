using WaterLily, StaticArrays, CUDA, TypedTables, ForwardDiff, Random, JLD2
Random.seed!(1234)

function OscCyl(;Re=4000,U=1,Vr=5,Ay_D=0.4,Ax_D=0.0,θ=0.,    # physical
                n=2^7,T=Float32,mem=Array,                   # numerical
                t₀=4Vr,duration=4Vr,step=Vr/100,verbose=true, # printing
                run=true)
    # Set up sizes, frequency and mapping
    D,H = n/6,n; ω = 2π*U/(Vr*D)
    verbose && println("n=$n, D=$D, H=$H, ω=$ω")
    sdf(x,t) = hypot(x[1]-n,x[2]-n) - D/2
    map(x,t) = x - SA[Ax_D*D*sin(2ω*t + θ), Ay_D*D*sin(ω*t), 0]

    # Make Simulation and run first convective length
    sim = Simulation((4n,2n,n),(U,0,0),D;ν=U*D/Re,body=AutoBody(sdf,map),T,mem,exitBC=true,perdir=())
    !run && return nothing, sim
    sim_step!(sim,1)
    verbose && println("2D: t=1, Cf=$(2WaterLily.total_force(sim)/(D*H*U^2))")

    # Perturb into 3D flow, checking and running to t=t₀
    a = sim.flow; R = WaterLily.inside_u(size(a.p))
    a.u[R] += 0.1*(0.5 .- rand(T,size(R))|>mem)
    WaterLily.BC!(a.u,a.U,a.exitBC,a.perdir)
    for t in 2:t₀-step
        sim_step!(sim,t)
        verbose && println("3D?: t=$t, Cf=$(2WaterLily.total_force(sim)/(D*H*U^2))")
    end

    # Run and save force and motion
    function get_force!(sim,tᵢ)
        sim_step!(sim,tᵢ;remeasure=true)
        Cf = 2WaterLily.total_force(sim)/(D*H*U^2)
        t = WaterLily.time(sim)
        x = -map([0,0,0],t)
        u = ForwardDiff.derivative(t->-map([0,0,0],t),t)
        verbose && println("tᵢ=$tᵢ, Cf=$Cf")
        (t=t*U/D,x=x[1],y=x[2],u=u[1],v=u[2],Cx=Cf[1],Cy=Cf[2])
    end
    data=[get_force!(sim,tᵢ) for tᵢ ∈ t₀:step:t₀+duration] |> Table
    data,sim
end
function read!(a::AbstractArray, fname; dir="./")
    obj = jldopen(joinpath(dir, fname*".jld2"))
    f = typeof(a).name.wrapper
    a .= obj["a"] |> f
end
function read(fname; dir="./")
    obj = jldopen(joinpath(dir, fname*".jld2"))
    obj["a"] |> Table
end
write!(fname, a::AbstractArray; dir="./") = jldsave(
    joinpath(dir, fname*".jld2");
    a=Array(a)
)
write!(fname, a; dir="./") = jldsave(
    joinpath(dir, fname*".jld2");
    a=a
)

p = 7
Vr = 5.4
run = true
data, sim = OscCyl(n=2^p,Re=7620,Ay_D=1.6,Ax_D=0.4,θ=π/6,Vr=Vr,mem=CUDA.CuArray,run=run,t₀=8Vr,duration=8Vr);
if run
    write!("data_p$p", data)
    write!("u_p$p", sim.flow.u)
else
    data = read("data_p$p")
    read!(sim.flow.u, "u_p$p")
end

# using CSV
# CSV.write("OscCyl_Ay1p6_Ax0p4_Vr5p4_T30.csv",data)

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
Cpow(p) = p.Cx*p.u+p.Cy*p.v
mCpow(data) = sum(Cpow(data[i])*(data.t[i+1]-data.t[i]) for i in eachindex(data)[1:end-1])/(data.t[end]-data.t[1])
hist(data) = plot(data.t,[data.Cy data.Cx Cpow.(data)],label=[L"$C_y$" L"$C_x$" L"$\overline{C}_P$"*" = $(-round(mCpow(data),digits=3))"],xlabel=L"$tU/D$")
hist(data);plot!(xlims=(45,80), ylims=(-10,10), legend=:bottomright, background_color_legend=RGBA{Float64}(1, 1, 1, 0.9),
    grid=true)
savefig("OscCyl_hist_p$p.pdf")

function flow_path_geom(sim,clim=5)
    # Vorticity of z-averaged velocity
    n = size(sim.flow.u,3)-2
    ω = view(sim.flow.σ|>Array,:,:,2)
    u_bar = sum(sim.flow.u,dims=3)./(n+2) |> Array
    u = view(u_bar,:,:,1,1:2)
    @inside ω[I] = WaterLily.curl(3,I,u)*sim.L/sim.U
    contourf(clamp.(ω[inside(ω)]',-clim,clim),aspect_ratio=:equal,
        linewidth=0,linecolor=RGBA{Float64}(1,1,1,0),c=:RdBu_11,clims=(-clim,clim),levels=10,
        border=:none,legend=false,show=false)

    # Add path and geometry
    cen(t) = SA[n+1.5,n+1.5]-sim.body.map(SA[0,0,0],t*sim.L/sim.U)[1:2]
    path = cen.(0:0.1:6)
    plot!(first.(path),last.(path),c=:black,legend=nothing)
    x = 0.5sim.L*sin.(0:π/20:2π) .+ cen(sim_time(sim))[1]
    y = 0.5sim.L*cos.(0:π/20:2π) .+ cen(sim_time(sim))[2]
    plot!(Shape(x,y),fill=(true,:grey))
end
flow_path_geom(sim, 3)
savefig("OscCyl_flow_p$p.pdf")