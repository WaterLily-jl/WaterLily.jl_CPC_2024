using CUDA,WaterLily,StaticArrays
CUDA.allowscalar(false)

function jelly(p=5;Re=5e2,mem=Array,U=1)
    # Define simulation size, geometry dimensions, & viscosity
    n = 2^p; R = 2n/3; h = 4n-2R; ν = U*R/Re

    # Motion functions
    ω = 2U/R
    @fastmath @inline A(t) = 1 .- SA[1,1,0]*0.1*cos(ω*t)
    @fastmath @inline B(t) = SA[0,0,1]*((cos(ω*t)-1)*R/4-h)
    @fastmath @inline C(t) = SA[0,0,1]*sin(ω*t)*R/4

    # Build jelly from a mapped sphere and plane
    sphere = AutoBody((x,t)->abs(√sum(abs2,x)-R)-1, # sdf
                      (x,t)->A(t).*x+B(t)+C(t))     # map
    plane = AutoBody((x,t)->x[3]-h,(x,t)->x+C(t))
    body =  sphere-plane

    # Return initialized simulation
    Simulation((n,n,4n),(0,0,-U),R;ν,body,mem,T=Float32)
end

include("ThreeD_plots.jl")
function mirrorto!(a,b)
    n = size(b,1)
    a[reverse(1:n),reverse(1:n),:].=b
    a[reverse(n+1:2n),1:n,:].=a[1:n,1:n,:]
    a[:,reverse(n+1:2n),:].=a[:,1:n,:]
    return a
end

begin
    # Define geometry and motion on GPU
    sim = jelly(mem=CUDA.CuArray);
    sim_step!(sim,3);

    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(2,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.1, color=:red)

    #Set up flow viz
    ω = ω!(md,d,sim) |> Observable;
    volume!(ω, algorithm=:mip, colormap=:algae, colorrange=(1,10))
    fig
end

foreach(1:6) do frame
    sim_step!(sim,sim_time(sim)+0.5);
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
    save("Jelly_$frame.png",fig)
end
