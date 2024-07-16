using CUDA,WaterLily, StaticArrays
CUDA.allowscalar(false)
include("PlanarBodies.jl")
include("ThreeD_plots.jl")
function mirrorto!(a,b)
    n = size(b,2)
    a[:,n+1:2n,:].=b
    a[:,reverse(1:n),:].=b
    return a
end

# Define simulation
function monarch(;L=32,U=1,Re=500,T=Float32,mem=Array)
    points = SA{Float32}[0 24 76 132 169 177 176 163 122 88 89  80  64  42  21   4   2  4 0
                         0 -4 -8  -4   9  23  35  46  60 82 88 107 122 130 128 117 103 40 0] |> reverse
    planform = ParametricBodies.interpNurbs(points,p=3)
    function map(x,t)
        θ = π/8-3π/8*sin(U*t/L)
        Rx = SA[1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
        Rx*100*(x/L-SA[0.75,0.1,1.5])
    end
    body=PlanarParametricBody(planform,(0,1);map,mem)
    Simulation((3L,2L,4L),(0.2,0.,-0.2),L;U,ν=U*L/Re,body,T,mem)
end

begin
    # Define geometry and motion on CPU
    # sim = monarch(mem=CuArray); # figures for paper
    # sim_step!(sim,4);
    sim = monarch(L=24,Re=250,mem=CuArray);  # closer to real-time
    sim_step!(sim,0.1);

    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(1,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.3, color=:fuchsia)

    # #Set up flow viz
    ω = ω!(md,d,sim) |> Observable;
    volume!(ω, algorithm=:mip, colormap=:dense, colorrange=(1,30))
    fig
end

# foreach(1:6) do frame
#     sim_step!(sim,sim_time(sim)+1);
foreach(1:100) do frame
    println(frame)
    sim_step!(sim,sim_time(sim)+0.05);
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
    # save("Butterfly_$frame.png",fig)
end