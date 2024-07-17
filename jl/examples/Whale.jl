using CUDA, WaterLily, StaticArrays
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
function whale(sweep=100;L=24,U=1,Re=1e6,T=Float32,mem=Array)
    pnts = SA{Float32}[0 40      190      200      190   170      100   0 -10 0
            0  0 0.8sweep 0.8sweep+40 0.5sweep+70 0.5sweep+50 0.5sweep+70 100  80 0]
    planform = BSplineCurve(reverse(pnts),degree=3)
    function map(x,t)
        θ,h = π/6*sin(U*t/L),1+0.5cos(U*t/L)
        Ry = SA[cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
        Ry*100*(x/L-SA[0.75,0,h])
    end
    body=PlanarParametricBody(planform,(0,1);map,mem)
    Simulation((5L,3L,2L),(1.,0.,0.),L;U,ν=U*L/Re,body,T,mem)
end

begin
    # Define geometry and motion on GPU
    sim = whale(0,mem=CuArray);  # closer to real-time
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
    volume!(ω, algorithm=:mip, colormap=:dense, colorrange=(1,20))
    fig
end

# simulate real-time
for frame ∈ 1:100
    println(frame)
    sim_step!(sim,sim_time(sim)+0.05);
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
end