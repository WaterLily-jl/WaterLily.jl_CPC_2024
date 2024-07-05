# using CUDA
using GLMakie, Meshing, GeometryBasics
using WaterLily, StaticArrays
using ParametricBodies

# Define struct for a planar ParametricBody
struct PlanarParametricBody{T,P<:ParametricBody,F<:Function} <: AbstractParametricBody
    body::P
    map::F
    scale::T
end
function PlanarParametricBody(surf,uv_bounds::Tuple,map,T=Float32)
    body = ParametricBody(surf,uv_bounds)
    scale = ParametricBodies.get_scale(map,SA{T}[0,0,0])
    PlanarParametricBody(body,map,T(scale))
end
function ParametricBodies.surf_props(body::PlanarParametricBody,x::SVector{3},t;ϵ=1)
    # Get properties on ζ=0 plane
    ξ = body.map(x,t)
    d,n,_ = ParametricBodies.surf_props(body.body,ξ[1:2],t)

    # Add out of plane contribution & thickness
    p = SA[max(d,0)*n[1],max(d,0)*n[2],ξ[3]]
    n =  p/(eps(d)+√(p'*p))
    thk = typeof(d)(√3/2+ϵ)
    return body.scale*p'*n-thk,n
end
using ForwardDiff
function WaterLily.measure(body::PlanarParametricBody,x::SVector{3},t)
    # Surf props and velocity in ξ-frame
    d,n = ParametricBodies.surf_props(body,x,t)
    dξdt = ForwardDiff.derivative(t->-body.map(x,t),t)

    # Convert to x-frame with dξ/dx⁻¹ (d has already been scaled)
    dξdx = ForwardDiff.jacobian(x->body.map(x,t),x)
    return (d,dξdx\n/body.scale,dξdx\dξdt)
end

# Plotting functions
function geom!(md,d,sim,t=WaterLily.time(sim))
    a = sim.flow.σ
    WaterLily.measure_sdf!(a,sim.body,t)
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror
    normal_mesh(GeometryBasics.Mesh(md,Meshing.MarchingCubes(),origin=Vec(0,0,0),widths=size(md)))
end

function ω!(md,d,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror
end

function mirrorto!(a,b)
    n = size(b,2)
    a[:,n+1:2n,:].=b
    a[:,reverse(1:n),:].=b
    return a
end

# Define simulation
function monarch(L=32,U=1,Re=500,T=Float32,mem=Array)
    points = SA{T}[0 24 76 132 169 177 176 163 122 88 89  80  64  42  21   4   2  4 0
                         0 -4 -8  -4   9  23  35  46  60 82 88 107 122 130 128 117 103 40 0] |> reverse
    planform = ParametricBodies.interpNurbs(points,p=3)
    function map(x::SVector{3,T},t) where T
        θ = -π/4
        Rx = SA{T}[1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
        Rx*100*(x/T(L)-SA{T}[0.25,0.07,1.5])
    end
    body=PlanarParametricBody(planform,(0,1),map,T)

    Simulation((2L,2L,4L),(0,0,0),L;U,ν=U*L/Re,body,T,mem)
end

Makie.inline!(false)
CUDA.allowscalar(false)
begin
    # Define geometry and motion on CPU
    sim = monarch();
    # sim_step!(sim,3);

    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(1,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.1, color=:orange)

    # #Set up flow viz
    # ω = ω!(md,d,sim) |> Observable;
    # volume!(ω, algorithm=:mip, colormap=:algae, colorrange=(1,10))
    fig
end
