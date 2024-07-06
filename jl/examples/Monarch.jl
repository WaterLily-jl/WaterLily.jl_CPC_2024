using CUDA
using GLMakie, Meshing, GeometryBasics
using WaterLily, StaticArrays
using ParametricBodies
Makie.inline!(false)
CUDA.allowscalar(false)

# Define simple NewtonLocator
struct NewtonLocator{T,F<:Function,G<:Function} <: AbstractLocator
    refine::F
    surf::G
    lims::NTuple{2,T} 
    scale::T
end
function NewtonLocator(curve::Function,lims::NTuple{2},scale=1)
    T = promote_type(eltype(lims),typeof(scale),Float32) # Need a floating point
    f = ParametricBodies.refine(curve,T.(lims),curve(first(lims),0)≈curve(last(lims),0))
    NewtonLocator(f,curve,T.(lims),T(scale))
end
ParametricBodies.notC¹(l::NewtonLocator,uv) = any(uv.≈l.lims)
function (l::NewtonLocator{T})(x,t) where T
    # Grid search for uv within bounds
    @inline dis2(uv) = (q=x-l.surf(uv,t); q'*q)
    uv = first(l.lims); d = dis2(uv)
    for uvᵢ in LinRange(l.lims...,128)
        dᵢ = dis2(uvᵢ)
        dᵢ<d && (uv=uvᵢ; d=dᵢ)
    end
    d*l.scale^2>100 && return uv

    # If close, refine estimate with two Newton steps
    uv = l.refine(x,uv,t)
    return l.refine(x,uv,t)
end

# Define struct for a planar ParametricBody
struct PlanarParametricBody{T,P<:ParametricBody,F<:Function} <: AbstractParametricBody
    body::P
    map::F
    scale::T
end
function PlanarParametricBody(curve,lims::Tuple;T=Float32,map=(x,t)->x)
    # Wrap in type safe functions (GPUs are picky)
    wcurve(u::U,t::T) where {U,T} = SVector{2,promote_type(U,T)}(curve(u,t))
    wmap(x::SVector{n,X},t::T) where {n,X,T} = SVector{n,promote_type(X,T)}(map(x,t))

    scale = T(ParametricBodies.get_scale(map,SA{T}[0,0,0]))
    locate = NewtonLocator(wcurve,T.(lims),scale)
    body = ParametricBody(wcurve,locate)
    PlanarParametricBody(body,wmap,scale)
end

function ParametricBodies.surf_props(body::PlanarParametricBody,x::SVector{3},t;ϵ=1)
    # Get properties on ζ=0 plane
    ξ = body.map(x,t)
    d,n,_ = ParametricBodies.surf_props(body.body,SA[ξ[1],ξ[2]],t)

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
WaterLily.loc(i,I::CartesianIndex{N},T=Float32) where N = SVector{N,T}(I.I .- 1.5 .- 0.5 .* δ(i,I).I)

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
    body=PlanarParametricBody(planform,(0,1);map)
    Simulation((3L,2L,4L),(0.2,0,-0.2),L;U,ν=U*L/Re,body,T,mem)
end

begin
    # Define geometry and motion on CPU
    sim = monarch(mem=CuArray);
    sim_step!(sim,3);

    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(1,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.3, color=:cyan)

    # #Set up flow viz
    ω = ω!(md,d,sim) |> Observable;
    volume!(ω, algorithm=:mip, colormap=:dense, colorrange=(1,30))
    fig
end

foreach(1:6) do frame
    sim_step!(sim,sim_time(sim)+1.25);
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
    save("Butterfly_$frame.png",fig)
end