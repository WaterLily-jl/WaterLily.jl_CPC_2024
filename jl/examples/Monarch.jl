# using CUDA
using GLMakie, Meshing, GeometryBasics
using WaterLily, StaticArrays
using ParametricBodies

points = SA{Float32}[0 24 76 132 169 177 176 163 122 88 89  80  64  42  21   4   2  4 0
                     0 -4 -8  -4   9  23  35  46  60 82 88 107 122 130 128 117 103 40 0] |> reverse

planform = ParametricBodies.interpNurbs(points,p=3)

# Check geometry
lines(planform.(range(0,1,100),0))
scatter!(points)

wing2D = ParametricBody(planform,(0,1));
measure(wing2D,SA{Float32}[50,50],0)
measure(wing2D,SA{Float32}[100,100],0)

struct PlanarParametricBody{T,P<:ParametricBody,F<:Function} <: AbstractParametricBody
    body::P
    map::F
    scale::T
end
function surf_props(body::PlanarParametricBody,x::SVector{3},t)
    # Get properties on ζ=0 plane
    ξ = body.map(x,t)
    d,n,_ = ParametricBodies.surf_props(body.body,ξ[1:2],t)

    # Add out of plane contribution
    p = SA[max(d,0)*n[1],max(d,0)*n[2],ξ[3]]
    n =  p/(eps(d)+√(p'*p))
    return body.scale*p'*n,n
end
using ForwardDiff
function WaterLily.measure(body::PlanarParametricBody,x::SVector{3},t)
    # Surf props and velocity in ξ-frame
    d,n = surf_props(body,x,t)
    dξdt = ForwardDiff.derivative(t->-body.map(x,t),t)

    # Convert to x-frame with dξ/dx⁻¹ (d has already been scaled)
    dξdx = ForwardDiff.jacobian(x->body.map(x,t),x)
    return (d,dξdx\n/body.scale,dξdx\dξdt)
end

wing3D = PlanarParametricBody(wing2D,(x,t)->x,1f0)
measure(wing3D,SA{Float32}[50,50,0],0)
measure(wing3D,SA{Float32}[50,50,50],0)
measure(wing3D,SA{Float32}[100,100,0],0)
measure(wing3D,SA{Float32}[100,100,10],0)
