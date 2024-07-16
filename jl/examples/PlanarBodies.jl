using ParametricBodies

# Define struct for a planar ParametricBody
struct PlanarParametricBody{T,P<:ParametricBody,F<:Function} <: AbstractParametricBody
    body::P
    map::F
    scale::T
end
function PlanarParametricBody(curve,lims::Tuple;T=Float32,map=(x,t)->x,mem=Array)
    # Wrap in type safe functions (GPUs are picky)
    wcurve(u::U,t::T) where {U,T} = SVector{2,promote_type(U,T)}(curve(u,t))
    wmap(x::SVector{n,X},t::T) where {n,X,T} = SVector{n,promote_type(X,T)}(map(x,t))

    scale = T(ParametricBodies.get_scale(map,SA{T}[0,0,0]))
    locate = HashedLocator(wcurve,T.(lims);step=inv(scale),T,mem)
    body = ParametricBody(wcurve,locate)
    PlanarParametricBody(body,wmap,scale)
end

using Adapt
Adapt.adapt_structure(to, x::PlanarParametricBody) = PlanarParametricBody(adapt(to,x.body),x.map,x.scale)

function ParametricBodies.surf_props(body::PlanarParametricBody,x::SVector{3},t;ϵ=1)
    # Get point and thickness offset
    ξ = body.map(x,t); thk = eltype(ξ)(√3/2+ϵ)

    # Get vector to point
    if body.scale*abs(ξ[3])<2thk # might be close to planar body
        d,n,_ = ParametricBodies.surf_props(body.body,SA[ξ[1],ξ[2]],t)
        p = SA[max(d,0)*n[1],max(d,0)*n[2],ξ[3]]
    else
        p = SA[0,0,ξ[3]] # simple planar approximation
    end

    # return scaled distance and normal
    n = p/(eps(eltype(p))+√(p'*p))
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

WaterLily.loc(i,I::CartesianIndex{N},T=Float32) where N = SVector{N,T}(I.I .- 1.5 .- 0.5 .* δ(i,I).I)
