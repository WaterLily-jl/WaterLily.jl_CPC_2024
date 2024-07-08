using ParametricBodies
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
    for uvᵢ in LinRange(l.lims...,64)
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
