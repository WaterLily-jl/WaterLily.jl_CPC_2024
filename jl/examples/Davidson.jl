# Davidon minimizer (should be in Optim, but isn't yet)
using ForwardDiff: Dual,Tag,value,partials
using TypedTables
function davidon(f::F,ax::R,bx::R;tol=1e-2,verbose=false,itmx=1000) where {F,R}
    function eval(x)
        T = typeof(Tag(f,R))
        fx = f(Dual{T}(x,one(R)))
        (x=x,f=value(fx),∂=partials(T,fx,1))
    end
    a,b,Δ = eval(ax),eval(bx),bx-ax
    verbose && (points = [a,b])
    for _ in 1:itmx
        v = a.∂+b.∂-3*(b.f-a.f)/Δ; w = copysign(√(v^2-a.∂*b.∂),Δ)
        x = b.x-Δ*(b.∂+w-v)/(b.∂-a.∂+2w)
        c = eval(clamp(x,min(a.x,b.x)+max(Δ/8,tol),max(a.x,b.x)-max(Δ/8,tol)))
        verbose && push!(points,c)
        c.f > max(a.f,b.f) && return ErrorException("Process failed")
        (c.f < min(a.f,b.f) ?  c.∂*Δ < 0 : a.f > b.f) ? (a=c) : (b=c)
        Δ = b.x-a.x; abs(Δ) ≤ 2tol && break
    end
    !verbose && return (a.f < b.f ? a : b)
    points |> Table
end
# test minimizer
@assert davidon(x->(x+3)*(x-1)^2,-2.,2.,tol=1e-6).x ≈ 1
@assert davidon(x->-log(x)/x,1.,10.,tol=1e-6).x ≈ exp(1)
@assert davidon(x->cos(x)+cos(3x)/3,0.,1.75π,tol=1e-7).x ≈ π