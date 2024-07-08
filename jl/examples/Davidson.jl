# Davidson minimizer (should be in Optim, but isn't yet)
using ForwardDiff: Dual,Tag,value,partials
using TypedTables
function davidson(f::F,ax::R,bx::R;tol=1e-2,verbose=true) where {F,R}
    function eval(f,x)
        T = typeof(Tag(f,R))
        fx = f(Dual{T}(x,one(R)))
        (x=x,f=value(fx),∂=partials(T,fx,1))
    end
    a,b,Δ = eval(f,ax),eval(f,bx),bx-ax
    points = [a,b]
    verbose && println("1: ",a,"\n2: ",b)
    while abs(Δ) > 2tol
        v = a.∂+b.∂-3*(b.f-a.f)/Δ; w = copysign(√(v^2-a.∂*b.∂),Δ)
        x = b.x-Δ*(b.∂+w-v)/(b.∂-a.∂+2w)
        (x+tol>max(a.x,b.x) || x-tol<min(a.x,b.x)) && (x = 0.5*(a.x+b.x))
        c = eval(f,x); push!(points,c)
        verbose && println(length(points),": ",c)
        c.f > max(a.f,b.f) && return ErrorException("Process failed")
        c.∂*a.∂ > 0 ? (a=c) : (b=c)
        Δ = b.x-a.x
    end
    (a.f < b.f ? a : b), points |> Table
end
opt,_=davidson(x->-log(x)/x,1.,10.,tol=1e-6); @assert opt.x ≈ exp(1) # test minimizer