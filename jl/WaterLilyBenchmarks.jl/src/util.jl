using KernelAbstractions: synchronize, get_backend
using CUDA: CuArray
# using GLMakie


function parse_cla(args; cases=["tgv"], log2p=[(6,7)], max_steps=[100], ftype=[Float32], backend=Array)
    iarg(arg) = occursin.(arg, args) |> findfirst
    arg_value(arg) = split(args[iarg(arg)], "=")[end]
    metaparse(x) = eval(Meta.parse(x))

    cases = !isnothing(iarg("cases")) ? arg_value("cases") |> metaparse : cases
    log2p = !isnothing(iarg("log2p")) ? arg_value("log2p") |> metaparse : log2p
    max_steps = !isnothing(iarg("max_steps")) ? arg_value("max_steps") |> metaparse : max_steps
    ftype = !isnothing(iarg("ftype")) ? arg_value("ftype") |> metaparse : ftype
    backend = !isnothing(iarg("backend")) ? arg_value("backend") |> x -> eval(Symbol(x)) : backend
    return cases, log2p, max_steps, ftype, backend
end

macro add_benchmark(args...)
    ex, b, suite, label = args
    return quote
        $suite[$label] = @benchmarkable begin
            $ex
            synchronize($b)
        end
    end |> esc
end

function add_to_suite!(suite, sim_function; p=(3,4,5), s=100, ft=Float32, backend=Array, bstr="CPU")#x$(Threads.nthreads())")
    suite[bstr] = BenchmarkGroup([bstr])
    for n in p
        sim = sim_function(n, backend; T=ft)
        sim_step!(sim, typemax(ft); max_steps=1, verbose=false, remeasure=false) # warm up
        suite[bstr][repr(n)] = BenchmarkGroup([repr(n)])
        KA_backend = get_backend(sim.flow.p)
        @add_benchmark sim_step!($sim, $typemax($ft); max_steps=$s, verbose=false, remeasure=false) $KA_backend suite[bstr][repr(n)] "sim_step!"
    end
end

waterlily_dir = get(ENV, "WATERLILY_DIR", "")
git_hash = read(`git -C $waterlily_dir rev-parse --short HEAD`, String) |> x -> strip(x, '\n')
getf(str) = eval(Symbol(str))

# function flow_ωmag!(dat,sim)
#     a, dt = sim.flow.σ, sim.L/sim.U
#     @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
#     copyto!(dat,a[inside(a)])
# end

# function visualize!(sim; Δt=0.1, nt=100, remeasure=false)
#     dat = sim.flow.σ[inside(sim.flow.σ)] |> Array
#     obs = dat |> Observable

#     f = contour(obs, levels=[-5,5], colormap=:balance)
#     display(f)
#     # plot!(body_mesh(sim))
#     for _ in range(1, nt)
#         sim_step!(sim, sim_time(sim) + Δt; remeasure=remeasure, verbose=true)
#         obs[] = flow_ωmag!(dat, sim)
#     end
# end