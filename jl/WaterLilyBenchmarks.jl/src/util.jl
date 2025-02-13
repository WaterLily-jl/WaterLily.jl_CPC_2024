using KernelAbstractions
using CUDA
using AMDGPU

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
            KernelAbstractions.synchronize($b)
        end
    end |> esc
end

function add_to_suite!(suite, sim_function; p=(3,4,5), s=100, ft=Float32, backend=Array, bstr="CPU", remeasure=false)
    suite[bstr] = BenchmarkGroup([bstr])
    for n in p
        sim = sim_function(n, backend; T=ft)
        sim_step!(sim, typemax(ft); max_steps=1, verbose=false, remeasure=remeasure) # warm up
        suite[bstr][repr(n)] = BenchmarkGroup([repr(n)])
        KA_backend = KernelAbstractions.get_backend(sim.flow.p)
        @add_benchmark sim_step!($sim, $typemax($ft); max_steps=$s, verbose=false, remeasure=$remeasure) $KA_backend suite[bstr][repr(n)] "sim_step!"
    end
end

waterlily_dir = get(ENV, "WATERLILY_DIR", "")
git_hash = read(`git -C $waterlily_dir rev-parse --short HEAD`, String) |> x -> strip(x, '\n')
hostname_dict = Dict{String, String}("alogin" => "MN5", "uan" => "LUMI")
hostname = gethostname()
for (k,v) in hostname_dict
    if occursin(k, hostname)
        hostname = v
    end
end
getf(str) = eval(Symbol(str))

# Fancy logarithmic scale ticks for plotting
# https://github.com/JuliaPlots/Plots.jl/issues/3318
using Plots
"""
    get_tickslogscale(lims; skiplog=false)
Return a tuple (ticks, ticklabels) for the axis limit `lims`
where multiples of 10 are major ticks with label and minor ticks have no label
skiplog argument should be set to true if `lims` is already in log scale.
"""
function get_tickslogscale(lims::Tuple{T, T}; skiplog::Bool=false) where {T<:AbstractFloat}
    mags = if skiplog
        # if the limits are already in log scale
        floor.(lims)
    else
        floor.(log10.(lims))
    end
    rlims = if skiplog; 10 .^(lims) else lims end

    total_tickvalues = []
    total_ticknames = []

    rgs = range(mags..., step=1)
    for (i, m) in enumerate(rgs)
        if m >= 0
            tickvalues = range(Int(10^m), Int(10^(m+1)); step=Int(10^m))
            ticknames  = vcat([string(round(Int, 10^(m)))],
                              ["" for i in 2:9],
                              [string(round(Int, 10^(m+1)))])
        else
            tickvalues = range(10^m, 10^(m+1); step=10^m)
            ticknames  = vcat([string(10^(m))], ["" for i in 2:9], [string(10^(m+1))])
        end

        if i==1
            # lower bound
            indexlb = findlast(x->x<rlims[1], tickvalues)
            if isnothing(indexlb); indexlb=1 end
        else
            indexlb = 1
        end
        if i==length(rgs)
            # higher bound
            indexhb = findfirst(x->x>rlims[2], tickvalues)
            if isnothing(indexhb); indexhb=10 end
        else
            # do not take the last index if not the last magnitude
            indexhb = 9
        end

        total_tickvalues = vcat(total_tickvalues, tickvalues[indexlb:indexhb])
        total_ticknames = vcat(total_ticknames, ticknames[indexlb:indexhb])
    end
    return (total_tickvalues, total_ticknames)
end

"""
    fancylogscale!(p; forcex=false, forcey=false)
Transform the ticks to log scale for the axis with scale=:log10.
forcex and forcey can be set to true to force the transformation
if the variable is already expressed in log10 units.
"""
function fancylogscale!(p::Plots.Subplot; forcex::Bool=false, forcey::Bool=false)
    kwargs = Dict()
    for (ax, force, lims) in zip((:x, :y), (forcex, forcey), (xlims, ylims))
        axis = Symbol("$(ax)axis")
        ticks = Symbol("$(ax)ticks")

        if force || p.attr[axis][:scale] == :log10
            # Get limits of the plot and convert to Float
            ls = float.(lims(p))
            ts = if force
                (vals, labs) = get_tickslogscale(ls; skiplog=true)
                (log10.(vals), labs)
            else
                get_tickslogscale(ls)
            end
            kwargs[ticks] = ts
        end
    end

    if length(kwargs) > 0
        plot!(p; kwargs...)
    end
    p
end
fancylogscale!(p::Plots.Plot; kwargs...) = (fancylogscale!(p.subplots[1]; kwargs...); return p)
fancylogscale!(; kwargs...) = fancylogscale!(plot!(); kwargs...)
