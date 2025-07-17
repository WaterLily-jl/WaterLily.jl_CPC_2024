using Printf, StaticArrays, CUDA, JLD2, Plots, LaTeXStrings, WaterLily

# Utils
iarg(arg) = occursin.(arg, ARGS) |> findfirst
iarg(arg, args) = occursin.(arg, args) |> findfirst
arg_value(arg) = split(ARGS[iarg(arg)], "=")[end]
arg_value(arg, args) = split(args[iarg(arg, args)], "=")[end]
metaparse(x) = eval(Meta.parse(x))
Plots.default(
    fontfamily = "Computer Modern",
    linewidth = 1,
    framestyle = :box,
    grid = false,
    left_margin = Plots.Measures.Length(:mm, 5),
    right_margin = Plots.Measures.Length(:mm, 5),
    bottom_margin = Plots.Measures.Length(:mm, 5),
    top_margin = Plots.Measures.Length(:mm, 5),
    titlefontsize = 15,
    legendfontsize = 14,
    tickfontsize = 14,
    labelfontsize = 14,
)

# Params
T = Float32
mem = !isnothing(iarg("mem", ARGS)) ? arg_value("mem", ARGS) |> x -> eval(Symbol(x)) : CuArray
U = 1
Re = 3700
Ds = [88, 128, 168]
D = !isnothing(iarg("D", ARGS)) ? arg_value("D", ARGS) |> metaparse : 88 # diameter resolution [88,128,168]
Lc = !isnothing(iarg("Lc", ARGS)) ? arg_value("Lc", ARGS) |> metaparse : 3 # domain side length [3,6]
L = (7,Lc,Lc) # domain size in D: (7,3,3), (7,6,6)
center = (1.5,L[2]/2,L[3]/2)
probe_loc = (4.5,center[2]+1/2,center[3]+1/2)
half_domain = D * L[2] ÷ 2

time_max = 400 # in CTU
stats_init = 100 # in CTU
stats_interval = 0.1 # in CTU
save_interval = 50 # in CTU

data_dir = !isnothing(iarg("datadir", ARGS)) ? arg_value("datadir", ARGS) |> String : "data/sphere/"
load_time = !isnothing(iarg("loadtime", ARGS)) ? arg_value("loadtime", ARGS) |> metaparse : nothing
fname_save = joinpath(data_dir,"D$(D)_Lc$(Lc)")
pdf_file = "../../tex/img/sphere_validation.pdf"
restart = nothing
restart_stats = nothing
restart_signals = nothing

run_flag = !isnothing(iarg("run", ARGS)) ? arg_value("run", ARGS) |> metaparse : true
plot_flag = !isnothing(iarg("plot", ARGS)) ? arg_value("plot", ARGS) |> metaparse : false

# I/O
function save_sim(sim, meanflow, force, probe, t)
    t_str = @sprintf("%i", sim_time(sim))
    save!(fname_save * "_t$(t_str).jld2", sim)
    save!(fname_save * "_t$(t_str)_meanflow.jld2", meanflow)
    jldsave(fname_save * "_t$(t_str)_signals.jld2"; force, probe, t)
end
function read_forces(fname)
    obj = jldopen(fname)
    return obj["force"], obj["probe"], obj["t"]
end
# Metrics
recirculation_length(U; Cu=1.5) = findfirst(axes(U,1)) do i
    (U[i,half_domain,half_domain,1] < 0) && (U[i+1,half_domain,half_domain,1] >= 0)
end |> i -> i/D-Cu

# Sphere
function sphere(D; Re=3700, U=1, L=(7,3,3), center=(1.5,1.5,1.5), mem=CuArray, T=Float32)
    center = SA{T}[center...]
    body = AutoBody((x,t) -> √sum(abs2, x .- (center .* D)) - D÷2)
    Simulation(L.*D, (U, 0, 0), D; U, ν=U*D/Re, body, mem, T, exitBC=true)
end

function run_sim(D; Re=3700, U=1, L=(7,3,3), center=(1.5,1.5,1.5), mem=CuArray, T=Float32,
    restart=nothing, restart_stats=nothing, restart_signals=nothing)

    sim = sphere(D; Re, U, L, center, mem, T)

    force, probe, t = Vector{T}[], T[], T[] # force coefficients, velocity probe, time
    probe_loc_i = @. (probe_loc * sim.L) |> floor |> Int

    if !isnothing(restart)
        load!(sim; fname=restart)
        println("Loaded: $restart")
    end
    sim_step!(sim, stats_init; remeasure=false, verbose=true)

    meanflow = MeanFlow(sim.flow; uu_stats=true)
    if !isnothing(restart_stats)
        load!(meanflow; fname=restart_stats)
        println("Loaded: $restart_stats")
        force, probe, t = read_forces(restart_signals)
        println("Loaded: $restart_signals")
    end

    println("\nComputing mean flow statistics from: T=$(sim_time(sim))\n\
    Total accumulated: T=$(WaterLily.time(meanflow)*sim.U/sim.L)\n\
    Remaining: T=$(time_max-stats_init-WaterLily.time(meanflow)*sim.U/sim.L)\n")

    next_save = sim_time(sim) + save_interval
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=false)
        sim_info(sim)

        WaterLily.update!(meanflow, sim.flow)
        push!(force, WaterLily.total_force(sim))
        push!(probe, view(sim.flow.u,probe_loc_i...,2) |> Array |> x->x[])
        push!(t, sim_time(sim))

        if WaterLily.sim_time(sim) > next_save || sim_time(sim) > time_max
            save_sim(sim, meanflow, force, probe, t)
            next_save = sim_time(sim) + save_interval
            println("Saved simulation and mean flow statistics.")
        end
    end
    return sim, meanflow, mapreduce(permutedims, vcat, force), probe, t
end

function run_sphere(D; Re=3700, U=1, L=(7,3,3), center=(1.5,1.5,1.5), mem=CuArray, T=Float32,
    restart=nothing, restart_stats=nothing, restart_signals=nothing, load_time=nothing)
    if isnothing(load_time)
        mkpath(data_dir)
        return run_sim(D; Re, U, L, center, mem, T, restart, restart_stats, restart_signals)
    else
        t_str = @sprintf("%i", load_time)
        sim = sphere(D; Re, U, L, center, mem, T)
        fname = "$(fname_save)_t$(t_str).jld2"
        load!(sim; fname)
        println("Loaded: $fname")

        # meanflow = MeanFlow(sim.flow; uu_stats=true)
        meanflow = MeanFlow(L.*D; uu_stats=false)
        fname = "$(fname_save)_t$(t_str)_meanflow.jld2"
        load!(meanflow; fname)
        println("Loaded: $fname")

        fname = "$(fname_save)_t$(t_str)_signals.jld2"
        force, probe, t = read_forces(fname)
        println("Loaded: $fname")

        return sim, meanflow, mapreduce(permutedims, vcat, force), probe, t
    end
end

function main()
    if run_flag
        println("Running: D=$(D), L=$(L), datadir=$data_dir")
        sim, meanflow, force, _, t = run_sphere(D; Re, U, L, center, mem, T, restart, restart_stats, restart_signals, load_time)
        println("▷ ΔT [CTU] = "*@sprintf("%.2f", WaterLily.time(meanflow)/sim.L*sim.U))
        println("▷ L/D = "*@sprintf("%.2f", recirculation_length(meanflow.U|>Array; Cu=center[1])))
        println("▷ CD_mean = "*@sprintf("%.2f", -sum(force[2:end,1].*diff(t))/sum(diff(t))/(0.5*U^2*π*(D/2)^2)))
    end
    if plot_flag
        p_cd = plot()
        for d in Ds
            force, _, t = read_forces(joinpath(data_dir,"D$(d)_Lc3") * "_t400_signals.jld2")
            CD_mean = -sum(force[2:end,1].*diff(t))/sum(diff(t))/(0.5*U^2*π*(d/2)^2)
            scatter!(p_cd, [d], [CD_mean], label=d==168 ? label=L"$\mathrm{Present}, 7D\times3D\times3D$" : "", color=:grey, ms=6, msc=:grey)
        end
        force, _, t = read_forces(joinpath(data_dir,"D88_Lc6") * "_t400_signals.jld2")
        CD_mean = -sum(force[2:end,1].*diff(t))/sum(diff(t))/(0.5*U^2*π*(88/2)^2)
        scatter!(p_cd, [88], [CD_mean], label=L"$\mathrm{Present}, 7D\times6D\times6D$", color=:black, ms=6, ma=1)
        scatter!(p_cd, grid=true, ms=6, ma=1, ylims=(0.28,0.5001), xlims=(80,176), xticks=Ds, xlabel=L"$D$",
            lw=0, framestyle=:box, size=(600, 600), legend=:bottomright, color=:black, ylabel=L"$\overline{C_D}$",
            legendfontsize=14, tickfontsize=18, labelfontsize=18, left_margin=Plots.Measures.Length(:mm, 5),
        )

        hline!(p_cd, [0.394], linestyle=:dash, color=:blue, label=L"\mathrm{Rodriguez}\,\,et\,\,al\mathrm{.\,\,(DNS)}")
        hline!(p_cd, [0.355], linestyle=:dashdot, color=:green, label=L"\mathrm{Yun}\,\,et\,\,al\mathrm{.\,\,(LES)}")
        fig_path = joinpath(string(@__DIR__), pdf_file)
        println("Figure stored in $(fig_path)")
        savefig(fig_path)
    end
end

main()
