
function makie_video!(makie_plot,sim,dat,obs_update!;remeasure=false,name="file.mp4",duration=1,step=0.1,framerate=30,compression=20)
    # Set up viz data and figure
    obs = obs_update!(dat,sim) |> Observable;
    f = makie_plot(obs)

    # Run simulation and update figure data
    t₀ = round(sim_time(sim))
    t = range(t₀,t₀+duration;step)

    GLMakie.record(f, name, t; framerate, compression) do tᵢ
        sim_step!(sim,tᵢ;remeasure)
        obs[] = obs_update!(dat,sim)
        println("simulation ",round(Int,(tᵢ-t₀)/duration*100),"% complete")
    end

    # GLMakie.record(f, name, t; framerate, compression) do tᵢ
    #     sim_step!(sim,tᵢ;remeasure)
    #     obs[] = obs_update!(dat,sim)
    #     println("simulation ",round(Int,(tᵢ-t₀)/duration*100),"% complete")
    # end
    return f
end

using Meshing, GeometryBasics
function body_mesh(sim,t=0)
    a = sim.flow.σ; R = inside(a)
    WaterLily.measure_sdf!(a,sim.body,t)
    normal_mesh(GeometryBasics.Mesh(a[R]|>Array,MarchingCubes(),origin=Vec(0,0,0),widths=size(R)))
end;
function flow_λ₂!(dat,sim)
    a = sim.flow.σ
    @inside a[I] = max(0,log10(-min(-1e-6,WaterLily.λ₂(I,sim.flow.u)*(sim.L/sim.U)^2))+.25)
    copyto!(dat,a[inside(a)])                  # copy to CPU
end
function flow_λ₂(sim)
    dat = sim.flow.σ[inside(sim.flow.σ)] |> Array
    flow_λ₂!(dat,sim)
    dat
end

function flow_ω_mag!(dat,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(dat,a[inside(a)])                  # copy to CPU
end

using Meshing, GeometryBasics
function geom!(md,d,sim,t=WaterLily.time(sim))
    a = sim.flow.σ
    WaterLily.measure_sdf!(a,sim.body,t)
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
    normal_mesh(GeometryBasics.Mesh(md,Meshing.MarchingCubes(),origin=Vec(0,0,0),widths=size(md)))
end

function ω!(md,d,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
end

function ω_θ!(dat,sim,center=center)
    dt, a = sim.L/sim.U, sim.flow.σ
    @inside a[I] = WaterLily.ω_θ(I,(1,0,0),center,sim.flow.u)*dt
    copyto!(dat,a[inside(a)])
end