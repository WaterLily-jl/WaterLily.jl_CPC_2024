using GLMakie,Meshing

function geom!(md,d,sim,t=WaterLily.time(sim))
    a = sim.flow.σ
    WaterLily.measure_sdf!(a,sim.body,t)
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
    alg = Meshing.MarchingCubes()
    ranges = range.((0, 0, 0), size(md))
    points, faces = Meshing.isosurface(md, alg, ranges...)
    p3f = Point3f.(points)
    gltriangles = GLMakie.GLTriangleFace.(faces)
    return GLMakie.normal_mesh(p3f, gltriangles)
end

function ω!(md,d,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
end

Makie.inline!(false)
