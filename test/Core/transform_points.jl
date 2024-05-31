using Test
using FrameTransformations
using Tempo
using StaticArrays

g = FrameSystem{4, Float64, BarycentricDynamicalTime}()

@test_nowarn let 
    fps = FrameTransformations.FramePointFunctions{4, Float64}()
    fp = FrameTransformations.FramePointNode{4, Float64}(:test, 0, 1, 1, 0, fps)
    @test_throws MethodError FrameTransformations.FramePointNode{1, Float64}(:test, 0, 1, 1, 0, fps) 
end

@test_nowarn let 
    fps1 = FrameTransformations.FramePointFunctions{1, Float64}()
    fp1 = FrameTransformations.FramePointNode{1, Float64}(:test, 0, 1, 1, 0, fps1)
    @test_throws MethodError add_point!(g, fp1)
end

fps = FrameTransformations.FramePointFunctions{4, Float64}()
fp = FrameTransformations.FramePointNode{4, Float64}(:P0, 0, 1, 1, 1, fps)

@test fp.name == :P0
@test fp.class == 0 
@test fp.id == 1
@test fp.parentid == 1
@test fp.axesid == 1
@test fp.f === fps

f(t) = SA[cos(2π*t), sin(2π*t), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

fps1 = FrameTransformations.FramePointFunctions{4, Float64}(f)
fp1 = FrameTransformations.FramePointNode{4, Float64}(:P1, 0, 2, 1, 1, fps1)

fps2 = FrameTransformations.FramePointFunctions{4, Float64}(t->-f(t))
fp2 = FrameTransformations.FramePointNode{4, Float64}(:P2, 0, 3, 2, 1, fps2)

for i in 1:2
    @test fp1.f[i](0.1) == -fp2.f[i](0.1)
end

for p in [fp, fp1, fp2]
    add_point!(g, p) 
    FrameTransformations.add_edge!(g.points, p.parentid, p.id)
end

# ---
# Root axes

faxs0 = FrameTransformations.FrameAxesFunctions{4, Float64}()
fax0 = FrameTransformations.FrameAxesNode{4, Float64}(:AX0, 0, 1, 1, faxs0)

add_axes!(g, fax0)
FrameTransformations.add_edge!(g.axes, fax0.parentid, fax0.id)

@test vector3(g, 1, 3, 1, 1.0) ≈ @SVector zeros(Float64, 3)
@test vector6(g, 1, 3, 1, 1.0) ≈ @SVector zeros(Float64, 6)
@test vector9(g, 1, 3, 1, 1.0) ≈ @SVector zeros(Float64, 9)
@test vector12(g, 1, 3, 1, 1.0) ≈ @SVector zeros(Float64, 12)

# Dummy constructor 
funs = FrameTransformations.FramePointFunctions{3, Int}()
@test typeof(funs) == FrameTransformations.FramePointFunctions{3, Int, 9}
@test length(funs.fun) == 3

for i in 1:3
    @test typeof(funs.fun[i](0)) == SVector{9, Int}
end