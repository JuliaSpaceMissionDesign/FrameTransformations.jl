using Test 
using StaticArrays
using ReferenceFrameRotations
using FrameTransformations

@testset "Core" verbose=true begin
    frames = FrameSystem{4, Float64}()
    add_axes_icrf!(frames)
    vec = SA[1., 0., 0.]
    
    fun = t->vec
    dfun = t->zeros(3)

    add_direction!(frames, :Dir1, 1, fun)
    add_direction!(frames, :Dir2, 1, fun, fun)
    add_direction!(frames, :Dir3, 1, fun, fun, fun)
    add_direction!(frames, :Dir4, 1, fun, fun, fun, fun)

    @test_throws ArgumentError add_direction!(frames, :Dir1, 1, fun)

    @test direction12(frames, :Dir4, 1, 1.0) == vcat([1., 0., 0.], zeros(9))

    add_direction_fixed!(frames, :DirFixed, 1, vec)
    @test direction12(frames, :DirFixed, 1, 1.0) == vcat([1., 0., 0.], zeros(9))
end;

