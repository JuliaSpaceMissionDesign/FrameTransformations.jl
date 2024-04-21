using Test
using FrameTransformations
using Tempo

g = FrameSystem{2, Float64, BarycentricDynamicalTime}()

@test get_order(g) == 2 
@test get_timescale(g) == BarycentricDynamicalTime
@test get_points(g) === g.points
@test get_axes(g) === g.axes

@test !has_point(g, 1)
@test !has_axes(g, 1)
