using Test
using FrameTransformations
using Tempo

g = FrameSystem{2, Float64, BarycentricDynamicalTime}()

@test order(g) == 2 
@test FrameTransformations.timescale(g) == BarycentricDynamicalTime
@test points_graph(g) === g.points
@test axes_graph(g) === g.axes

@test !has_point(g, 1)
@test !has_axes(g, 1)
