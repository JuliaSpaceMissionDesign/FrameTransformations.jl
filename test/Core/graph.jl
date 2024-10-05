using FrameTransformations

using Tempo
using StaticArrays
using Test 

@testset "Point node" begin
    begin
        fps = FrameTransformations.FramePointFunctions{1, Float64}()
        fp = FrameTransformations.FramePointNode{1, Float64}(:P0, 1, 0, 1, fps)

        @test fp.name == :P0
        @test fp.id == 1 
        @test fp.parentid == 0
        @test fp.axesid == 1
        @test fp.f === fps
    end

    begin 
        f1(t) = Translation{2}(SA[cos(t), sin(t), 1.0])
        f2(t) = Translation{2}(SA[cos(t), sin(t), 1.0, -sin(t), cos(t), 0.0])
        fps = FrameTransformations.FramePointFunctions{2, Float64}(f1, f2)

        @test f1(π/3) == fps[1].fw[1](π/3)
        @test f2(π/3) == fps[2].fw[1](π/3)
        @test f1(0)[2] == zeros(3)
    end
end;

@testset "Axes node" begin
    begin 
        faxs = FrameTransformations.FrameAxesFunctions{1, Float64}()
        fax = FrameTransformations.FrameAxesNode{1, Float64}(:Ax0, 1, 1, faxs)

        @test fax.name == :Ax0 
        @test fax.id == 1
        @test fax.parentid == 1
        @test fax.f === faxs
    end
end;

@testset "FrameSystem" begin 
    g = FrameSystem{2, Float64}()

    @test order(g) == 2 
    @test FrameTransformations.timescale(g) == BarycentricDynamicalTime
    @test points_graph(g) === g.points.graph
    @test axes_graph(g) === g.axes.graph
    @test points_alias(g) === g.points.alias
    @test axes_alias(g) === g.axes.alias

    @test !has_point(g, 1)
    @test !has_axes(g, 1)
end;