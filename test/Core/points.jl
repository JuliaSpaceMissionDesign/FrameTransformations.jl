using Test 
using StaticArrays
using ReferenceFrameRotations
using FrameTransformations

@testset "Root" verbose = false begin
    frames = FrameSystem{4, Float64}()

    # test unknown axes 
    @test_throws ArgumentError add_point_root!(frames, :Root, 1, 1)

    # test unknown axes 
    add_axes_root!(frames, :Ax, 1)
    add_point_root!(frames, :Root, 1, 1)

    # test point properties
    node = get_points(frames).nodes[1]
    @test node.class == 0
    @test node.name == :Root
    @test node.axesid == 1
    @test node.id == 1

    # test root point already exists
    @test_throws ArgumentError add_point_root!(frames, :Root, 1, 1)
    @test_throws ArgumentError add_point_root!(frames, :Root, 2, 1)
end

@testset "Fixed-offset" begin
    frames = FrameSystem{4, Float64}()

    add_axes_root!(frames, :Ax, 1)
    add_point_root!(frames, :Root, 1, 1)

    # test wrong offset dimension (should be 3)
    offset = rand(5)
    @test_throws DimensionMismatch add_point_fixedoffset!(frames, :Point, 2, 1, 1, offset)

    # test unknown parent point 
    offset = rand(3)
    @test_throws ArgumentError add_point_fixedoffset!(frames, :Point, 2, 0, 1, offset)
    @test_throws ArgumentError add_point_fixedoffset!(frames, :Point, 2, 1, 0, offset)

    add_point_fixedoffset!(frames, :Fixed, 2, 1, 1, offset)

    # test point properties 
    node = get_points(frames).nodes[2]
    @test node.class == 1
    @test node.name == :Fixed
    @test node.axesid == 1
    @test node.id == 2

    # test point transformation 
    x = vector3(frames, 2, 1, 1, rand())
    @test x ≈ -offset atol = 1e-12 rtol = 1e-12

    x = vector6(frames, 2, 1, 1, rand())
    @test x[1:3] ≈ -offset atol = 1e-12 rtol = 1e-12
    @test x[4:6] == zeros(3)

    x = vector9(frames, 2, 1, 1, rand())
    @test x[1:3] ≈ -offset atol = 1e-12 rtol = 1e-12
    @test x[4:9] == zeros(6)

    x = vector12(frames, 2, 1, 1, rand())
    @test x[1:3] ≈ -offset atol = 1e-12 rtol = 1e-12
    @test x[4:12] == zeros(9)
end

@testset "Dynamical" begin

    # Numerical functions
    f = t -> SA[sin(t) * ones(3)...]
    df = t -> SA[sin(t) * ones(3)..., cos(t) * ones(3)...]
    ddf = t -> SA[sin(t) * ones(3)..., cos(t) * ones(3)..., -sin(t) * ones(3)...]
    dddf = t -> SA[sin(t) * ones(3)..., cos(t) * ones(3)..., -sin(t) * ones(3)..., -cos(t) * ones(3)...]

    frames = FrameSystem{4, Float64}()
    add_axes_root!(frames, :Ax, 1)
    add_point_root!(frames, :Root, 1, 1)
    add_point_dynamical!(frames, :Dyn, 2, 1, 1, f)

    # test point properties 
    node = get_points(frames).nodes[2]
    @test node.class == 2
    @test node.name == :Dyn
    @test node.axesid == 1
    @test node.id == 2

    R = angle_to_dcm(π/4, :Z)

    atol, rtol = 1e-12, 1e-12
    # test AD derivatives for all combinations of specified functions
    for funs in ((f,), (f, df), (f, df, ddf), (f, df, ddf, dddf))
        frames = FrameSystem{4,Float64}()
        add_axes_root!(frames, :Ax, 1)
        add_point_root!(frames, :Root, 1, 1)
        add_point_dynamical!(frames, :Dyn, 2, 1, 1, funs...)
        add_axes_fixedoffset!(frames, :Fox, 2, 1, R)

        node = get_points(frames).nodes[2]
        
        for _ in 1:10
            t = rand()
            
            # test wrappers
            y = node.f.fun[1](t)
            @test y[1:3] ≈ f(t) atol = atol rtol = rtol

            y = node.f.fun[2](t)
            @test y[1:6] ≈ df(t) atol = atol rtol = rtol

            y = node.f.fun[3](t)
            @test y[1:9] ≈ ddf(t) atol = atol rtol = rtol

            y = node.f.fun[4](t)
            @test y[1:12] ≈ dddf(t) atol = atol rtol = rtol

            # test high-level
            x = vector3(frames, 2, 1, 2, t)
            @test x[1:3] ≈ -R * sin(t) * ones(3) atol = atol rtol = rtol

            x = vector6(frames, 2, 1, 2, t)
            @test x[1:3] ≈ -R * sin(t) * ones(3) atol = atol rtol = rtol
            @test x[4:6] ≈ -R * cos(t) * ones(3) atol = atol rtol = rtol

            x = vector9(frames, 2, 1, 2, t)
            @test x[1:3] ≈ -R * sin(t) * ones(3) atol = atol rtol = rtol
            @test x[4:6] ≈ -R * cos(t) * ones(3) atol = atol rtol = rtol
            @test x[7:9] ≈ R * sin(t) * ones(3) atol = atol rtol = rtol

            x = vector12(frames, 2, 1, 2, t)
            @test x[1:3] ≈ -R * sin(t) * ones(3) atol = atol rtol = rtol
            @test x[4:6] ≈ -R * cos(t) * ones(3) atol = atol rtol = rtol
            @test x[7:9] ≈ R * sin(t) * ones(3) atol = atol rtol = rtol
            @test x[10:12] ≈ R * cos(t) * ones(3) atol = atol rtol = rtol
        end
    end

end