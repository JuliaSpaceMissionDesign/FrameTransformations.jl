
@point SSB 0 SolarSystemBarycenter
@point Sun 10
@point Earth 399

@testset "Frame types" verbose = true begin
    atol = 1e-10

    @testset "AxesNodes" begin

        # ComputableAxesVector 
        # ------------------------ 
        x = Frames.ComputableAxesVector(Earth, Sun, 1)
        @test x.from == 399
        @test x.to == 10
        @test x.order == 1

        x = Frames.ComputableAxesVector(20, 10, 3)
        @test x.from == 20
        @test x.to == 10
        @test x.order == 3

        @test_throws ArgumentError Frames.ComputableAxesVector(1, 2, 4)
        @test_throws ArgumentError Frames.ComputableAxesVector(1, 2, 0)
        @test_throws ArgumentError Frames.ComputableAxesVector(1, 1, 2)

        # ComputableAxesProperties
        # ------------------------
        x = Frames.ComputableAxesVector(10, 20, 3)
        y = Frames.ComputableAxesVector(2, 4, 2)

        cap = Frames.ComputableAxesProperties(x, y)

        @test cap.v1 == x
        @test cap.v2 == y

        # FrameAxesFunctions 
        # ------------------

        fcn_rot1(t, x, y) = Rotation{1}(angle_to_dcm(t, :Z))
        fcn_rot2(t, x, y) = Rotation{2}(angle_to_dcm(t, :Y))
        fcn_rot3(t, x, y) = Rotation{3}(angle_to_dcm(t, :X))

        s1 = @SVector zeros(Float64, 3)
        s2 = @SVector zeros(Float64, 6)
        s3 = @SVector zeros(Float64, 9)

        # Dummy constructor 
        funs = FrameAxesFunctions{Int,3}()
        @test typeof(funs) == FrameAxesFunctions{Int,3,9}
        @test length(funs.fun) == 3

        for i in 1:3
            @test funs.fun[i](0, s3, s3) == Rotation{3}(1I)
        end

        # Default constructor 
        faf2 = FrameAxesFunctions{Float64}(_get_fixedrot, fcn_rot2)
        @test typeof(faf2) == FrameAxesFunctions{Float64,2,6}

        @test faf2[2](π / 3, s2, s2)[1] ≈ fcn_rot2(π / 3, s2, s2)[1] atol = atol
        @test faf2[2](π / 3, s2, s2)[2] ≈ fcn_rot2(π / 3, s2, s2)[2] atol = atol

        # Filtering constructor 
        @test_throws ArgumentError FrameAxesFunctions{Float64,2}(_get_fixedrot)

        faf3 = FrameAxesFunctions{Float64,3}(fcn_rot3, fcn_rot3, _get_fixedrot)
        @test typeof(faf3) == FrameAxesFunctions{Float64,3,9}

        faf1 = FrameAxesFunctions{Float64,1}(fcn_rot1, _get_fixedrot, fcn_rot1)
        @test typeof(faf1) == FrameAxesFunctions{Float64,1,3}
        @test faf1[1](π / 3, s1, s1)[1] ≈ angle_to_dcm(π / 3, :Z) atol = atol
    end

    @testset "PointNodes" begin

        # FramePointFunctions 
        # -------------------

        fcn_stv1!(x, t) = x[1:3] .= [1.0, 2.0, t]
        fcn_stv3!(x, t) = x[[4, 2, 8]] .= [6, t, -1.0]

        # Dummy constructor 
        funs = FramePointFunctions{Int,3}()
        @test typeof(funs) == FramePointFunctions{Int,3,9}
        @test length(funs.fun) == 3

        z = @MVector zeros(Int, 9)

        for i in 1:3
            @test typeof(funs.fun[i](z, 0)) == Nothing
        end

        # Default constructor 
        fpf2 = FramePointFunctions{Float64}(_empty_stv_update!, fcn_stv1!)
        @test typeof(fpf2) == FramePointFunctions{Float64,2,6}

        x = @MVector zeros(Float64, 6)

        # Stored functions testing 
        fpf2[1](x, 4.0)
        @test x[1:3] == zeros(3)

        fpf2[2](x, 4.0)
        @test x[1:3] == [1.0, 2.0, 4.0]

        # Filtering constructor 
        @test_throws ArgumentError FramePointFunctions{Float64,2}(_empty_stv_update!)

        fpf3 = FramePointFunctions{Float64,3}(fcn_stv1!, fcn_stv3!, _empty_stv_update!)
        @test typeof(fpf3) == FramePointFunctions{Float64,3,9}

        fpf1 = FramePointFunctions{Float64,1}(fcn_stv1!, _empty_stv_update!, fcn_stv3!)
        @test typeof(fpf1) == FramePointFunctions{Float64,1,3}

        y = @MVector zeros(Float64, 3)
        fpf1[1](y, 2.0)
        @test y[1:3] == [1.0, 2.0, 2.0]
    end

    @testset "FrameSystem" begin

        # FrameSystem Constructors 
        # ------------------------

        # Default constructors
        fs1 = FrameSystem{1,Int64}()
        @test typeof(fs1) ==
            FrameSystem{1,Int64,BarycentricDynamicalTime,NullEphemerisProvider,3}

        @test_throws ArgumentError FrameSystem{0,Float64}()
        @test_throws ArgumentError FrameSystem{5,Float64}()

        fs2 = FrameSystem{3,Float64,Tempo.InternationalAtomicTime}()
        @test typeof(fs2) ==
            FrameSystem{3,Float64,InternationalAtomicTime,NullEphemerisProvider,9}

        # Utilities 
        # ---------

        fs3 = FrameSystem{4,Float64}()
        @test frames_order(fs3) == 4
        @test frames_timescale(fs3) == BarycentricDynamicalTime

        @test has_point(fs3, 20) == false
        @test has_axes(fs3, 3) == false

        @test ephemeris_points(fs3) == Int64[]
    end
end;
