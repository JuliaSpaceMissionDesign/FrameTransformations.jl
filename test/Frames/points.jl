kclear()

@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 22
@axes ECLIPJ2000 17

@point SSB 0 SolarSystemBarycenter
@point Sun 10
@point EMB 3
@point Earth 399
@point Moon 301
@point Jupiter 599

@testset "Points Definitions" verbose = true begin
    @testset "Macro @point" verbose = false begin
        @test Frames.point_id(SSB) == 0

        @test point_alias(SSB) == 0
        @test point_alias(0) == 0
        @test point_alias(Earth) == 399

        @test typeof(SSB) == SolarSystemBarycenterPoint
        @test typeof(Sun) == SunPoint

        @test isa(SSB, Frames.AbstractFramePoint)
        @test isa(Sun, Frames.AbstractFramePoint)
    end

    # -- Testing ROOT POINT
    @testset "Root" verbose = false begin
        frames = FrameSystem{3,Float64}()

        # test unknown axes 
        @test_throws ArgumentError add_point_root!(frames, SSB, ICRF)
        add_axes_inertial!(frames, ICRF)
        add_point_root!(frames, SSB, ICRF)

        # test point properties
        node = frames_points(frames).nodes[1]
        @test node.class == :RootPoint
        @test length(node.stv) == 1
        @test length(node.nzo) == 1

        @test node.name == :SSB
        @test node.axesid == 1
        @test node.NAIFId == 0

        # test root point already exists
        @test_throws ArgumentError add_point_root!(frames, Sun, ICRF)
    end

    # -- Testing Fixed points!
    @testset "Fixed" verbose = false begin
        frames = FrameSystem{3,Float64}()
        add_axes_inertial!(frames, ICRF)
        add_point_root!(frames, SSB, ICRF)

        # test wrong offset dimension (should be 3)
        offset = rand(5)
        @test_throws DimensionMismatch add_point_fixed!(frames, Earth, SSB, ICRF, offset)

        offset = rand(3)

        # test unknown parent point 
        @test_throws ArgumentError add_point_fixed!(frames, Earth, Sun, ICRF, offset)
        @test_throws ArgumentError add_point_fixed!(frames, Earth, SSB, MEME2000, offset)

        add_point_fixed!(frames, Earth, SSB, ICRF, offset)

        # test point properties 
        node = frames_points(frames).nodes[2]
        @test node.class == :FixedPoint
        @test length(node.stv) == 1
        @test length(node.nzo) == 1
        @test node.name == :Earth
        @test node.axesid == 1
        @test node.NAIFId == 399

        # test point transformation 
        x = vector3(frames, Earth, SSB, ICRF, rand())
        @test x ≈ -offset atol = 1e-12 rtol = 1e-12

        x = vector6(frames, Earth, SSB, ICRF, rand())
        @test x[1:3] ≈ -offset atol = 1e-12 rtol = 1e-12
        @test x[4:6] == zeros(3)

        x = vector9(frames, Earth, SSB, ICRF, rand())
        @test x[1:3] ≈ -offset atol = 1e-12 rtol = 1e-12
        @test x[4:9] == zeros(6)
    end

    # -- Testing Ephemeris points!
    @testset "Ephemeris" verbose = false begin
        furnsh(path(KERNELS[:LEAP]), path(KERNELS[:DE432]))

        frames = FrameSystem{3,Float64}()
        # test ephemeris data is not available 
        @test_throws ErrorException add_point_ephemeris!(frames, Sun)

        eph = EphemerisProvider(path(KERNELS[:DE432]))
        frames = FrameSystem{3,Float64}(eph)

        add_axes_inertial!(frames, ICRF)
        add_axes_fixedoffset!(frames, ECLIPJ2000, ICRF, Orient.DCM_J2000_TO_ECLIPJ2000)

        add_point_root!(frames, SSB, ICRF)
        add_point_ephemeris!(frames, Sun)

        # test point already registered
        @test_throws ArgumentError add_point_ephemeris!(frames, Sun)
        # test cannot add point because Earth is defined wrt the EMB (parent not registered)
        @test_throws ErrorException add_point_ephemeris!(frames, Earth)

        add_point_ephemeris!(frames, EMB)
        add_point_ephemeris!(frames, Earth)
        add_point_ephemeris!(frames, Moon)

        add_point_fixed!(frames, Jupiter, SSB, ICRF, rand(3))

        # test Sun point properties 
        node = frames_points(frames).nodes[2]
        @test node.class == :EphemerisPoint
        @test node.axesid == 1
        @test node.NAIFId == 10
        @test node.parentid == 0
        @test node.name == :Sun

        nth = Threads.nthreads()
        @test length(node.nzo) == nth
        @test length(node.epochs.du) == nth
        @test length(node.stv) == nth

        # test Earth point properties
        node = frames_points(frames).nodes[4]
        @test node.axesid == 1
        @test node.NAIFId == 399
        @test node.parentid == 3
        @test node.name == :Earth

        # test Moon point properties 
        node = frames_points(frames).nodes[5]
        @test node.axesid == 1
        @test node.parentid == 3
        @test node.name == :Moon

        # test transformation 
        for _ in 1:25
            et = rand(0.0:1e8)
            # test transform with pos 
            x = vector3(frames, Sun, SSB, ICRF, et)
            y = spkpos("SSB", et, "J2000", "NONE", "SUN")[1]
            @test x ≈ y atol = 1e-12 rtol = 1e-12

            # test transform pos\vel 
            xx = vector6(frames, SSB, Sun, ECLIPJ2000, et)
            yy = spkezr("SUN", et, "ECLIPJ2000", "NONE", "SSB")[1]
            @test xx ≈ yy atol = 1e-12 rtol = 1e-12

            # test multiple transform 
            xx = vector6(frames, Moon, Sun, ICRF, et)
            yy = spkezr("SUN", et, "J2000", "NONE", "MOON")[1]
            @test xx ≈ yy atol = 1e-12 rtol = 1e-12

            xx = vector6(frames, EMB, Sun, ECLIPJ2000, et)
            yy = spkezr("SUN", et, "ECLIPJ2000", "NONE", "EMB")[1]
            @test xx ≈ yy atol = 1e-12 rtol = 1e-12
        end
        
        # Check that the epoch must have the same timescale
        tai = Epoch("$(rand(0.0:1e8)) TAI")

        @test_throws ArgumentError vector3(frames, Moon, Sun, ICRF, tai)

        frames2 = FrameSystem{3,Float64}(eph)
        add_axes_inertial!(frames2, MEME2000)
        add_point_root!(frames2, SSB, MEME2000)

        # test ephemeris axes have not been registered
        @test_throws ErrorException add_point_ephemeris!(frames2, Sun)

        # test that derivatives up to order 4 work 
        frames = FrameSystem{4,Float64}(eph)

        add_axes_inertial!(frames, ICRF)
        add_axes_fixedoffset!(frames, ECLIPJ2000, ICRF, Orient.DCM_J2000_TO_ECLIPJ2000)

        add_point_root!(frames, SSB, ICRF)
        add_point_ephemeris!(frames, Sun)

        # check that the derivatives work up to jerk!
        vector12(frames, SSB, Sun, ICRF, 0.0)
    end

    # -- Testing Dynamical points!
    @testset "Dynamical" verbose = false begin
        R = Orient.DCM_J2000_TO_ECLIPJ2000

        nth = Threads.nthreads()
        tid = Threads.threadid()

        # Numerical functions
        f = t -> SA[sin(t) * ones(3)...]
        df = t -> SA[sin(t) * ones(3)..., cos(t) * ones(3)...]
        ddf = t -> SA[sin(t) * ones(3)..., cos(t) * ones(3)..., -sin(t) * ones(3)...]
        dddf =
            t -> SA[
                sin(t) * ones(3)...,
                cos(t) * ones(3)...,
                -sin(t) * ones(3)...,
                -cos(t) * ones(3)...,
            ]

        frames = FrameSystem{4,Float64}()
        add_axes_inertial!(frames, ICRF)
        add_point_root!(frames, EMB, ICRF)
        add_point_dynamical!(frames, Earth, EMB, ICRF, f)

        # test point properties 
        node = frames_points(frames).nodes[2]

        @test node.class == :DynamicalPoint
        @test node.name == :Earth
        @test node.NAIFId == 399
        @test node.parentid == 3
        @test node.axesid == 1

        @test length(node.stv) == nth
        @test length(node.epochs.du) == nth
        @test length(node.nzo) == nth
        @test node.nzo[1][1] == -1
        @test node.nzo[1][2] == -1

        atol, rtol = 1e-12, 1e-12
        # test AD derivatives for all combinations of specified functions
        for funs in ((f,), (f, df), (f, df, ddf), (f, df, ddf, dddf))
            frames = FrameSystem{4,Float64}()
            add_axes_inertial!(frames, ICRF)
            add_axes_fixedoffset!(frames, ECLIPJ2000, ICRF, R)
            add_point_root!(frames, EMB, ICRF)

            add_point_dynamical!(frames, Earth, EMB, ICRF, funs...)

            for _ in 1:5
                ep = rand()

                y = @MVector zeros(12)
                node.f.fun[1](y, ep)
                @test y[1:3] ≈ f(ep) atol = atol rtol = rtol

                node.f.fun[2](y, ep)
                @test y[1:6] ≈ df(ep) atol = atol rtol = rtol

                node.f.fun[3](y, ep)
                @test y[1:9] ≈ ddf(ep) atol = atol rtol = rtol

                node.f.fun[4](y, ep)
                @test y[1:12] ≈ dddf(ep) atol = atol rtol = rtol

                # test transformation with rotations to new frame!
                x = vector6(frames, Earth, EMB, ECLIPJ2000, ep)
                @test x[1:3] ≈ -R * sin(ep) * ones(3) atol = atol rtol = rtol
                @test x[4:6] ≈ -R * cos(ep) * ones(3) atol = atol rtol = rtol

                # test transformation with rotations to existing frame! 
                x = vector12(frames, EMB, Earth, ICRF, ep)
                @test x[10:12] ≈ -cos(ep) * ones(3) atol = atol rtol = rtol
            end
        end
    end

    # -- Testing Updatable points!
    @testset "Updatable" verbose = false begin
        R = Orient.DCM_J2000_TO_ECLIPJ2000

        frames = FrameSystem{3,Float64}()
        add_axes_inertial!(frames, ICRF)
        add_axes_fixedoffset!(frames, ECLIPJ2000, ICRF, R)

        add_point_root!(frames, EMB, ICRF)
        add_point_updatable!(frames, Earth, EMB, ECLIPJ2000)

        node = frames_points(frames).nodes[2]

        nth = Threads.nthreads()
        tid = Threads.threadid()

        # test point properties 
        @test node.class == :UpdatablePoint
        @test node.NAIFId == 399
        @test node.parentid == 3
        @test node.name == :Earth
        @test node.axesid == 17

        @test length(node.stv) == nth
        @test length(node.epochs.du) == nth
        @test length(node.nzo) == nth
        @test node.nzo[1][1] == -1
        @test node.nzo[1][2] == -1

        # test point not updated
        @test_throws ErrorException vector3(frames, Earth, EMB, ICRF, rand())
        # test point not in frame system 
        @test_throws ArgumentError update_point!(frames, Sun, rand(3), rand())
        # test state vector is greater than order 
        @test_throws ArgumentError update_point!(frames, Sun, rand(12), rand())
        # test state vector is not multiple of 3 
        @test_throws ArgumentError update_point!(frames, Earth, rand(7), rand())
        # test point is not updatable 
        @test_throws ArgumentError update_point!(frames, EMB, rand(6), rand())

        r = rand(3)
        rr = rand(6)

        update_point!(frames, Earth, r, 0.0)
        @test node.nzo[tid][1] == 1
        @test get_tmp(node.epochs, 0.0)[tid] == 0.0
        @test get_tmp(node.stv[tid], 0.0)[1:3] == r

        # test transformation at different epoch 
        @test_throws ErrorException vector3(frames, Earth, EMB, ECLIPJ2000, 1.0)
        # test transformation at different order 
        @test_throws ErrorException vector6(frames, Earth, EMB, ICRF, 0.0)

        # test transformation 
        @test vector3(frames, Earth, EMB, ECLIPJ2000, 0.0) == -r

        update_point!(frames, Earth, rr, 1.0)

        @test node.nzo[tid][1] == 2
        @test get_tmp(node.epochs, 1.0)[tid] == 1.0

        # test epoch has been updated correctly 
        @test_throws ErrorException vector3(frames, EMB, Earth, ICRF, 0.0)

        # test transformation with rotation! 
        x = vector6(frames, EMB, Earth, ICRF, 1.0)
        y = vcat(R' * rr[1:3], R' * rr[4:6])
        @test x ≈ y atol = 1e-12 rtol = 1e-12
    end
end;

kclear()
