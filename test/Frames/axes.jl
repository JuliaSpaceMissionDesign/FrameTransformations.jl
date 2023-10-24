kclear()

@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 22
@axes ECLIPJ2000 17

@axes AXES_ROT 2
@axes AXES_EPHEM 31008
@axes AXES_COMP 3

@point SSB 0 SolarSystemBarycenter
@point Sun 10
@point Earth 399

@testset "Axes Definitions" verbose = true begin
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600
    v = rand(BigFloat, 3)

    # -- Testing @axes macro and type definition
    @testset "Macro @axes" verbose = true begin
        @test Frames.axes_id(ICRF) == 1
        @test axes_alias(ICRF) == 1
        @test axes_alias(414) == 414

        @test typeof(ICRF) == InternationalCelestialReferenceFrameAxes
        @test typeof(AXES_ROT) == AxesRotAxes
        @test typeof(ECLIPJ2000) == ECLIPJ2000Axes

        @test isa(ICRF, Frames.AbstractFrameAxes)
    end

    # -- Testing INERTIAL AXES 
    frames = FrameSystem{3,Float64}()
    R = angle_to_dcm(π / 7, :Y)

    @testset "Inertial" verbose = false begin
        @test_throws ArgumentError add_axes_inertial!(frames, ICRF; parent=MEME2000)

        add_axes_inertial!(frames, ICRF)
        @test frames_axes(frames).nodes[1].name == :ICRF
        @test is_timefixed(frames, ICRF)
        @test is_inertial(frames, ICRF)

        # test missing parent 
        @test_throws ArgumentError add_axes_inertial!(frames, MEME2000)
        # test missing parent but with DCM 
        @test_throws ArgumentError add_axes_inertial!(frames, MEME2000; dcm=DCM(1.0I))
        # test parent not in frames 
        @test_throws ArgumentError add_axes_inertial!(frames, MEME2000; parent=AXES_ROT)
        # test missing DCM 
        @test_throws ArgumentError add_axes_inertial!(frames, MEME2000; parent=ICRF)

        G = FrameSystem{2, Float64}()
        add_axes_inertial!(G, ICRF)
        add_axes_rotating!(G, AXES_ROT, ICRF, t->angle_to_dcm(t, :Z))

        # test insufficient frame system order 
        @test_throws ErrorException rotation9(G, ICRF, AXES_ROT, 0.0)
        
        # test parent must be inertial 
        @test_throws ArgumentError add_axes_inertial!(G, MEME2000; parent=AXES_ROT, dcm=DCM(1.0I))

        # Check that if the axes are not registered an error is thrown 
        @test_throws ErrorException rotation9(frames, ICRF, MEME2000, 0.0)
        add_axes_inertial!(frames, MEME2000; parent=ICRF, dcm=R)

        # Test actual rotation 
        Rb = rotation9(frames, ICRF, MEME2000, rand())
        @test v2as(Rb[1] * v, R * v) ≤ 1e-6

        for i in 2:3
            @test maximum(abs.(Rb[i])) == 0.0
        end

        # Test axes properties 
        for i in 1:2
            node = frames_axes(frames).nodes[i]

            @test node.class == :InertialAxes
            @test length(node.angles) == 1

        end

        @test is_timefixed(frames, MEME2000)
        @test is_inertial(frames, MEME2000)
    end

    A = angle_to_dcm(π / 3, :Z)

    # -- Testing FIXEDOFFSET AXES
    @testset "FixedOffset" verbose = false begin

        # test axes are already registered 
        @test_throws ArgumentError add_axes_fixedoffset!(frames, MEME2000, ICRF, R)
        @test_throws ArgumentError add_axes_fixedoffset!(frames, ECLIPJ2000, AXES_ROT, R)

        add_axes_fixedoffset!(frames, ECLIPJ2000, MEME2000, A)

        # Test actual rotation 
        Rb = rotation9(frames, ICRF, ECLIPJ2000, rand())
        @test v2as(A * R * v, Rb[1] * v) ≤ 1e-6

        for i in 2:3
            @test maximum(abs.(Rb[i])) == 0.0
        end

        @test is_timefixed(frames, ECLIPJ2000)
        @test is_inertial(frames, ECLIPJ2000)

        node = frames_axes(frames).nodes[end]
        @test node.class == :FixedOffsetAxes
        @test node.id == 17
        @test node.name == :ECLIPJ2000
    end

    # -- Testing ROTATING AXES 
    @testset "Rotating" verbose = false begin
        smb = rand([:X, :Y, :Z])
        f = t -> angle_to_dcm(t, smb)

        df = t -> (angle_to_dcm(t, smb), angle_to_δdcm([t, 1], smb))

        ddf =
            t -> (
                angle_to_dcm(t, smb),
                angle_to_δdcm([t, 1], smb),
                angle_to_δ²dcm([t, 1, 0], smb),
            )

        dddf =
            t -> (
                angle_to_dcm(t, smb),
                angle_to_δdcm([t, 1], smb),
                angle_to_δ²dcm([t, 1, 0], smb),
                angle_to_δ³dcm([t, 1, 0, 0], smb),
            )

        nth = Threads.nthreads()
        tid = Threads.threadid()

        G = FrameSystem{4,Float64}()
        add_axes_inertial!(G, ICRF)

        # test axes properties 
        add_axes_rotating!(G, AXES_ROT, ICRF, f)
        node = frames_axes(G).nodes[2]

        @test node.class == :RotatingAxes
        @test node.name == :AXES_ROT
        @test node.id == 2
        @test node.parentid == 1

        # Check axes ineritality
        @test is_inertial(G, AXES_ROT) == false
        @test is_inertial(G, 2) == false
        @test is_timefixed(G, AXES_ROT) == false
        @test is_timefixed(G, 2) == false

        B = Orient.DCM_J2000_TO_ECLIPJ2000

        atol, rtol = 1e-12, 1e-12
        # test AD derivatives for all combinations of specified functions
        for funs in ((f,), (f, df), (f, df, ddf), (f, df, ddf, dddf))
            G = FrameSystem{4,Float64}()

            add_axes_inertial!(G, ICRF)
            add_axes_rotating!(G, AXES_ROT, ICRF, funs...)
            add_axes_inertial!(G, ECLIPJ2000; parent=ICRF, dcm=B)

            x = @SVector zeros(12)

            for _ in 1:5
                ep = rand()
                tdb = Epoch("$(ep/86400) TDB")
                tai = Epoch("$ep TAI")

                R = node.f.fun[1](ep, x, x)
                @test typeof(R) == Rotation{4,Float64}
                @test R[1] ≈ angle_to_dcm(ep, smb) atol = atol rtol = rtol

                R = node.f.fun[2](ep, x, x)
                @test R[1] ≈ angle_to_dcm(ep, smb) atol = atol rtol = rtol
                @test R[2] ≈ angle_to_δdcm([ep, 1], smb) atol = atol rtol = rtol

                R = node.f.fun[3](ep, x, x)
                @test R[1] ≈ angle_to_dcm(ep, smb) atol = atol rtol = rtol
                @test R[2] ≈ angle_to_δdcm([ep, 1], smb) atol = atol rtol = rtol
                @test R[3] ≈ angle_to_δ²dcm([ep, 1, 0], smb) atol = atol rtol = rtol

                R = node.f.fun[4](ep, x, x)
                @test R[1] ≈ angle_to_dcm(ep, smb) atol = atol rtol = rtol
                @test R[2] ≈ angle_to_δdcm([ep, 1], smb) atol = atol rtol = rtol
                @test R[3] ≈ angle_to_δ²dcm([ep, 1, 0], smb) atol = atol rtol = rtol
                @test R[4] ≈ angle_to_δ³dcm([ep, 1, 0, 0], smb) atol = atol rtol = rtol

                # test throws error for different epoch timescale 
                @test_throws ArgumentError rotation3(G, AXES_ROT, ICRF, tai)
                
                # test transformations 
                A = rotation6(G, AXES_ROT, ICRF, ep)
                A2 = rotation6(G, AXES_ROT, ICRF, tdb)
                
                @test A[1] ≈ angle_to_dcm(ep, smb)' atol = atol rtol = rtol
                @test A[2] ≈ angle_to_δdcm([ep, 1], smb)' atol = atol rtol = rtol
                @test A[1] ≈ A2[1] atol=atol rtol=rtol 
                @test A[2] ≈ A2[2] atol=atol rtol=rtol

                A = rotation12(G, AXES_ROT, ICRF, ep)
                @test A[3] ≈ angle_to_δ²dcm([ep, 1, 0], smb)' atol = atol rtol = rtol
                @test A[4] ≈ angle_to_δ³dcm([ep, 1, 0, 0], smb)' atol = atol rtol = rtol

                # test transformation with another frame 
                A = rotation6(G, AXES_ROT, ECLIPJ2000, ep)
                @test A[1] ≈ B * angle_to_dcm(ep, smb)' atol = atol rtol = rtol
                @test A[2] ≈ B * angle_to_δdcm([ep, 1], smb)' atol = atol rtol = rtol
            end
        end
    end

    # -- Testing COMPUTABLE AXES 
    @testset "Computable" verbose = false begin
        nth = Threads.nthreads()
        tid = Threads.threadid()

        rfun = t -> SA[sin(t), cos(t), t]
        vfun = t -> SA[cos(t), -sin(t), 1]
        afun = t -> SA[-sin(t), -cos(t), 0]
        jfun = t -> SA[-cos(t), sin(t), 0]

        eph = EphemerisProvider(path(KERNELS[:DE432]))

        B = Orient.DCM_J2000_TO_ECLIPJ2000

        F = FrameSystem{4,Float64}(eph)
        add_axes_inertial!(F, ICRF)
        add_axes_fixedoffset!(F, ECLIPJ2000, ICRF, B)

        v1 = ComputableAxesVector(Earth, Sun, 1)
        v2 = ComputableAxesVector(Earth, Sun, 2)

        # test wrong rotation sequence
        @test_throws ArgumentError add_axes_computable!(
            F, AXES_COMP, ECLIPJ2000, v1, v2, :Z
        )

        # test points not defined 
        @test_throws ArgumentError add_axes_computable!(
            F, AXES_COMP, ECLIPJ2000, v1, v2, :XY
        )

        add_point_root!(F, Earth, ICRF)
        add_point_dynamical!(F, Sun, Earth, ICRF, rfun)

        seq = rand([:XY, :YX, :XZ, :ZX, :YZ, :ZY])

        add_axes_computable!(F, AXES_COMP, ECLIPJ2000, v1, v2, seq)

        # test axes properties 
        node = frames_axes(F).nodes[3]
        @test node.name == :AXES_COMP
        @test node.id == 3
        @test node.parentid == 17
        @test node.class == :ComputableAxes

        @test length(node.angles) == nth

        @test is_inertial(F, AXES_COMP) == false
        @test is_timefixed(F, AXES_COMP) == false

        # test order is not enough
        @test_throws ErrorException rotation12(F, AXES_COMP, ECLIPJ2000, 0.0)

        atol, rtol = 1e-12, 1e-12

        # test rotations
        for _ in 1:10
            ep = rand()
            # These are rotated first to ECLIPJ2000, because its the parent frame of AXES_COMP!
            r, v, a, j = B * rfun(ep), B * vfun(ep), B * afun(ep), B * jfun(ep)

            A = Frames.twovectors_to_dcm(r, v, seq)
            dA = Frames.twovectors_to_δdcm(vcat(r, v), vcat(v, a), seq)
            ddA = Frames.twovectors_to_δ²dcm(vcat(r, v, a), vcat(v, a, j), seq)

            R3 = rotation3(F, AXES_COMP, ICRF, ep)
            @test R3[1] ≈ (A * B)' atol = atol rtol = rtol

            R6 = rotation6(F, AXES_COMP, ICRF, ep)
            @test R3[1] ≈ (A * B)' atol = atol rtol = rtol
            @test R6[2] ≈ (dA * B)' atol = atol rtol = rtol

            R9 = rotation9(F, AXES_COMP, ICRF, ep)
            @test R3[1] ≈ (A * B)' atol = atol rtol = rtol
            @test R9[2] ≈ (dA * B)' atol = atol rtol = rtol
            @test R9[3] ≈ (ddA * B)' atol = atol rtol = rtol
        end
    end

    # -- Testing EPHEMERIS AXES 
    @testset "Ephemeris" verbose = false begin
        nth = Threads.nthreads()
        tid = Threads.threadid()

        F = FrameSystem{2,Float64}()

        # test parent axes not defined
        @test_throws ArgumentError add_axes_ephemeris!(F, AXES_EPHEM, :ZYX)

        add_axes_inertial!(F, ICRF)

        # test ephem data not available for that ID
        @test_throws ArgumentError add_axes_ephemeris!(F, MEME2000, :ZYX)

        # test invalid rotation sequence 
        @test_throws ArgumentError add_axes_ephemeris!(F, MEME2000, :ZY)

        # Load kernels!
        furnsh(path(KERNELS[:LEAP]), path(KERNELS[:PA440]), path(KERNELS[:FK_DE440]))

        eph = EphemerisProvider(path(KERNELS[:PA440]))
        G = FrameSystem{4,Float64}(eph)

        add_axes_inertial!(G, ICRF)
        add_axes_ephemeris!(G, AXES_EPHEM, :ZXZ)

        # test axes properties 
        node = frames_axes(G).nodes[2]
        @test node.name == :AXES_EPHEM
        @test node.id == 31008
        @test node.parentid == 1
        @test node.class == :EphemerisAxes

        @test length(node.angles) == nth

        @test is_inertial(G, AXES_EPHEM) == false
        @test is_timefixed(G, AXES_EPHEM) == false

        # test rotation 
        atol, rtol = 1e-12, 1e-12
        y = @MVector zeros(12)

        for _ in 1:10
            et = rand(0.0:1e5)
            Rs = sxform("MOON_PA", "J2000", et)

            angles = ephem_orient!(y, eph, DJ2000, et / Tempo.DAY2SEC, 31008, 1, 3)
            ddR = Math._3angles_to_δ²dcm(y, :ZXZ)
            dddR = Math._3angles_to_δ³dcm(y, :ZXZ)

            R3 = rotation3(G, AXES_EPHEM, ICRF, et)
            @test R3[1] ≈ Rs[1:3, 1:3] atol = atol rtol = rtol

            R6 = rotation6(G, AXES_EPHEM, ICRF, et)
            @test R6[1] ≈ Rs[1:3, 1:3] atol = atol rtol = rtol
            @test R6[2] ≈ Rs[4:6, 1:3] atol = atol rtol = rtol

            R9 = rotation9(G, AXES_EPHEM, ICRF, et)
            @test R9[1] ≈ Rs[1:3, 1:3] atol = atol rtol = rtol
            @test R9[2] ≈ Rs[4:6, 1:3] atol = atol rtol = rtol
            @test R9[3] ≈ ddR' atol = atol rtol = rtol

            R12 = rotation12(G, ICRF, AXES_EPHEM, et)
            @test R12[1] ≈ Rs[1:3, 1:3]' atol = atol rtol = rtol
            @test R12[2] ≈ Rs[4:6, 1:3]' atol = atol rtol = rtol
            @test R12[3] ≈ ddR atol = atol rtol = rtol
            @test R12[4] ≈ dddR atol = atol rtol = rtol
        end
    end

    # -- Testing PROJECTED AXES 
    @testset "Projected" verbose = false begin
        nth = Threads.nthreads()
        tid = Threads.threadid()

        fun = t -> angle_to_dcm(t, t / 2, :XY)
        fun2 = t -> angle_to_dcm(t, t / 3, t^2, :ZYX)

        G = FrameSystem{4,Float64}()
        add_axes_inertial!(G, ICRF)
        add_axes_projected!(G, MEME2000, ICRF, fun)
        add_axes_projected!(G, ECLIPJ2000, MEME2000, fun2)

        add_axes_fixedoffset!(G, AXES_ROT, ECLIPJ2000, fun(π / 3))

        node = frames_axes(G).nodes[2]

        # test axes properties 
        @test node.name == :MEME2000
        @test node.class == :ProjectedAxes
        @test node.id == 22
        @test node.parentid == 1

        @test length(node.angles) == nth

        @test is_inertial(G, MEME2000)
        @test is_timefixed(G, MEME2000) == false

        @test is_inertial(G, AXES_ROT)
        @test is_timefixed(G, AXES_ROT) == false

        # test rotation
        atol, rtol = 1e-12, 1e-12
        for _ in 1:5
            ep = rand()
            Rₑ = (fun(π / 3) * fun2(ep) * fun(ep))'

            for (i, rot) in enumerate((rotation3, rotation6, rotation9, rotation12))

                # test rotation 
                R = rot(G, AXES_ROT, ICRF, ep)
                @test R[1] ≈ Rₑ atol = atol rtol = rtol

                # test that all rotation derivatives are null
                for j in 2:i
                    @test maximum(abs.(R[i])) == 0.0
                end
            end
        end
    end
end;

kclear()

