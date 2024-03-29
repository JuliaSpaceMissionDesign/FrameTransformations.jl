kclear()

@axes ICRF 1 InternationalCelestialReferenceFrame
@axes ECI 200
@axes ME421 31007
@axes PA421 31006
@axes PA440 31008
@axes PA_TEST 23412

@testset "Moon PA and ME Axes" verbose = true begin
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    @testset "DE421" verbose = false begin

        # Check that if you haven't loaded the kernel you get an error 
        frames = FrameSystem{2,Float64}()
        add_axes_inertial!(frames, ICRF)
        @test_throws ArgumentError add_axes_pa421!(frames, PA421)

        for kernel in (:LEAP, :PA421, :FK_DE421)
            furnsh(path(KERNELS[kernel]))
        end

        # Test orientation between PA421 and ICRF
        eph = EphemerisProvider(path(KERNELS[:PA421]))
        frames = FrameSystem{2,Float64}(eph)

        add_axes_inertial!(frames, ICRF)

        # Test that you cant add a set of axes with wrong ID
        @test_throws ArgumentError add_axes_pa421!(frames, PA_TEST)

        add_axes_pa421!(frames, PA421)

        # Test that you can't add them wrt to the ICRF
        @test_throws ArgumentError add_axes_me421!(frames, ME421, ICRF)

        add_axes_me421!(frames, ME421, PA421)

        for _ in 1:100
            et = rand(0.0:1e8)

            v = rand(BigFloat, 3)
            v /= norm(v)

            # Test PA421!
            Rb = rotation6(frames, PA421, ICRF, et)
            Rs = sxform("MOON_PA", "J2000", et)

            @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
            @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6

            # Test ME421!
            Rb = rotation6(frames, ICRF, ME421, et)
            Rs = sxform("J2000", "MOON_ME", et)

            @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
            @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6
        end

        # Test that you must add the PA421 with respect to the ICRF 
        frames = FrameSystem{2,Float64}(eph)
        add_axes_inertial!(frames, ECI)

        # Test that you can only define it with respect to the ICRF (not registered)
        @test_throws ErrorException add_axes_pa421!(frames, PA421)

        frames = FrameSystem{2,Float64}(eph)
        add_axes_icrf!(frames)

        # Test can you must match the ID of the PA440 axes
        @test_throws ArgumentError add_axes_pa421!(frames, ME421)


    end

    @testset "DE440" verbose = false begin

        # Check that if you haven't loaded the kernel you get an error 
        frames = FrameSystem{2,Float64}()
        add_axes_inertial!(frames, ICRF)
        @test_throws ArgumentError add_axes_pa440!(frames, PA440)

        for kernel in (:LEAP, :PA440, :FK_DE440)
            furnsh(path(KERNELS[kernel]))
        end

        # Test orientation between PA440 and ICRF
        eph = EphemerisProvider(path(KERNELS[:PA440]))
        frames = FrameSystem{2,Float64}(eph)

        add_axes_inertial!(frames, ICRF)

        # Test that you cant add a set of axes with wrong ID
        @test_throws ArgumentError add_axes_pa440!(frames, PA_TEST)

        add_axes_pa440!(frames, PA440)

        # Test that you can't add them wrt to the ICRF
        @test_throws ArgumentError add_axes_me421!(frames, ME421, ICRF)
        add_axes_me421!(frames, ME421, PA440)

        for _ in 1:100
            et = rand(0.0:1e8)

            v = rand(BigFloat, 3)
            v /= norm(v)

            # Test PA421!
            Rb = rotation6(frames, PA440, ICRF, et)
            Rs = sxform("MOON_PA", "J2000", et)

            @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
            @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6

            # Test ME421!
            Rb = rotation6(frames, ICRF, ME421, et)
            Rs = sxform("J2000", "MOON_ME", et)

            @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
            @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6
        end

        # Test that you must add the PA440 with respect to the ICRF 
        frames = FrameSystem{2,Float64}(eph)
        add_axes_inertial!(frames, ECI)

        # Test that you can only define it with respect to the ICRF (not registered)
        @test_throws ErrorException add_axes_pa440!(frames, PA421)

        frames = FrameSystem{2,Float64}(eph)
        add_axes_icrf!(frames)

        # Test can you must match the ID of the PA440 axes
        @test_throws ArgumentError add_axes_pa440!(frames, ME421)

    end
end;

kclear()
