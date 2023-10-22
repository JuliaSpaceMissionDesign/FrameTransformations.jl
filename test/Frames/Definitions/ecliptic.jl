
# insert ICRF
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 22
@axes ECLIPJ2000 17
@axes ICRF_TEST 1
@axes MEME_TEST -10000000
@axes MOD 18

@testset "MEME 2000" verbose = false begin
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, MEME_TEST)

    # Check that you can't add meme2000 to a random axes
    @test_throws ArgumentError add_axes_meme2000!(frames, MEME2000, MEME_TEST)

    # Check that you can add MEME to ICRF with true ID but different name 
    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, ICRF_TEST)
    add_axes_meme2000!(frames, MEME2000, ICRF_TEST)

    # Test rotation matrix from ICRF to MEME2000
    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, ICRF)
    add_axes_meme2000!(frames, MEME2000, ICRF)

    R = rotation9(frames, ICRF, MEME2000, rand())

    v = rand(BigFloat, 3)
    v /= norm(v)

    @test v2as(R[1] * v, Orient.DCM_ICRF_TO_J2000_BIAS * v) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[2])) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[3])) ≈ 0.0 atol = 1e-14 rtol = 1e-14

    # Test rotation matrix from ECLIPJ2000 to MEME2000
    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, ECLIPJ2000)
    add_axes_meme2000!(frames, MEME2000, ECLIPJ2000)

    R = rotation9(frames, MEME2000, ECLIPJ2000, rand())

    v = rand(BigFloat, 3)
    v /= norm(v)

    @test v2as(R[1] * v, Orient.DCM_J2000_TO_ECLIPJ2000 * v) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[2])) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[3])) ≈ 0.0 atol = 1e-14 rtol = 1e-14

    frames = FrameSystem{3,Float64}()
    add_axes_icrf!(frames)

    # Check warning on axesid not being the true MEME2000 ID
    @test_logs (:warn, "ECLIPJ2000 is aliasing an ID that is not the standard MEME2000 ID" *
    " ($(Orient.AXESID_MEME2000)).") add_axes_meme2000!(frames, ECLIPJ2000, ICRF)

end;

@testset "Ecliptic Equinox at J2000" verbose = false begin
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, MEME_TEST)

    # Check that you can't add eclipj2000 to a random axes
    @test_throws ArgumentError add_axes_eclipj2000!(frames, ECLIPJ2000, MEME_TEST)

    f1 = FrameSystem{3,Float64}()
    f2 = FrameSystem{3,Float64}()

    # Test Rotation from MEME2000 to ECLIPJ2000
    # Note: SPICE defines this rotation using the obliquity of 
    # the ecliptic taken from the IAU 1976/1980 Theory of Nutation
    add_axes_inertial!(f1, ICRF)
    add_axes_meme2000!(f1, MEME2000, ICRF)
    add_axes_eclipj2000!(f1, ECLIPJ2000, ICRF)

    add_axes_inertial!(f2, MEME2000)
    add_axes_eclipj2000!(f2, ECLIPJ2000, MEME2000)

    # Test that ECLIPJ2000 is defined correctly with both ICRF and 
    # MEME2000 as parents!
    for frames in (f1, f2)
        ep = rand(0.0:1e7)

        R = sxform("J2000", "ECLIPJ2000", ep)
        R_ = rotation6(frames, MEME2000, ECLIPJ2000, ep)

        v = rand(BigFloat, 3)
        v /= norm(v)

        @test v2as(R[1:3, 1:3] * v, R_[1] * v) ≈ 0.0 atol = 1e-6
        @test R_[2] ≈ zeros(3, 3) atol = 1e-14
    end

    frames = FrameSystem{3,Float64}()
    add_axes_icrf!(frames)

    # Check warning on axesid not being the true MEME2000 ID
    @test_logs (:warn, "MEME2000 is aliasing an ID that is not the standard ECLIPJ2000 ID" *
    " ($(Orient.AXESID_ECLIPJ2000)).") add_axes_eclipj2000!(frames, MEME2000, ICRF)

end;

@testset "Mean of Date Ecliptic Equinox" verbose=false begin 


    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, MEME_TEST)

    # Check that you can add MOD only with respect to the ICRF 
    @test_throws ArgumentError add_axes_mememod!(frames, MOD, MEME_TEST)

    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, ICRF)
    add_axes_mememod!(frames, MOD, ICRF)

    # Test that the MEMEMOD transformation is defined correctly 
    ep = rand(0.0:1e7)
    R = rotation6(frames, MOD, ICRF, ep)
    @test R[2] ≈ zeros(3, 3) atol=1e-14

    R_ =  Orient.orient_rot3_icrf_to_mememod(ep);

    v = rand(BigFloat, 3)
    v /= norm(v)

    @test v2as(R[1] * v, R_' * v) ≈ 0.0 atol = 1e-14 rtol = 1e-14

end;
