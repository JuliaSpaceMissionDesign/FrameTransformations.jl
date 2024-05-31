using Test
using SPICE 
using LinearAlgebra
using FrameTransformations
using ReferenceFrameRotations

v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

ICRF = 1
EME2000 = 22
ECL2000 = 17
ICRF_TEST = 1
EME_TEST = -10000000

@testset "EME 2000" verbose = false begin
    frames = FrameSystem{3, Float64}()
    add_axes_root!(frames, :Test, EME_TEST)

    # Check that you can't add meme2000 to a random axes
    @test_throws ArgumentError add_axes_eme2000!(frames, :EME, EME_TEST)

    # Test rotation matrix from ICRF to EME2000
    frames = FrameSystem{3,Float64}()
    add_axes_icrf!(frames)
    add_axes_eme2000!(frames, :EME)

    R = rotation9(frames, ICRF, EME2000, rand())
    v = rand(BigFloat, 3)
    v /= norm(v)

    @test v2as(R[1] * v, FrameTransformations.DCM_ICRF_TO_EME2000 * v) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[2])) ≈ 0.0 atol = 1e-14 rtol = 1e-14
    @test maximum(abs.(R[3])) ≈ 0.0 atol = 1e-14 rtol = 1e-14
end

@testset "ECL 2000" verbose = false begin

    frames = FrameSystem{3, Float64}()
    add_axes_root!(frames, :Test, EME_TEST)

    # Check that you can't add eclipj2000 to a random axes
    @test_throws ArgumentError add_axes_ecl2000!(frames, :ECL, EME_TEST)

    f1 = FrameSystem{3,Float64}()
    f2 = FrameSystem{3,Float64}()

    # Test Rotation from EME2000 to ECL2000
    # Note: SPICE defines this rotation using the obliquity of 
    # the ecliptic taken from the IAU 1976/1980 Theory of Nutation
    add_axes_icrf!(f1)
    add_axes_eme2000!(f1, :EME)
    add_axes_ecl2000!(f1, :ECL)

    add_axes_root!(f2, :EME, EME2000) 
    add_axes_ecl2000!(f2, :ECL, EME2000, ECL2000)

    # Test that ECL2000 is defined correctly with both ICRF and 
    # EME2000 as parents!
    for frames in (f1, f2)
        ep = rand(0.0:1e7)

        R = sxform("J2000", "ECLIPJ2000", ep)
        R_ = rotation6(frames, EME2000, ECL2000, ep)

        v = rand(BigFloat, 3)
        v /= norm(v)

        @test v2as(R[1:3, 1:3] * v, R_[1] * v) ≈ 0.0 atol = 1e-6
        @test R_[2] ≈ zeros(3, 3) atol = 1e-14
    end

end