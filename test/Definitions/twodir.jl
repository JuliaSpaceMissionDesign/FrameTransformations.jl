using FrameTransformations
using Ephemerides
using RemoteFiles
using Test
using LinearAlgebra
using StaticArrays

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

atol = 1e-8

# Define two non parallel vectors and their 1st, 2nd and 3rd order time derivatives! 
function get_v1(t)
    return SA[cos(3t), t*sin(t), t^2*cos(t)]
end
get_v2(t) = SA[t^3, t^2, t]

@testset "Sequence assembly" begin
    for _ = 1:100

        θ = rand()
        a, b = get_v1(θ)[1:3], get_v2(θ)[1:3]
        c = cross(a, b)

        aᵤ, bᵤ, cᵤ = a / norm(a), b / norm(b), c / norm(c)

        # XY 
        R = FrameTransformations.twodir_to_dcm(a, b, :XY)

        @test R' * [1, 0, 0] ≈ aᵤ atol = atol
        @test R' * [0, 0, 1] ≈ cᵤ atol = atol
        @test dot(R' * [0, 1, 0], aᵤ) ≈ 0 atol = atol

        # YX 
        R = FrameTransformations.twodir_to_dcm(a, b, :YX)

        @test R' * [0, 1, 0] ≈ aᵤ atol = atol
        @test R' * [0, 0, 1] ≈ -cᵤ atol = atol
        @test dot(R' * [1, 0, 0], aᵤ) ≈ 0 atol = atol

        # XZ 
        R = FrameTransformations.twodir_to_dcm(a, b, :XZ)

        @test R' * [1, 0, 0] ≈ aᵤ atol = atol
        @test R' * [0, 1, 0] ≈ -cᵤ atol = atol
        @test dot(R' * [0, 1, 0], aᵤ) ≈ 0 atol = atol

        # ZX 
        R = FrameTransformations.twodir_to_dcm(a, b, :ZX)

        @test R' * [0, 0, 1] ≈ aᵤ atol = atol
        @test R' * [0, 1, 0] ≈ cᵤ atol = atol
        @test dot(R' * [1, 0, 0], aᵤ) ≈ 0 atol = atol

        # YZ 
        R = FrameTransformations.twodir_to_dcm(a, b, :YZ)

        @test R' * [0, 1, 0] ≈ aᵤ atol = atol
        @test R' * [1, 0, 0] ≈ cᵤ atol = atol
        @test dot(R' * [0, 0, 1], aᵤ) ≈ 0 atol = atol

        # ZY
        R = FrameTransformations.twodir_to_dcm(a, b, :ZY)

        @test R' * [0, 0, 1] ≈ aᵤ atol = atol
        @test R' * [1, 0, 0] ≈ -cᵤ atol = atol
        @test dot(R' * [0, 1, 0], aᵤ) ≈ 0 atol = atol

    end
end;

@testset "Frames" begin

    frames = FrameSystem{4,Float64}()
    add_axes_icrf!(frames)

    eph = EphemerisProvider(path(KERNELS[:DE432]))
    add_point!(frames, :SSB, 0, 1)
    add_point_ephemeris!(frames, eph, :Sun, 10)
    add_point_ephemeris!(frames, eph, :EarthB, 3)
    add_point_ephemeris!(frames, eph, :Earth, 399)
    add_direction_position!(frames, :SunEarthPos, :Sun, :Earth, :ICRF)
    add_direction_velocity!(frames, :SunEarthVel, :Sun, :Earth, :ICRF)
    add_axes_twodir!(frames, :SunEarthRot, 2, :ICRF, :SunEarthPos, :SunEarthVel, :XY)

    @test_throws ArgumentError add_axes_twodir!(frames, :SunEarthRot, 2, :ICRF, :A, :SunEarthVel, :XY)
    @test_throws ArgumentError add_axes_twodir!(frames, :SunEarthRot, 2, :ICRF, :SunEarthPos, :B, :XY)

    for _ = 1:100
        e = rand(0:1e8)
        v = vector3(frames, :Sun, :Earth, :SunEarthRot, e)
        v /= norm(v)
        @test v[1] ≈ 1.0
    end

end;