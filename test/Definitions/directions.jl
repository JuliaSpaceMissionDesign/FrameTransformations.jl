using FrameTransformations
using Ephemerides
using RemoteFiles
using Test
using LinearAlgebra

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

frames = FrameSystem{4,Float64}()
add_axes_icrf!(frames)

eph = EphemerisProvider(path(KERNELS[:DE432]))

add_point!(frames, :SSB, 0, 1)
add_point_alias!(frames, :SSB, :SolarSystemB)

add_point_ephemeris!(frames, eph, :Sun, 10)
add_point_ephemeris!(frames, eph, :EarthB, 3)
add_point_ephemeris!(frames, eph, :Earth, 399)

@testset "Position" verbose = false begin
    add_direction_position!(frames, :SunEarth, :Sun, :Earth, :ICRF)

    for _ = 1:100
        e = rand(0:1e8)

        ref = vector3(frames, :Sun, :Earth, :ICRF, e)
        ref /= norm(ref)

        val = direction3(frames, :SunEarth, :ICRF, e)
        @test dot(val, ref) ≈ 1.0
    end
end

@testset "Velocity" verbose = false begin
    add_direction_velocity!(frames, :SunEarthVel, :Sun, :Earth, :ICRF)

    for _ = 1:100
        e = rand(0:1e8)

        @views ref = vector6(frames, :Sun, :Earth, :ICRF, e)[4:end]
        ref /= norm(ref)

        val = direction3(frames, :SunEarthVel, :ICRF, e)
        @test dot(val, ref) ≈ 1.0
    end
end

@testset "Orthogonal" verbose = false begin
    add_direction_orthogonal!(frames, :SunEarthMom, :SunEarth, :SunEarthVel, :ICRF)

    for _ = 1:100
        e = rand(0:1e8)

        p, v = Translation(vector6(frames, :Sun, :Earth, :ICRF, e)...)
        p /= norm(p)
        v /= norm(v)
        ref = cross(p, v)
        ref /= norm(ref)

        val = direction3(frames, :SunEarthMom, :ICRF, e)
        @test dot(val, ref) ≈ 1.0
    end
end

@testset "Fixed" verbose = false begin
    add_direction_fixed!(frames, :Fix, :ICRF, [1.0, 0.0, 0.0])

    for _ = 1:100
        e = rand(0:1e8)

        val = direction3(frames, :Fix, :ICRF, e)
        @test val[1] ≈ 1.0
    end
end
