using FrameTransformations
using StaticArrays
using LinearAlgebra
using SPICE
using RemoteFiles
using Test

using Tempo
using JSMDInterfaces.Ephemeris

using Ephemerides

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

@testset "Interface" verbose = false begin
    frames = FrameSystem{2,Float64}()
    add_axes_icrf!(frames)

    @test_throws Exception add_point_ephemeris!(frames, eph, :Sun, 10)
    @test_throws Exception add_axes_ephemeris!(frames, eph, :Sun, 10, :XYZ, 1)

    eph = EphemerisProvider(path(KERNELS[:DE432]))

    add_point!(frames, :SolarSystemBarycenter, 0, 1)

    add_point_ephemeris!(frames, eph, :Sun, 10)
    add_point_ephemeris!(frames, eph, :EarthB, 3)
    add_point_ephemeris!(frames, eph, :Earth, 399)

    frames2 = FrameSystem{2,Float64}()
    add_axes_icrf!(frames2)

    book = Dict(
        0 => :SolarSystemBarycenter, 10 => :Sun, 3 => :EarthB, 399 => :Earth
    )
    @test_nowarn add_point_ephemeris!(frames2, eph, book)
    @test points_alias(frames) == points_alias(frames2)
    @test_throws Exception FrameTransformations.check_point_ephemeris(frames, eph, -1)

    frames = FrameSystem{2,Float64}()
    add_axes_icrf!(frames)
    @test_throws Exception FrameTransformations.check_point_ephemeris(frames, eph, 10)

    frames = FrameSystem{2,Float64}()
    @test_throws Exception FrameTransformations.check_point_ephemeris(frames, eph, 0)

    frames = FrameSystem{2,Float64}()
    @test_throws Exception FrameTransformations.check_axes_ephemeris(frames, eph, -1)
end;


@testset "EphemeridesExt" verbose = false begin
    kclear()

    frames = FrameSystem{4,Float64}()
    add_axes_icrf!(frames)

    eph = EphemerisProvider(path(KERNELS[:DE432]))

    add_point!(frames, :SolarSystemBarycenter, 0, 1)

    add_point_ephemeris!(frames, eph, :Sun, 10)
    add_point_ephemeris!(frames, eph, :EarthB, 3)
    add_point_ephemeris!(frames, eph, :Earth, 399)

    furnsh(path(KERNELS[:DE432]))

    for _ in 1:100
        e = rand(0:1e9)
        pv_, lt_ = spkezr("EARTH", e, "J2000", "NONE", "SUN")
        pv = vector6(frames, 10, 399, 1, e)

        @test norm(pv[1:3] - pv_[1:3]) ≈ 0.0 atol = 1e-6
        @test norm(pv[4:end] - pv_[4:end]) ≈ 0.0 atol = 1e-9
    end

    @test_nowarn vector3(frames, 10, 399, 1, 1.0)
    @test_nowarn vector9(frames, 10, 399, 1, 1.0)
    @test_nowarn vector12(frames, 10, 399, 1, 1.0)

    kclear()
end;
