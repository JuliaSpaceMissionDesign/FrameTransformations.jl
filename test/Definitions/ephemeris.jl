using FrameTransformations
using StaticArrays
using LinearAlgebra
using SPICE 
using RemoteFiles 
using Test

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

@testset "EphemeridesExt" verbose=false begin 
    using Ephemerides
    kclear()

    frames = FrameSystem{4, Float64}()
    add_axes_icrf!(frames)

    eph = EphemerisProvider(path(KERNELS[:DE432]))

    add_point_root!(frames, :SolarSystemBarycenter, 0, 1)

    add_point_ephemeris!(frames, eph, :Sun, 10)
    add_point_ephemeris!(frames, eph, :EarthB, 3)
    add_point_ephemeris!(frames, eph, :Earth, 399)

    furnsh(path(KERNELS[:DE432]))

    for _ in 1:100
        e = rand(0:1e9)
        pv_, lt_ = spkezr("EARTH", e, "J2000", "NONE", "SUN")
        pv = vector6(frames, 10, 399, 1, e)

        @test norm(pv[1:3] - pv_[1:3]) ≈ 0. atol=1e-6
        @test norm(pv[4:end] - pv_[4:end]) ≈ 0. atol=1e-9
    end

    kclear()
end;