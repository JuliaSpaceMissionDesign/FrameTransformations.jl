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
    PA421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    PA440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    FK_DE421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    FK_DE440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

@testset "DE421" verbose = false begin

    for kernel in (:LEAP, :PA421, :FK_DE421)
        furnsh(path(KERNELS[kernel]))
    end


    frames = FrameSystem{2,Float64}()
    add_axes_icrf!(frames)

    eph = EphemerisProvider(path(KERNELS[:PA421]))
    @test_throws ArgumentError add_axes_pa421!(frames, eph, :PA421, -1)
    @test_throws ArgumentError add_axes_me421!(frames, :ME421, -1)
    @test_throws ArgumentError add_axes_me421!(frames, :ME421, 31006)

    add_axes_pa421!(frames, eph, :PA421)
    add_axes_me421!(frames, :ME421, :PA421)

    for _ in 1:100
        et = rand(0.0:1e8)

        v = rand(BigFloat, 3)
        v /= norm(v)

        # Test PA421!
        Rb = rotation6(frames, :PA421, :ICRF, et)
        Rs = sxform("MOON_PA", "J2000", et)

        @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
        @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6

        # Test ME421!
        Rb = rotation6(frames, :ICRF, :ME421, et)
        Rs = sxform("J2000", "MOON_ME", et)

        @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
        @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6
    end

end

@testset "DE440" verbose = false begin

    for kernel in (:LEAP, :PA440, :FK_DE440)
        furnsh(path(KERNELS[kernel]))
    end

    frames = FrameSystem{2,Float64}()
    add_axes_icrf!(frames)

    eph = EphemerisProvider(path(KERNELS[:PA440]))
    @test_throws ArgumentError add_axes_pa440!(frames, eph, :PA440, -1)

    add_axes_pa440!(frames, eph, :PA440)

    for _ in 1:100
        et = rand(0.0:1e8)

        v = rand(BigFloat, 3)
        v /= norm(v)

        # Test PA421!
        Rb = rotation6(frames, :PA440, :ICRF, et)
        Rs = sxform("MOON_PA", "J2000", et)

        @test v2as(Rb[1] * v, Rs[1:3, 1:3] * v) ≤ 1e-6
        @test v2as(Rb[2] * v, Rs[4:6, 1:3] * v) ≤ 1e-6
    end

end
