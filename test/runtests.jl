using FrameTransformations
using Test

using CalcephEphemeris
using ERFA
using ForwardDiff
using ReferenceFrameRotations
using RemoteFiles
using SPICE
using StaticArrays
using Tempo

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Math: interpolate
using JSMDUtils.Math: arcsec2rad, D¹, D², D³
using JSMDUtils: NullEphemerisProvider

using FrameTransformations.Frames
using FrameTransformations.Utils: angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm

import FrameTransformations.Frames:
    FrameAxesFunctions, FramePointFunctions, _get_fixedrot, _empty_stv_update!

import LinearAlgebra: cross, dot, norm

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "assets"
    )
    GMS = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc" dir = joinpath(
        @__DIR__, "assets"
    )
    PCK10 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc" dir = joinpath(
        @__DIR__, "assets"
    )
    PCK11 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc" dir = joinpath(
        @__DIR__, "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "assets"
    )
    ITRF = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_000101_230424_230129.bpc" dir = joinpath(
        @__DIR__, "assets"
    )
    PA421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc" dir = joinpath(
        @__DIR__, "assets"
    )
    PA440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc" dir = joinpath(
        @__DIR__, "assets"
    )
    FK_DE421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf" dir = joinpath(
        @__DIR__, "assets"
    )
    FK_DE440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf" dir = joinpath(
        @__DIR__, "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

@eval begin
    modules = [:Utils, :Orient, :Frames]
    for m in modules
        @testset "$m" verbose = true begin
            include("$m/$m.jl")
        end
    end
end;
