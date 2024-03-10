using FrameTransformations

using Test
using SPICE
using StaticArrays
using ReferenceFrameRotations
using RemoteFiles
using PreallocationTools

import LinearAlgebra: cross, dot, norm

using Ephemerides
using IERSConventions
using Tempo

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Math: interpolate

using JSMDUtils.Math
using JSMDUtils.Math: angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm, arcsec2rad, D¹, D², D³
using JSMDUtils: NullEphemerisProvider

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Math: interpolate

@RemoteFileSet KERNELS "Spice Kernels Set" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    PA440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    PA421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    FK_DE440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    PCK10 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

download(KERNELS; verbose=true, force=false)

@eval begin 
    @testset "Core" verbose=true begin
        include("types.jl")
        include("axes.jl")
        include("points.jl")
        include("rotation.jl")
        include("lightime.jl")
        include("twovectors.jl")
    end
end;