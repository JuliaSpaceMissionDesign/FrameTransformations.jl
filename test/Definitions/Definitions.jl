using FrameTransformations

using Test
using LinearAlgebra
using SPICE
using StaticArrays
using ReferenceFrameRotations
using Ephemerides
using RemoteFiles
using PreallocationTools
using Tempo

import LinearAlgebra: cross, dot, norm

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Math: interpolate

using JSMDUtils.Math
using JSMDUtils.Math: angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm, arcsec2rad, D¹, D², D³
using JSMDUtils: NullEphemerisProvider

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Math: interpolate

using IERSConventions

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
    FK_DE421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf" dir = joinpath(
        @__DIR__, "..", "assets"
    )
    FK_DE440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf" dir = joinpath(
        @__DIR__, "..", "assets"
    )
end;

EOP_DATA_FILE = @RemoteFile "https://datacenter.iers.org/data/csv/finals2000A.data.csv" dir = joinpath(
    @__DIR__, "..", "assets"
);

download(KERNELS; verbose=true, force=false)
download(EOP_DATA_FILE; verbose=true, force=false)

@info "Prepare EOP data"
let
    eopfile = joinpath(@__DIR__, "..", "assets", "iau2000a")
    eop_generate_from_csv(iers2010b, path(EOP_DATA_FILE), eopfile)
    eop_load_data!(iers2010b, eopfile*".eop.dat")
end;

@eval begin 
    @testset "Definitions" verbose=true begin
        include("celestial.jl")
        include("ecliptic.jl")
        include("terrestrial.jl")
        include("moon.jl")
    end
end;