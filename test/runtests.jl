using Basic
using Test

using ERFA 
using ForwardDiff
using ReferenceFrameRotations
using RemoteFiles
using SPICE 
using StaticArrays

import LinearAlgebra: cross, dot, norm

@RemoteFileSet KERNELS "Generic Spice Kernels" begin
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir=joinpath(@__DIR__, "assets")
    GMS = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc" dir=joinpath(@__DIR__, "assets")
    PCK10 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc" dir=joinpath(@__DIR__, "assets")
    PCK11 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc" dir=joinpath(@__DIR__, "assets")
    DE432 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir=joinpath(@__DIR__, "assets")
    ITRF = @RemoteFile  "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_000101_230424_230129.bpc" dir=joinpath(@__DIR__, "assets")
    PA421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc" dir=joinpath(@__DIR__, "assets")
    PA440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc" dir=joinpath(@__DIR__, "assets")
    FK_DE421 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf" dir=joinpath(@__DIR__, "assets")
    FK_DE440 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf" dir=joinpath(@__DIR__, "assets")
end;
 
download(KERNELS, verbose=true, force=true)

@testset "Basic" verbose=true begin
    @eval begin
        modules = [:Utils, :Tempo, :Orient, :Frames]
        for m in modules
            @testset "$m" verbose=true begin 
                include("$m/$m.jl")
            end         
        end
    end
end;