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
    IAU10 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc" dir=joinpath(@__DIR__, "assets")
    IAU11 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc" dir=joinpath(@__DIR__, "assets")
    GMS = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc" dir=joinpath(@__DIR__, "assets")
    SPK = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp" dir=joinpath(@__DIR__, "assets")
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