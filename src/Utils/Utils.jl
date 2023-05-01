module Utils

using ForwardDiff
using ForwardDiff: derivative
using LinearAlgebra
using ReferenceFrameRotations: DCM
using StaticArrays: SMatrix, SA, SVector

using InterfacesUtils
import InterfacesUtils.IO: load
using InterfacesUtils.IO: @filetype, AbstractFile

using PrecompileTools: PrecompileTools

# IO 
include("IO/tpc.jl")

include("constants.jl")
include("geodesy.jl")

# Math
include("Math/vectors.jl")
include("Math/rotation.jl")

# Precompilation routines 
PrecompileTools.@setup_workload begin
    x3 = rand(3)
    x6 = rand(6)
    x9 = rand(9)
    x12 = rand(12)

    x3s = SA[rand(3)...]
    x6s = SA[rand(6)...]
    x9s = SA[rand(9)...]
    x12s = SA[rand(12)...]

    vfcns = (
        (normalize, cross3),
        (δnormalize, cross6),
        (δ²normalize, cross9),
        (δ³normalize, cross12),
    )

    vects = ((x3, x3s), (x6, x6s), (x9, x9s), (x12, x12s))

    PrecompileTools.@compile_workload begin

        # Precompile Geodesy routines
        for x in (x3, x3s)
            geoc2pos(x)
            pos2geoc(x)
            pos2geod(x, 231.0, 0.0123)
        end

        for x in (x6, x6s)
            geoc2pv(x6)
            pv2geoc(x6s)
        end

        geod2pos(1000.3, 1.0, 1.0, 231.0, 0.0123)

        # Precompile vector routines
        for (vfcn, vx) in zip(vfcns, vects)
            for x in vx
                vfcn[1](x)
                vfcn[2](x, x)
            end
        end

        # Precompile rotation routines 
        skew(x3)
        skew(x3s)

        angle_to_δdcm(x3, :Z)
        angle_to_δdcm(x3, x3, :ZX)
        angle_to_δdcm(x3, x3, x3, :ZXZ)

        angle_to_δ²dcm(x3, :Z)
        angle_to_δ²dcm(x3, x3, :ZX)
        angle_to_δ²dcm(x3, x3, x3, :ZXZ)

        angle_to_δ³dcm(x6, :Z)
        angle_to_δ³dcm(x6, x6, :ZX)
        angle_to_δ³dcm(x6, x6, x6, :ZXZ)

        _3angles_to_δdcm(x6, :ZXZ)
        _3angles_to_δdcm(x6s, :ZXZ)

        _3angles_to_δ²dcm(x9, :ZXZ)
        _3angles_to_δ²dcm(x9s, :ZXZ)

        _3angles_to_δ³dcm(x12, :ZXZ)
        _3angles_to_δ³dcm(x12s, :ZXZ)
    end
end

end
