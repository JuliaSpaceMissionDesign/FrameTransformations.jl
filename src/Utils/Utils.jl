module Utils

    using JSON3
    using JLD2

    using ForwardDiff: derivative
    using ReferenceFrameRotations: DCM
    using StaticArrays: SMatrix, SA, SVector
    using LinearAlgebra
    
    import YAML as YAMLLib
    
    # IO 
    include("IO/file.jl")
    include("IO/load.jl")
    include("IO/write.jl")
    
    include("constants.jl")

    include("format.jl")
    include("angles.jl")

    include("geodesy.jl")
    
    # Interpolation 
    include("Math/Interpolation/abstract.jl")
    include("Math/Interpolation/akima.jl")
    include("Math/Interpolation/splines.jl")

    # Math
    include("Math/derivatives.jl")
    include("Math/vectors.jl")
    include("Math/rotation.jl")
    
end