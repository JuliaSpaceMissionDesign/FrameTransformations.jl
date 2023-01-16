module Utils

    using JSON3
    using JLD2

    using ForwardDiff: derivative
    using ReferenceFrameRotations: DCM
    using StaticArrays: SMatrix 
    
    import YAML as YAMLLib

    include("IO/file.jl")
    include("IO/load.jl")
    include("IO/write.jl")
    include("format.jl")
    include("angles.jl")
    
    # Interpolation 
    include("Math/derivatives.jl")
    include("Math/akima.jl")
    include("Math/rotation.jl")
    
end