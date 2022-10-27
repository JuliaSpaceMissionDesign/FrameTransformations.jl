module Frames 

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations
    using Basic.Tempo
    using Basic.Orient

    export Rotation

    include("abstract.jl")
    include("graph.jl")
    include("rotation.jl")
    
    include("definitions/types.jl")
    include("definitions/singletons.jl")
    include("definitions/ecliptic.jl")

end