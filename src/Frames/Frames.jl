module Frames 

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations
    
    using Basic.Tempo
    using Basic.Orient

    include("abstract.jl")
    include("graph.jl")
    include("transform.jl")

    export Rotation

    include("Definitions/types.jl")
    include("Definitions/default.jl")
    include("Definitions/planets.jl")
end