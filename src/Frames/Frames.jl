module Frames 
    using Logging
    import ForwardDiff.derivative

    using Basic.Utils: format_camelcase
    using Basic.MappedGraphs
    using Basic.Tempo: Epoch

    include("rotation.jl")

    # Frame system and types
    include("types.jl")
    include("axes.jl")
    include("points.jl")
    include("transform.jl")

    # Rotations definitions 
    include("Definitions/twovectors.jl")
    include("Definitions/fixedoffset.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/planets.jl")
end