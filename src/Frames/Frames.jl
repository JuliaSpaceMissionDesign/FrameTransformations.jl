module Frames 

    using Logging
    using ReferenceFrameRotations
    using StaticArrays
    
    import ForwardDiff.derivative
    import FunctionWrappers: FunctionWrapper

    using Basic.Ephemeris: AbstractEphemerisProvider, 
                           NullEphemerisProvider, 
                           ephem_position_records,
                           ephem_timescale,
                           ephem_compute_order!
                           
    using Basic.Tempo: AbstractTimeScale, 
                       BarycentricDynamicalTime, 
                       Epoch, J2000, DJ2000, 
                       CENTURY2DAY, CENTURY2SEC, 
                       DAY2SEC

    using Basic.Utils: format_camelcase
    using Basic.Orient
    using Basic.MappedGraphs

    include("rotation.jl")

    # Frame system and types
    include("types.jl")
    include("axes.jl")
    include("points.jl")
    include("transform.jl")

    # Rotations definitions 
    include("Definitions/default.jl")
    include("Definitions/twovectors.jl")
    include("Definitions/fixedoffset.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/planets.jl")
end