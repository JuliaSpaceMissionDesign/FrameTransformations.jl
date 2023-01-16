module Frames 

    using Logging
    using ReferenceFrameRotations
    using StaticArrays

    import FunctionWrappers: FunctionWrapper

    using Basic.Ephemeris 
    using Basic.Orient

    using Basic.Ephemeris: AbstractEphemerisProvider, 
                           NullEphemerisProvider

    using Basic.Tempo: AbstractTimeScale, 
                       BarycentricDynamicalTime, 
                       Epoch, J2000, DJ2000, 
                       CENTURY2DAY, CENTURY2SEC, 
                       DAY2SEC

    using Basic.Utils: format_camelcase
    using Basic.Utils: D¹, D², D³
    using Basic.Utils: angle_to_δdcm, angle_to_δ²dcm
    using Basic.Utils: _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm

    using Basic.MappedGraphs # TODO: remove in future

    include("rotation.jl")

    # Frame system and types
    include("types.jl")
    include("axes.jl")
    include("points.jl")
    include("transform.jl")

    # Rotations definitions 
    include("Definitions/ecliptic.jl")
    include("Definitions/planets.jl")
    include("Definitions/earth.jl")
    include("Definitions/moon.jl")
    include("Definitions/twovectors.jl")
    
end