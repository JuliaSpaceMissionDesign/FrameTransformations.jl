module Frames 

    using Logging
    using ReferenceFrameRotations
    using StaticArrays

    import LinearAlgebra: dot, norm
    import FunctionWrappers: FunctionWrapper

    using Basic.Ephemeris 
    using Basic.Orient
    using Basic.Orient: AXESID_ICRF

    using Basic.Ephemeris: AbstractEphemerisProvider, 
                           NullEphemerisProvider

    using Basic.Tempo: AbstractTimeScale, 
                       BarycentricDynamicalTime, 
                       Epoch, J2000, DJ2000, 
                       CENTURY2DAY, CENTURY2SEC, 
                       DAY2SEC, j2000

    using Basic.Utils: light_speed
    using Basic.Utils: format_camelcase
    using Basic.Utils: D¹, D², D³
    
    using Basic.Utils: normalize, δnormalize, δ²normalize, δ³normalize
    using Basic.Utils: cross3, cross6, cross9, cross12

    using Basic.Utils: angle_to_δdcm, angle_to_δ²dcm
    using Basic.Utils: _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm

    using Basic.MappedGraphs # TODO: remove in future

    include("rotation.jl")

    # Frame system and types
    include("types.jl")
    include("axes.jl")
    include("points.jl")
    include("lightime.jl")
    include("transform.jl")

    # Rotations definitions 
    include("Definitions/topocentric.jl")
    include("Definitions/twovectors.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/planets.jl")
    include("Definitions/earth.jl")
    include("Definitions/moon.jl")
    
end