module Frames

using FunctionWrappersWrappers: FunctionWrappersWrapper, 
                                FunctionWrappers.FunctionWrapper
using Logging
using PreallocationTools
using ReferenceFrameRotations
using StaticArrays

using SMDGraphs:
    MappedNodeGraph,
    AbstractGraphNode,
    SimpleGraph,
    MappedGraph,
    get_path,
    get_mappedid,
    get_mappednode,
    get_node,
    get_node_id,
    has_vertex,
    add_vertex!,
    add_edge!

using JSMDInterfaces.Ephemeris
using JSMDUtils: format_camelcase, NullEphemerisProvider
using JSMDUtils.Math: D¹, D², D³
using JSMDUtils.Autodiff

using Tempo
using Tempo:
    AbstractTimeScale,
    BarycentricDynamicalTime,
    Epoch,
    J2000,
    DJ2000,
    CENTURY2DAY,
    CENTURY2SEC,
    DAY2SEC,
    j2000

using FrameTransformations.Orient
using FrameTransformations.Orient: AXESID_ICRF
using FrameTransformations.Utils: light_speed, geod2pos
using FrameTransformations.Utils: normalize, δnormalize, δ²normalize, δ³normalize
using FrameTransformations.Utils: cross3, cross6, cross9, cross12
using FrameTransformations.Utils: angle_to_δdcm, angle_to_δ²dcm
using FrameTransformations.Utils: _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm

import LinearAlgebra: dot, norm, matprod, UniformScaling
import StaticArrays: similar_type, Size, MMatrix, SMatrix

import SMDGraphs: get_node_id

include("rotation.jl")

# Frame system and types
include("types.jl")
include("axes.jl")
include("points.jl")
include("lightime.jl")
include("transform.jl")

# Rotations definitions 
include("Definitions/celestial.jl")
include("Definitions/topocentric.jl")
include("Definitions/twovectors.jl")
include("Definitions/ecliptic.jl")
include("Definitions/planets.jl")
include("Definitions/earth.jl")
include("Definitions/moon.jl")

end
