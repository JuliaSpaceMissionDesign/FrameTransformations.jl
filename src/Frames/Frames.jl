module Frames

using Logging
using ReferenceFrameRotations
using StaticArrays

import LinearAlgebra: dot, norm
import FunctionWrappers: FunctionWrapper

using Basic.Ephemeris
using Basic.Orient
using Basic.Tempo
using Basic.Orient: AXESID_ICRF

using Basic.Ephemeris: AbstractEphemerisProvider, NullEphemerisProvider

using Basic.Tempo:
    AbstractTimeScale,
    BarycentricDynamicalTime,
    Epoch,
    J2000,
    DJ2000,
    CENTURY2DAY,
    CENTURY2SEC,
    DAY2SEC,
    j2000

using Basic.Utils: light_speed, geod2pos
using Basic.Utils: normalize, δnormalize, δ²normalize, δ³normalize
using Basic.Utils: cross3, cross6, cross9, cross12
using Basic.Utils: angle_to_δdcm, angle_to_δ²dcm
using Basic.Utils: _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm

using InterfacesUtils.Utils: format_camelcase
using InterfacesUtils.Math: D¹, D², D³

using MultiGraphs:
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

import MultiGraphs: get_node_id

import PrecompileTools

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

# Precompilation routines 
PrecompileTools.@setup_workload begin 

    x12 = rand(12)

    x3s = SA[rand(3)...]
    x6s = SA[rand(6)...]
    x9s = SA[rand(9)...]
    x12s = SA[rand(12)...]

    PrecompileTools.@compile_workload begin 

        # Precompile twovectors routines 
        twovectors_to_dcm(x3s, x3s, :XZ)
        twovectors_to_δdcm(x6s, x6s, :XZ)
        twovectors_to_δ²dcm(x9s, x9s, :XZ)
        twovectors_to_δ³dcm(x12s, x12s, :XZ)

        twovectors_to_dcm(x12, x12, :XZ)
        twovectors_to_δdcm(x12, x12, :XZ)
        twovectors_to_δ²dcm(x12, x12, :XZ)
        twovectors_to_δ³dcm(x12, x12, :XZ)

        _two_vectors_to_rot6(x6s, x6s, :XZ)
        _two_vectors_to_rot9(x9s, x9s, :XZ)
        _two_vectors_to_rot12(x12s, x12s, :XZ)

        _two_vectors_to_rot6(x12, x12, :XZ)
        _two_vectors_to_rot9(x12, x12, :XZ)
        _two_vectors_to_rot12(x12, x12, :XZ)

    end
end


end
