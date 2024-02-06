module Frames

# External dependencies
# ===========================

using FunctionWrappersWrappers: 
    FunctionWrappersWrapper, 
    FunctionWrappers.FunctionWrapper

using Logging
using PreallocationTools
using ReferenceFrameRotations
using StaticArrays

# JSMD ecosystem dependencies 
# ===========================

using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Graph: 
    AbstractJSMDGraphNode,
    add_edge!, 
    add_vertex!, 
    get_path,
    has_vertex

using JSMDUtils.Autodiff
using JSMDUtils.Math: 
    D¹, 
    D², 
    D³, 
    unitvec, 
    δunitvec, 
    δ²unitvec, 
    δ³unitvec, 
    cross3, 
    cross6, 
    cross9, 
    cross12, 
    angle_to_δdcm, 
    angle_to_δ²dcm, 
    _3angles_to_δdcm, 
    _3angles_to_δ²dcm, 
    _3angles_to_δ³dcm

using JSMDUtils: format_camelcase, NullEphemerisProvider

using SMDGraphs:
    MappedNodeGraph,
    SimpleGraph,
    MappedGraph,
    get_mappedid,
    get_mappednode,
    get_node,
    get_node_id

using Tempo
using Tempo: AbstractTimeScale

using FrameTransformations.Orient
using FrameTransformations.Orient: AXESID_ICRF

# Imports
# ===========================

import LinearAlgebra: dot, norm, matprod, UniformScaling
import StaticArrays: similar_type, Size, MMatrix, SMatrix

import SMDGraphs: get_node_id

# ===========================

include("rotation.jl")

# Frame system and types
include("types.jl")
include("axes.jl")
include("points.jl")
include("lightime.jl")
include("transform.jl")

# Frames definitions
include("Definitions/moon.jl")

end
