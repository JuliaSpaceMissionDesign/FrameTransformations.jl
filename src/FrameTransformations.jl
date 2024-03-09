module FrameTransformations


# ----
# External dependencies

using PrecompileTools: PrecompileTools
using PreallocationTools
using Reexport
using Logging

using FunctionWrappersWrappers: 
    FunctionWrappersWrapper, 
    FunctionWrappers.FunctionWrapper

using ReferenceFrameRotations
using StaticArrays

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)


import LinearAlgebra: dot, norm, matprod, UniformScaling
import StaticArrays: similar_type, Size, MMatrix, SMatrix
import SMDGraphs: get_node_id


# ----
# JSMD ecosystem dependencies

using JSMDInterfaces
using JSMDInterfaces.Ephemeris
using JSMDInterfaces.Graph: 
    AbstractJSMDGraphNode, 
    add_edge!, 
    add_vertex!, 
    get_path, 
    has_vertex

using JSMDUtils 
using JSMDUtils.Autodiff
using JSMDUtils: format_camelcase, NullEphemerisProvider
using JSMDUtils.Math: 
    D¹, D², D³, 
    unitvec, δunitvec, δ²unitvec, δ³unitvec, 
    cross3, cross6, cross9, cross12, 
    angle_to_δdcm, angle_to_δ²dcm, 
    _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm,
    arcsec2rad

using SMDGraphs
using SMDGraphs:
    MappedNodeGraph, 
    SimpleGraph, 
    MappedGraph,
    get_mappedid, 
    get_mappednode, 
    get_node

using Tempo
using Tempo: AbstractTimeScale

using IERSConventions

import JSMDInterfaces.FilesIO: load
import SMDGraphs: get_node_id

@reexport using Tempo
export load


# Core Routines
# ==============================================

include("Core/index.jl")
include("Core/rotation.jl")
include("Core/twovectors.jl")
include("Core/types.jl")
include("Core/axes.jl")
include("Core/points.jl")
include("Core/lightime.jl")
include("Core/transform.jl")

# Standard axes definitions
# ==============================================

include("Definitions/celestial.jl")
include("Definitions/ecliptic.jl")
include("Definitions/terrestrial.jl")
include("Definitions/topocentric.jl")
include("Definitions/planetary.jl")
include("Definitions/moon.jl")

end
