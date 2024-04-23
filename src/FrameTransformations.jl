module FrameTransformations 

using LinearAlgebra
import LinearAlgebra: UniformScaling, matprod

using ForwardDiff
using StaticArrays
using ReferenceFrameRotations
using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

using JSMDUtils
using JSMDUtils.Autodiff
using JSMDUtils.Math: D¹, D², D³ 
using JSMDInterfaces.Graph: AbstractJSMDGraphNode, add_edge!, add_vertex!, get_path, has_vertex
using JSMDInterfaces.Ephemeris: AbstractEphemerisProvider, ephem_position_records, ephem_available_points
using JSMDInterfaces.Interface: @interface

using SMDGraphs: MappedNodeGraph, SimpleGraph, MappedGraph,
                 get_mappedid, get_mappednode,  get_node, get_path

import SMDGraphs: get_node_id

using Tempo: AbstractTimeScale, Epoch, j2000s, BarycentricDynamicalTime

# Core Routines
# ==============================================

# Low-level types and aliases
export Rotation

include("Core/rotation.jl")
include("Core/ad.jl")

# Frame system (graph & nodes)
export FrameSystem, get_order, get_timescale, get_points, get_axes,
       has_axes, has_point, add_axes!, add_point!,
       rotation3, rotation6, rotation9, rotation12, vector3, vector6, vector9, vector12

include("Core/nodes.jl")
include("Core/graph.jl")
include("Core/transform.jl")

# Helper functions 
export add_axes_root!, add_axes_inertial!, add_axes_rotating!, add_axes_fixedoffset!,
       add_point_root!, add_point_dynamical!, add_point_fixedoffset!

include("Core/axes.jl")
include("Core/points.jl")

# Standard axes definitions
# ==============================================

export add_point_ephemeris!

include("Definitions/ephemeris.jl")

end