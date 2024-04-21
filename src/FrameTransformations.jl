module FrameTransformations 

using LinearAlgebra
import LinearAlgebra: UniformScaling, matprod

using ForwardDiff
using StaticArrays
using PreallocationTools
using ReferenceFrameRotations
using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

using JSMDUtils
using JSMDUtils.Autodiff
using JSMDInterfaces.Graph: 
    AbstractJSMDGraphNode, add_edge!, add_vertex!, get_path, has_vertex

using SMDGraphs
using SMDGraphs:
    MappedNodeGraph, SimpleGraph, MappedGraph,
    get_mappedid, get_mappednode,  get_node

import SMDGraphs: get_node_id

using Tempo
using Tempo: AbstractTimeScale, Epoch

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

# # Frame system 
# export FrameSystem, order, timescale, points, axes, add_axes!, add_point!,
#        rotation3, rotation6, rotation9, rotation12
# include("Core/node.jl")
# include("Core/graph.jl")
# include("Core/transform.jl")

# # Helpers
# export add_axes_root!, add_axes_projected!, add_axes_fixedoffset! 
# include("Core/axes.jl")

end