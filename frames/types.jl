
using StaticArrays, LinearAlgebra 
using ReferenceFrameRotations

import FunctionWrappers: FunctionWrapper

# AXES SECTION 

struct AstroAxes{T} <: AbstractGraphNode
    name::Symbol            
    class::Int          
    id::Int                
    parent::Int             
    point::Int          # meaningful only for Computable Axes 
    R::Vector{DCM{T}}
    dR::Vector{DCM{T}}
    epochs::Vector{T}
    fun::FunctionWrapper{DCM{T}, Tuple{T, SVector{6, T}}} 
    dfun::FunctionWrapper{DCM{T}, Tuple{T, SVector{6, T}}}
end

get_node_id(axes::AstroAxes) = axes.id

function Base.show(io::IO, ax::AstroAxes{T}) where T
    println(io, "AstroAxes{$T}")
    println(io, "  name: $(ax.name)")
    println(io, "  id: $(ax.id)")
    println(io, "  class: $(ax.class)")

    ax.parent == ax.id && println(io, "  parent: $(ax.parent)")
    ax.class == 33 && println(io, " reference point: $(ax.point)")
end


# POINTS SECTION

struct AstroPoint{T} <: AbstractGraphNode
    name::Symbol
    heph::Bool      # if ephemerides are available for this point
    class::Int
    axes::Int      
    parent::Int
    NAIFId::Int 
    stv::Vector{MVector{6, T}}
    epochs::Vector{T}
    fun!::FunctionWrapper{Nothing, Tuple{MVector{6, T}, T}} 
    dfun!::FunctionWrapper{Nothing, Tuple{MVector{6, T}, T}}
end 

get_node_id(p::AstroPoint) = p.NAIFId

function Base.show(io::IO, p::AstroPoint{T}) where T 
    println(io, "AstroPoint{$T}")
    println(io, "  name: $(p.name)")
    println(io, "  id: $(p.NAIFId)")
    println(io, "  class: $(p.class)")
    println(io, "  axes: $(p.axes)")

    p.parent == p.NAIFId && println(io, "  parent: $(p.parent)")
end


# FRAMES SYSTEM SECTION

struct FrameSystem{T}
    eph::EphemerisKernels
    ephBodyIDs::Vector{Int}
    pnts_graph::GraphSystem{AstroPoint{T}, SimpleGraph{Int}}
    axes_graph::GraphSystem{AstroAxes{T}, SimpleGraph{Int}}
end

FrameSystem{T}() where T = FrameSystem{T}(EphemerisKernels(), Int[])

  
# Mettendo i centri abbiamo a disposizione anche il SSB
# Funziona sempre sta roba? (può capitare che non ci siano dati per ricostruire 
# # quel centro?)
function FrameSystem{T}(eph::EphemerisKernels) where T
    precords = get_position_records(eph)
    tids = map(x->x.target, precords)
    cids = map(x->x.center, precords)
    FrameSystem{T}(eph, unique([tids..., cids...]));
end

function FrameSystem{T}(eph::EphemerisKernels, ids::Vector{Int}) where T
    FrameSystem{T}(eph, ids, 
                   GraphSystem{AstroPoint{T}}(), 
                   GraphSystem{AstroAxes{T}}())
end

get_datatype(::FrameSystem{T}) where {T} = T
points_graph(f::FrameSystem) = f.pnts_graph; 
axes_graph(f::FrameSystem) = f.axes_graph;

available_ephemeris_bodies(frame::FrameSystem) = frame.ephBodyIDs

@inline get_axes(frame::FrameSystem) = _graph_tree(axes_graph(frame))
@inline get_points(frame::FrameSystem) = _graph_tree(points_graph(frame))

function _graph_tree(g::GraphSystem)
    if !isempty(g.nodes)
        println("\n", g.nodes[1].name)
        _graph_tree(g, get_node_id(g.nodes[1]), 2, 1)
    end
end

function _graph_tree(g::GraphSystem, pid::Int, idx::Int, del::Int=1)
    @inbounds for i = idx:length(g.nodes)
        if g.nodes[i].parent == pid 
            println(" "^(3del), g.nodes[i].name)
            _graph_tree(g, get_node_id(g.nodes[i]), i, del+1)
        end 
    end
end
