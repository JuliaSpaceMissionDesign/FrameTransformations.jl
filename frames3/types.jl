using StaticArrays, LinearAlgebra 
using ReferenceFrameRotations
import FunctionWrappers: FunctionWrapper

# AXES SECTION 


"""
    ComputableAxesVector 

Store the properties required to retrieve the i-th order components of a 
desired vector. 
"""
struct ComputableAxesVector 
    to::Int 
    from::Int 
    order::Int 
end

ComputableAxesVector() = ComputableAxesVector(0, 0, 0)
function ComputableAxesVector(to::T, from::T,  order::Int) where {T <: Union{Int, <:AbstractPoint}}
    order > 3 && throw(ArgumentError("Order must be <= 3."))
    ComputableAxesVector(get_alias(to), get_alias(from), order)
end

""" 
    ComputableAxesProperties 

Store the properties required to retrieve all the vectors required by  
a computable set of axes. 
"""
struct ComputableAxesProperties 
    v1::ComputableAxesVector 
    v2::ComputableAxesVector 
end

ComputableAxesProperties() = ComputableAxesProperties(ComputableAxesVector(),
                                                      ComputableAxesVector())


"""
    FrameAxes{T} <: AbstractGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `class` -- `Symbol` representing the class of the axes 
- `id` -- axes id (equivalent of NAIFId for axes)
- `parentid` -- id of the parent axes 
- ``
"""
struct FrameAxesNode{T} <: AbstractGraphNode
    name::Symbol            
    class::Symbol          
    id::Int                
    parentid::Int         
    comp::ComputableAxesProperties 
    R::Vector{Rotation{3, T}}
    epochs::Vector{T}
    nzo::Vector{Int} # last updated order
    fun::FunctionWrapper{Rotation{1, T}, Tuple{T, SVector{3, T}, SVector{3, T}}} 
    dfun::FunctionWrapper{Rotation{2, T}, Tuple{T, SVector{6, T}, SVector{6, T}}}
    ddfun::FunctionWrapper{Rotation{3, T}, Tuple{T, SVector{9, T}, SVector{9, T}}}
end

get_node_id(axes::FrameAxesNode) = axes.id

function Base.show(io::IO, ax::FrameAxesNode{T}) where T
    println(io, "FrameAxesNode{$T}")
    println(io, "  name: $(ax.name)")
    println(io, "  class: $(ax.class)")
    println(io, "  id: $(ax.id)")

    ax.parentid == ax.id || println(io, "  parent: $(ax.parentid)")
end


# POINTS SECTION

struct FramePointNode{T} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    axes::Int      
    parent::Int
    NAIFId::Int 
    heph::Bool      # if ephemerides are available for this point
    stv::Vector{MVector{9, T}}
    epochs::Vector{T}
    nzo::Vector{Int}
    fun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}} 
    dfun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
    ddfun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
end 

get_node_id(p::FramePointNode) = p.NAIFId

function Base.show(io::IO, p::FramePointNode{T}) where T 
    println(io, "FramePointNode{$T}")
    println(io, "  name: $(p.name)")
    println(io, "  class: $(p.class)")
    println(io, "  id: $(p.NAIFId)")
    println(io, "  axes: $(p.axes)")

    p.parent == p.NAIFId && println(io, "  parent: $(p.parent)")
end


# FRAMES SYSTEM SECTION

struct FrameSystem{T}
    eph::EphemerisKernels
    ephBodyIDs::Vector{Int}
    pnts_graph::GraphSystem{FramePointNode{T}, SimpleGraph{Int}}
    axes_graph::GraphSystem{FrameAxesNode{T}, SimpleGraph{Int}}
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
                   GraphSystem{FramePointNode{T}}(), 
                   GraphSystem{FrameAxesNode{T}}())
end

get_datatype(::FrameSystem{T}) where {T} = T

points_graph(f::FrameSystem) = f.pnts_graph; 
axes_graph(f::FrameSystem) = f.axes_graph;

has_point(f::FrameSystem, NAIFId::Int) = has_vertex(points_graph(f), NAIFId)
has_axes(f::FrameSystem, axesid::Int) = has_vertex(axes_graph(f), axesid)

add_axes!(f::FrameSystem, a::FrameAxesNode) = add_vertex!(axes_graph(f), a)

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
        if g.nodes[i].parentid == pid 
            println(" "^(3del), g.nodes[i].name)
            _graph_tree(g, get_node_id(g.nodes[i]), i, del+1)
        end 
    end
end
