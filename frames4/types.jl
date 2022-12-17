import FunctionWrappers: FunctionWrapper

using Basic.Ephemeris: AbstractEphemerisProvider, ephem_position_records

abstract type AbstractFramePoint end # TODO: move to proper location
abstract type AbstractFrameAxes end # TODO: move to proper location

# -------------------------------------
# AXES
# -------------------------------------

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
function ComputableAxesVector(to::T, from::T, order::Int) where {T <: Union{Int, <:AbstractFramePoint}}
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
    FrameAxesNode{T} <: AbstractGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `class` -- `Symbol` representing the class of the axes 
- `id` -- axes id (equivalent of NAIFId for axes)
- `parentid` -- id of the parent axes 
- `comp` -- properties for computable axes 
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

get_node_id(ax::FrameAxesNode) = ax.id

function Base.show(io::IO, ax::FrameAxesNode{T}) where T
    pstr = "FrameAxesNode{$T}(name=$(ax.name), class=$(ax.class), id=$(ax.id)"
    if !(ax.parentid == ax.id) 
        pstr *= ", parent=$(ax.parentid)"
    end
    pstr *= ")"
    println(io, pstr)
end

# -------------------------------------
# POINTS
# -------------------------------------

struct FramePointNode{T} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    axes::Int      
    parentid::Int
    NAIFId::Int 
    iseph::Bool      # TODO: se possibile, togliere
    stv::Vector{MVector{9, T}}
    epochs::Vector{T}
    nzo::Vector{Int}
    fun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}} 
    dfun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
    ddfun!::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
end 

get_node_id(p::FramePointNode) = p.NAIFId

function Base.show(io::IO, p::FramePointNode{T}) where T
    pstr = "FramePointNode{$T}(name=$(p.name), class=$(p.class), id=$(p.id), axes=$(p.axes)"
    if !(p.parentid == p.id)
        pstr *= ", parent=$(p.parentid)"
    end
    pstr *= ")"
    println(io, pstr)
end

# -------------------------------------
# FRAMES
# -------------------------------------

struct FrameSystemProperties{T}
    ebid::Vector{Int}  # ephemeris body ids
end
FrameSystemProperties() = FrameSystemProperties([])

@inline ephemeris_points(fsp::FrameSystemProperties) = fsp.ebid

struct FrameSystem{T, E<:AbstractEphemerisProvider}
    eph::E
    prop::FrameSystemProperties{T}
    points::MappedNodeGraph{FramePointNode{T}, SimpleGraph{Int}}
    axes::MappedNodeGraph{FrameAxesNode{T}, SimpleGraph{Int}}
end

function FrameSystem{T}(eph::E, points::Vector{Int}) where {T, E<:AbstractEphemerisProvider} 
    return FrameSystem{T, E}(
        eph, FrameSystemProperties{T}(points),
        MappedGraph(FramePointNode{T}),
        MappedGraph(FrameAxesNode{T})
    )
end

function FrameSystem{T}(eph::E) where {T, E}
    prec = ephem_position_records(eph)
    tids = map(x->x.target, prec)
    cids = map(x->x.center, prec)
    return FrameSystem{T, E}(eph, unique([tids..., cids...]))
end

frames_points(fs::FrameSystem) = fs.points 
frames_axes(fs::FrameSystem) = fs.axes

add_point!(fs::FrameSystem{T}, p::FramePointNode{T}) where {T} = add_vertex!(fs.points, p)
add_axes!(fs::FrameSystem{T}, ax::FrameAxesNode{T}) where {T} = add_vertex!(fs.axes, ax)

@inline has_point(f::FrameSystem, NAIFId::Int) = has_vertex(frames_points(f), NAIFId)
@inline has_axes(f::FrameSystem, axesid::Int) = has_vertex(frames_axes(f), axesid)
@inline ephemeris_points(fs::FrameSystem) = ephemeris_points(fs.prop)

show_points(frame::FrameSystem) = _graph_tree(frames_points(frame))
show_axes(frame::FrameSystem) = _graph_tree(frames_axes(frame))

# TODO: da fare piu' carino
function _graph_tree(g::MappedNodeGraph)
    if !isempty(g.nodes)
        println("\n", g.nodes[1].name)
        _graph_tree(g, get_node_id(g.nodes[1]), 2, 1)
    end
end

# TODO: da fare piu' carino
function _graph_tree(g::MappedNodeGraph, pid::Int, idx::Int, del::Int=1)
    @inbounds for i = idx:length(g.nodes)
        if g.nodes[i].parentid == pid 
            println(" "^(3del), g.nodes[i].name)
            _graph_tree(g, get_node_id(g.nodes[i]), i, del+1)
        end 
    end
end