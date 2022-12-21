import FunctionWrappers: FunctionWrapper

using Basic.Ephemeris: AbstractEphemerisProvider, ephem_position_records

abstract type AbstractFramePoint end 
abstract type AbstractFrameAxes end 

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
    f::FunctionWrapper{Rotation{3, T}, Tuple{T, SVector{3, T}, SVector{3, T}}} 
    δf::FunctionWrapper{Rotation{3, T}, Tuple{T, SVector{6, T}, SVector{6, T}}}
    δ²f::FunctionWrapper{Rotation{3, T}, Tuple{T, SVector{9, T}, SVector{9, T}}}
end

get_node_id(ax::FrameAxesNode) = ax.id

function Base.show(io::IO, ax::FrameAxesNode{T}) where T
    pstr = "FrameAxesNode{$T}(name=$(ax.name), class=$(ax.class), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    println(io, pstr)
end

# -------------------------------------
# POINTS
# -------------------------------------

struct FramePointNode{T} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    axesid::Int      
    parentid::Int
    NAIFId::Int 
    stv::Vector{MVector{9, T}}
    epochs::Vector{T}
    nzo::Vector{Int}
    f::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}} 
    δf::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
    δ²f::FunctionWrapper{Nothing, Tuple{MVector{9, T}, T}}
end 

get_node_id(p::FramePointNode) = p.NAIFId

function Base.show(io::IO, p::FramePointNode{T}) where T
    pstr = "FramePointNode{$T}(name=$(p.name), class=$(p.class), NAIFId=$(p.NAIFId), axes=$(p.axesid)"
    p.parentid == p.NAIFId || (pstr *= ", parent=$(p.parentid)")
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
    return FrameSystem{T}(eph, unique([tids..., cids...]))
end

frames_points(fs::FrameSystem) = fs.points 
frames_axes(fs::FrameSystem) = fs.axes

add_point!(fs::FrameSystem{T}, p::FramePointNode{T}) where {T} = add_vertex!(fs.points, p)
add_axes!(fs::FrameSystem{T}, ax::FrameAxesNode{T}) where {T} = add_vertex!(fs.axes, ax)

@inline has_point(f::FrameSystem, NAIFId::Int) = has_vertex(frames_points(f), NAIFId)
@inline has_axes(f::FrameSystem, axesid::Int) = has_vertex(frames_axes(f), axesid)
@inline ephemeris_points(fs::FrameSystem) = ephemeris_points(fs.prop)

show_points(frame::FrameSystem) = mappedgraph_tree(frames_points(frame))
show_axes(frame::FrameSystem) = mappedgraph_tree(frames_axes(frame))

function mappedgraph_tree(g::MappedNodeGraph)
    s = ""
    s = _mappedgraph_tree!(s, g)
    println(s)
    nothing
end

function _mappedgraph_tree!(s::String, g::MappedNodeGraph)
    if !isempty(g.nodes)
        s *= "\n$(g.nodes[1].name)\n"
        s = _mappedgraph_tree!(s, g, get_node_id(g.nodes[1]), 2, 1)
    end
    s
end

function _mappedgraph_tree!(s::String, g::MappedNodeGraph, pid::Int, idx::Int, del::Int=1)
    @inbounds for i = idx:length(g.nodes)
        if g.nodes[i].parentid == pid 
            s *= "$(" "^(del))├── $(g.nodes[i].name) \n"
            s = _mappedgraph_tree!(s, g, get_node_id(g.nodes[i]), i, del+1)
        end
    end
    s
end

function Base.show(io::IO, fs::FrameSystem{T, E}) where {T, E}
    println(io, "FrameSystem{$T, $E}(")
    println(io, "  eph: $(fs.eph),")

    spoints = ""
    spoints = _mappedgraph_tree!(spoints, frames_points(fs))
    if spoints != ""
        println(io, "  points: $(join(["\t "*si for si in split(spoints, "\n")], "\n"))")
    else
        println(io, "  points: NONE")
    end

    saxes = ""
    saxes = _mappedgraph_tree!(saxes, frames_axes(fs))
    if saxes != ""
        println(io, "  axes: $(join(["\t"*si for si in split(saxes, "\n")], "\n"))")
    else 
        println(io, "  axes: NONE")
    end
    println(io, ")")
end