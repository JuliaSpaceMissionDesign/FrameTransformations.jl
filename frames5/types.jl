import FunctionWrappers: FunctionWrapper

using Basic.Ephemeris: AbstractEphemerisProvider, NullEphemerisProvider, 
                       ephem_position_records

using Basic.Tempo: AbstractTimeScale, BarycentricDynamicalTime, Epoch

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
    ComputableAxesVector(point_alias(to), point_alias(from), order)
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


# Frame Axes Function signatures 
_FAxesFunIn{N, T} = Tuple{T, SVector{N, T}, SVector{N, T}}
_FAxesFunSig{O, T, N} = FunctionWrapper{Rotation{O, T}, _FAxesFunIn{N, T}}

# Container to store frame axes update functions 
struct FrameAxesFunctions{T, O, N}
    fun::NTuple{O, _FAxesFunSig{O, T, N}}
end

Base.getindex(af::FrameAxesFunctions, i) = af.fun[i]

# Default rotation function for axes that do not require updates
_get_fixedrot(::T, ::SVector{N, T}, ::SVector{N, T}) where {N, T} = Rotation{N/3}(T(1)I)       

# Constructors for FrameAxesFunctions 
@generated function FrameAxesFunctions{T}(funs::Function...) where T
    O = length(funs)
    expr = :(tuple($([Expr(:ref, :funs, i) for i in 1:O]...)))

    return quote 
        @inbounds FrameAxesFunctions{T, $O, 3*$O}($(expr))
    end
end

# Constructor to filter out some of the specified functions!
@generated function FrameAxesFunctions{T, O}(funs::Function...) where {T, O}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = :(tuple($([Expr(:ref, :funs, i) for i in 1:O]...)))
    return quote 
        @inbounds FrameAxesFunctions{T, O, 3*O}($(expr))
    end
end

# Default constructors for dummy axes function updates 
@generated function FrameAxesFunctions{T, O}() where {T, O}
    expr = :(tuple($([_get_fixedrot for i in 1:O]...)))
    return quote 
        FrameAxesFunctions{T, O, 3*O}($(expr))
    end
end

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
struct FrameAxesNode{O, T, N} <: AbstractGraphNode
    name::Symbol            
    class::Symbol          
    id::Int                
    parentid::Int         
    comp::ComputableAxesProperties 
    R::Vector{Rotation{O, T}}
    epochs::Vector{T}
    nzo::Vector{Int} # last updated order
    f::FrameAxesFunctions{T, O, N}
end

get_node_id(ax::FrameAxesNode) = ax.id

function Base.show(io::IO, ax::FrameAxesNode{O, T}) where {O, T}
    pstr = "FrameAxesNode{$O, $T}(name=$(ax.name), class=$(ax.class), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    println(io, pstr)
end


# -------------------------------------
# POINTS
# -------------------------------------

# Frame Point Function signatures 
_FPointFunIn{N, T} = Tuple{MVector{N, T}, T}
_FPointFunSig{T, N} = FunctionWrapper{Nothing, _FPointFunIn{N, T}}

# Container to store frame point update functions 
struct FramePointFunctions{T, O, N}
    fun::NTuple{O, _FPointFunSig{T, N}}
end

Base.getindex(pf::FramePointFunctions, i) = pf.fun[i]

# Default state-vector update function for points that do not require updates
_empty_stv_update!(::AbstractVector{T}, ::T) where {T} = nothing

# Constructors for FramePointFunctions 
@generated function FramePointFunctions{T}(funs::Function...) where T
    O = length(funs)
    expr = :(tuple($([(Expr(:ref, :funs, i)) for i in 1:O]...)))

    return quote 
        @inbounds FramePointFunctions{T, $O, 3*$O}($(expr))
    end
end

# Constructor to filter out some of the specified functions!
@generated function FramePointFunctions{T, O}(funs::Function...) where {T, O}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = :(tuple($([Expr(:ref, :funs, i) for i in 1:O]...)))
    return quote 
        @inbounds FramePointFunctions{T, O, 3*O}($(expr))
    end
end

# Default constructors for dummy point function updates 
@generated function FramePointFunctions{T, O}() where {T, O}
    expr = :(tuple($([_empty_stv_update! for i in 1:O]...)))
    return quote 
        FramePointFunctions{T, O, 3*O}($(expr))
    end
end


struct FramePointNode{O, T, N} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    axesid::Int      
    parentid::Int
    NAIFId::Int 
    stv::Vector{MVector{N, T}}
    epochs::Vector{T}
    nzo::Vector{Int}
    f::FramePointFunctions{T, O, N}
end 

get_node_id(p::FramePointNode) = p.NAIFId

function Base.show(io::IO, p::FramePointNode{O, T}) where {O, T}
    pstr = "FramePointNode{$O, $T}(name=$(p.name), class=$(p.class), NAIFId=$(p.NAIFId), axes=$(p.axesid)"
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
FrameSystemProperties() = FrameSystemProperties(Int64[])

@inline ephemeris_points(fsp::FrameSystemProperties) = fsp.ebid

struct FrameSystem{O, T <: Number, S <: AbstractTimeScale, E <: AbstractEphemerisProvider, N}
    eph::E
    prop::FrameSystemProperties{T}
    points::MappedNodeGraph{FramePointNode{O, T, N}, SimpleGraph{Int}}
    axes::MappedNodeGraph{FrameAxesNode{O, T, N}, SimpleGraph{Int}}
end

function FrameSystem{O, T, S}(eph::E, points::Vector{Int}) where {O, T <: Number, 
            S <: AbstractTimeScale, E <:AbstractEphemerisProvider} 
    
    if O < 1 || O > 4
        throw(ArgumentError("FrameSystem order must be between 1 and 4."))
    end

    return FrameSystem{O, T, S, E, 3*O}(
        eph, FrameSystemProperties{T}(points),
        MappedGraph(FramePointNode{O, T, 3O}),
        MappedGraph(FrameAxesNode{O, T, 3O})
    )
end

function FrameSystem{O, T}(eph::E) where {O, T, E}
    prec = ephem_position_records(eph)
    tids = map(x->x.target, prec)
    cids = map(x->x.center, prec)
    S = typeof(ephem_timescale(eph))
    return FrameSystem{O, T, S}(eph, unique([tids..., cids...]))
end

FrameSystem{O, T, S}() where {O, T, S} = FrameSystem{O, T, S}(NullEphemerisProvider(), Int64[])
FrameSystem{O, T}() where {O, T <: Number} = FrameSystem{O, T, BarycentricDynamicalTime}()

frames_order(::FrameSystem{O}) where O = O 
frames_timescale(::FrameSystem{O, T, S}) where {O, T, S} = S

frames_points(fs::FrameSystem) = fs.points 
frames_axes(fs::FrameSystem) = fs.axes

add_point!(fs::FrameSystem{O, T}, p::FramePointNode{O, T}) where {O, T} = add_vertex!(fs.points, p)
add_axes!(fs::FrameSystem{O, T}, ax::FrameAxesNode{O, T}) where {O, T} = add_vertex!(fs.axes, ax)

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

function Base.summary(io::IO, ::FrameSystem{O, T, S, E}) where {O, T, S, E}
    println(io, "FrameSystem{$O, $T, $S, $E}")
end

function Base.show(io::IO, fs::FrameSystem{O, T, S, E}) where {O, T, S, E}
    println(io, "FrameSystem{$O, $T, $S, $E}(")
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