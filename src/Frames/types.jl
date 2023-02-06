export FrameSystem, ComputableAxesVector,
       frames_order, frames_timescale,
       frames_axes, frames_points, 
       show_axes, show_points,
       add_axes!, add_point!,
       has_axes, has_point, 
       ephemeris_points


# -------------------------------------
# ABSTRACT 
# -------------------------------------
"""
    AbstractFramePoint

Abstract type for all reference frames points.
"""
abstract type AbstractFramePoint end 

"""
    AbstractFrameAxes

Abstract type for all reference frames axes.
"""
abstract type AbstractFrameAxes end 


# -------------------------------------
# AXES
# -------------------------------------

"""
    ComputableAxesVector(from, to, order::Int)

Store the properties required to retrieve the i-th order components of a 
desired vector. Arguments `from` and `to` are the NAIFIDs or `AbstractFramePoint` instances 
that define the observer and target points.

Only orders between 1 and 3 are supported.

### Example 
```jldoctest
julia> @point SSB 0 SolarSystemBarycenter

julia> @point Sun 10 

julia> ComputableAxesVector(SSB, Sun, 1)
ComputableAxesVector(0, 10, 1)

julia> ComputableAxesVector(0, 10, 1)
ComputableAxesVector(0, 10, 1)
```
"""
struct ComputableAxesVector 
    from::Int 
    to::Int 
    order::Int 

    function ComputableAxesVector(from::Int, to::Int, order::Int)
        if order < 1 || order > 3
            throw(ArgumentError("order must be between 1 and 3."))
        end

        if to == from 
            throw(ArgumentError("vector origin and target must be different."))
        end

        new(from, to, order)
    end

end
 
ComputableAxesVector() = ComputableAxesVector(1, 2, 1)
function ComputableAxesVector(from, to, order::Int)
    ComputableAxesVector(point_alias(from), point_alias(to), order)
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
    FrameAxesNode{O, T, N} <: AbstractGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `class` -- `Symbol` representing the class of the axes 
- `id` -- axes ID (equivalent of NAIFId for axes)
- `parentid` -- ID of the parent axes 
- `comp` -- properties for computable axes 
- `R` -- vector storing rotation matrices 
- `epochs` -- vector storing the epochs associated to `R`
- `nzo` -- last order at which `R` has been computed 
- `f` -- `FrameAxesFunctions` container 
- `angles` -- vector storing the libration angles retrived from ephemerides
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
    angles::Vector{MVector{N, T}}
end

MappedGraphs.get_node_id(ax::FrameAxesNode) = ax.id

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

"""
    FramePointNode{O, T, N} <: AbstractGraphNode

Define a frame system point.

### Fields
- `name` -- point name 
- `class` -- `Symbol` representing the class of the point 
- `axesid` -- ID of the axes in which the point coordinates are expressed 
- `parentid` -- NAIF ID of the parent point 
- `NAIFId` -- NAIF ID of the point
- `stv` -- vector storing the point state vectors
- `epochs` -- vector storing the epochs associated to `stv`
- `nzo` -- last order at which `stv` has been computed 
- `f` -- `FramePointFunctions` container 
"""
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

MappedGraphs.get_node_id(p::FramePointNode) = p.NAIFId

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
    eaid::Vector{Int}  # ephemeris axes ids 
end
FrameSystemProperties() = FrameSystemProperties(Int64[], Int64[])

@inline ephemeris_points(fsp::FrameSystemProperties) = fsp.ebid
@inline ephemeris_axes(fsp::FrameSystemProperties) = fsp.eaid

# TODO: add documentation!
"""
    FrameSystem
"""
struct FrameSystem{O, T <: Number, S <: AbstractTimeScale, E <: AbstractEphemerisProvider, N}
    eph::E
    prop::FrameSystemProperties{T}
    points::MappedNodeGraph{FramePointNode{O, T, N}, SimpleGraph{Int}}
    axes::MappedNodeGraph{FrameAxesNode{O, T, N}, SimpleGraph{Int}}
end

function FrameSystem{O, T, S}(eph::E, points::Vector{Int}, axes::Vector{Int}) where {O, 
            T <: Number, S <: AbstractTimeScale, E <:AbstractEphemerisProvider} 
    
    if O < 1 || O > 4
        throw(ArgumentError("FrameSystem order must be between 1 and 4."))
    end

    return FrameSystem{O, T, S, E, 3*O}(
        eph, FrameSystemProperties{T}(points, axes),
        MappedGraph(FramePointNode{O, T, 3O}),
        MappedGraph(FrameAxesNode{O, T, 3O})
    )
end

function FrameSystem{O, T}(eph::E) where {O, T, E}

    points = ephem_available_points(eph)
    axes = ephem_available_axes(eph) 

    S = typeof(ephem_timescale(eph))
    return FrameSystem{O, T, S}(eph, points, axes)
end

FrameSystem{O, T, S}() where {O, T, S} = FrameSystem{O, T, S}(NullEphemerisProvider(), Int64[], Int64[])
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
@inline ephemeris_axes(fs::FrameSystem) = ephemeris_axes(fs.prop)

show_points(frame::FrameSystem) = mappedgraph_tree(frames_points(frame))
show_axes(frame::FrameSystem) = mappedgraph_tree(frames_axes(frame))


# -------------------------------------
# UTILS
# -------------------------------------

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
        println(io, "  points: EMPTY")
    end

    saxes = ""
    saxes = _mappedgraph_tree!(saxes, frames_axes(fs))
    if saxes != ""
        println(io, "  axes: $(join(["\t"*si for si in split(saxes, "\n")], "\n"))")
    else 
        println(io, "  axes: EMPTY")
    end
    println(io, ")")
end