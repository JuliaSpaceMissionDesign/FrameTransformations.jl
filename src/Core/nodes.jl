
const SVectorNT{D, T} = SVector{D, T}

@generated function SVectorNT{D, T}() where {D, T}
    expr = Expr(:call, Expr(:curly, :SVector, D, T))
    for _ in 1:D 
        push!(expr.args, T(0))
    end
    return quote
        @inbounds $(expr)
    end
end

@generated function SVectorNT{D, T}(vec::SVector{N, T}) where {D, T, N}
    expr = Expr(:call, Expr(:curly, :SVector, D, T))
    for i in 1:min(D, N)
        push!(expr.args, Expr(:ref, :vec, i))
    end
    for _ in 1:(D-N) 
        push!(expr.args, Expr(:call, :zero, T))
    end
    return quote
        @inbounds $(expr)
    end
end

# ------------------------------------------------------------------------------------------
# POINTS
# ------------------------------------------------------------------------------------------

struct FramePointFunctions{O, N, D}
    fun::NTuple{O, FramePointFunWrapper{D, N}}
end

Base.getindex(pf::FramePointFunctions, i) = pf.fun[i]

@generated function FramePointFunctions{N}(funs::Function...) where {N}
    O = length(funs)
    D = 3O

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FramePointFunWrapper, D, N),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FramePointFunctions, O, N, D),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FramePointFunctions{O, N}(funs::Function...) where {O, N}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))
    D = 3O

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FramePointFunWrapper, D, N),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FramePointFunctions, O, N, D),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FramePointFunctions{O, N}(fun::Function) where {O, N}
    D = 3O
    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FramePointFunWrapper, D, N),
                :fun
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FramePointFunctions, O, N, D),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

function FramePointFunctions{O, N}() where {O, N}
    return FramePointFunctions{O, N}(t->SVectorNT{3*O, N}())
end

# ------
# Node

"""
    FramePointNode{O, T, N} <: AbstractJSMDGraphNode

Define a frame system point.

### Fields
- `name` -- point name 
- `class` -- `Symbol` representing the class of the point 
- `id` -- ID of the point
- `parentid` -- ID of the parent point 
- `axesid` -- ID of the axes in which the point coordinates are expressed 
- `f` -- `FramePointFunctions` container 
"""
struct FramePointNode{O, N <: Number, D} <: AbstractJSMDGraphNode
    name::Symbol
    class::Int
    id::Int
    parentid::Int
    axesid::Int

    # internals
    f::FramePointFunctions{O, N, D}
end

function FramePointNode{O, N}(
    name::Symbol, class::Int, id::Int, pid::Int, axid::Int, 
    fps::FramePointFunctions{O, N, D}
) where {O, N, D}
    return FramePointNode{O, N, 3O}(name, class, id, pid, axid, fps)
end

get_node_id(p::FramePointNode{O, N}) where {O, N} = p.id

function Base.show(io::IO, p::FramePointNode{O, N, D}) where {O, N, D}
    pstr = "FramePointNode{$O, $N}(name=$(p.name), class=$(p.class)"
    pstr *= ", id=$(p.id), axesid=$(p.axesid)"
    p.parentid == p.id || (pstr *= ", parent=$(p.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

const PointsGraph{O, N, D} = MappedNodeGraph{FramePointNode{O, N, D}, SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# AXES
# ------------------------------------------------------------------------------------------

struct FrameAxesFunctions{O, N, D}
    fun::NTuple{O, FrameAxesFunWrapper{O, N}}
end

Base.getindex(pf::FrameAxesFunctions, i) = pf.fun[i]

@generated function FrameAxesFunctions{N}(funs::Function...) where {N}
    O = length(funs)

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FrameAxesFunWrapper, O, N),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, N, 3O),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FrameAxesFunctions{O, N}(funs::Function...) where {O, N}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FrameAxesFunWrapper, O, N),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, N, 3O),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FrameAxesFunctions{O, N}(fun::Function) where {O, N}
    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FrameAxesFunWrapper, O, N),
                :fun
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, N, 3O),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

function FrameAxesFunctions{O, N}() where {O, N}
    return FrameAxesFunctions{O, N}(t->Rotation{O, N}(N(1)I))
end

# ------
# Node

"""
    FrameAxesNode{O, T, N} <: AbstractJSMDGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `class` -- `Symbol` representing the class of the axes 
- `id` -- axes ID (equivalent of NAIFId for axes)
- `parentid` -- ID of the parent axes 
- `f` -- `FrameAxesFunctions` container 
"""
struct FrameAxesNode{O, N <: Number, D} <: AbstractJSMDGraphNode
    name::Symbol
    class::Int
    id::Int
    parentid::Int

    # internals
    f::FrameAxesFunctions{O, N, D}
end

function FrameAxesNode{O, N}(
    name::Symbol, class::Int, id::Int, pid::Int, faxs::FrameAxesFunctions{O, N, D}
) where {O, N, D}
    return FrameAxesNode{O, N, 3O}(name, class, id, pid, faxs)
end

get_node_id(ax::FrameAxesNode{O, N, D}) where {O, N, D} = ax.id

function Base.show(io::IO, ax::FrameAxesNode{O, N, D}) where {O, N, D}
    pstr = "FrameAxesNode{$O, $N}(name=$(ax.name), class=$(ax.class), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

const AxesGraph{O, N, D} = MappedNodeGraph{FrameAxesNode{O, N, D}, SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# DIRECTIONs
# ------------------------------------------------------------------------------------------

const DirectionFunctions{O, N, D} = FramePointFunctions{O, N, D}

"""
    Direction{O, N, D}

Define a new direction.

### Fields
- `name` -- direction name 
- `id` -- direction ID
- `f` -- `DirectionFunctions` container 
"""
struct Direction{O, N, D}
    name::Symbol 
    id::Int 

    # internals
    f::DirectionFunctions{O, N, D}
end

function Direction{O, N}(
    name::Symbol, id::Int, dfun::DirectionFunctions{O, N, D}
) where {O, N, D}
    return Direction{O, N, 3O}(name, id, dfun)
end

function Base.show(io::IO, d::Direction{O, N, D}) where {O, N, D}
    return println(io, "Direction{$O, $N}(name=$(d.name), id=$(d.id))")
end