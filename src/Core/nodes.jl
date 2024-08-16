# ------------------------------------------------------------------------------------------
# POINTS
# ------------------------------------------------------------------------------------------

# ------
# Functions

struct FramePointFunctions{O, T}
    fun::NTuple{O, FramePointFunWrapper{O, T}}
end

Base.getindex(pf::FramePointFunctions, i) = pf.fun[i]

@generated function FramePointFunctions{T}(funs::Function...) where {T}
    O = length(funs)

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call, Expr(:curly, :FramePointFunWrapper, O, T), Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(:call, Expr(:curly, :FramePointFunctions, O, T), expr)

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FramePointFunctions{O, T}(funs::Function...) where {O, T}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call, Expr(:curly, :FramePointFunWrapper, O, T), Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(:call, Expr(:curly, :FramePointFunctions, O, T), expr)

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FramePointFunctions{O, T}(fun::Function) where {O, T}
    expr = Expr(:call, :tuple)
    for _ in 1:O 
        push!(
            expr.args,
            Expr(
                :call, Expr(:curly, :FramePointFunWrapper, O, T), :fun
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FramePointFunctions, O, T),
        expr
    )

    return quote
        Base.@_inline_meta
        $(pexpr)
    end
end

function FramePointFunctions{O, T}() where {O, T}
    return FramePointFunctions{O, T}(t->Translation{O, T}())
end


# ------
# Node

"""
    FramePointNode{O, T} <: AbstractJSMDGraphNode

Define a frame system point.

### Fields
- `name` -- point name 
- `id` -- ID of the point
- `parentid` -- ID of the parent point 
- `axesid` -- ID of the axes in which the point coordinates are expressed 
- `f` -- `FramePointFunctions` container 
"""
struct FramePointNode{O, T <: Number} <: AbstractJSMDGraphNode
    name::Symbol
    id::Int
    parentid::Int
    axesid::Int

    # internals
    f::FramePointFunctions{O, T}
end

get_node_id(p::FramePointNode{O, T}) where {O, T} = p.id

function Base.show(io::IO, p::FramePointNode{O, T}) where {O, T}
    pstr = "FramePointNode{$O, $T}(name=$(p.name)"
    pstr *= ", id=$(p.id), axesid=$(p.axesid)"
    p.parentid == p.id || (pstr *= ", parent=$(p.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

const PointsGraph{O, T} = MappedNodeGraph{FramePointNode{O, T}, SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# AXES
# ------------------------------------------------------------------------------------------

# ------
# Functions

struct FrameAxesFunctions{O, T}
    fun::NTuple{O, FrameAxesFunWrapper{O, T}}
end

Base.getindex(pf::FrameAxesFunctions, i) = pf.fun[i]

@generated function FrameAxesFunctions{T}(funs::Function...) where {T}
    O = length(funs)

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FrameAxesFunWrapper, O, T),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, T),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FrameAxesFunctions{O, T}(funs::Function...) where {O, T}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = Expr(:call, :tuple)
    for i in 1:O 
        push!(
            expr.args,
            Expr(
                :call,
                Expr(:curly, :FrameAxesFunWrapper, O, T),
                Expr(:ref, :funs, i)
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, T),
        expr
    )

    return quote
        @inbounds $(pexpr)
    end
end

@generated function FrameAxesFunctions{O, T}(fun::Function) where {O, T}
    expr = Expr(:call, :tuple)
    for _ in 1:O 
        push!(
            expr.args,
            Expr(
                :call, Expr(:curly, :FrameAxesFunWrapper, O, T), :fun
            )
        )
    end
    pexpr = Expr(
        :call, 
        Expr(:curly, :FrameAxesFunctions, O, T),
        expr
    )

    return quote
        Base.@_inline_meta
        $(pexpr)
    end
end

function FrameAxesFunctions{O, T}() where {O, T}
    return FrameAxesFunctions{O, T}(t->Rotation{O, T}(one(T)I))
end

# ------
# Node

"""
    FrameAxesNode{O, T} <: AbstractJSMDGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `id` -- axes ID (equivalent of NAIFId for axes)
- `parentid` -- ID of the parent axes 
- `f` -- `FrameAxesFunctions` container 
"""
struct FrameAxesNode{O, T <: Number} <: AbstractJSMDGraphNode
    name::Symbol
    id::Int
    parentid::Int

    # internals
    f::FrameAxesFunctions{O, T}
end

get_node_id(ax::FrameAxesNode{O, T}) where {O, T} = ax.id

function Base.show(io::IO, ax::FrameAxesNode{O, T}) where {O, T}
    pstr = "FrameAxesNode{$O, $T}(name=$(ax.name), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

const AxesGraph{O, T} = MappedNodeGraph{FrameAxesNode{O, T}, SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# DIRECTIONS
# ------------------------------------------------------------------------------------------

const DirectionFunctions{O, T} = FramePointFunctions{O, T}

"""
    Direction{O, T}

Define a new direction.

### Fields
- `name` -- direction name 
- `id` -- direction ID
- `f` -- `DirectionFunctions` container 
"""
struct Direction{O, T}
    name::Symbol 
    id::Int 

    # internals
    f::DirectionFunctions{O, T}
end

function Direction{O, T}(name::Symbol, id::Int, dfun::DirectionFunctions{O, T}) where {O, T}
    return Direction{O, T}(name, id, dfun)
end

function Base.show(io::IO, d::Direction{O, T}) where {O, T}
    return println(io, "Direction{$O, $T}(name=$(d.name), id=$(d.id))")
end