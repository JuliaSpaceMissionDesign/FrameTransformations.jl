export FrameSystem,
    ComputableAxesVector,
    frames_order,
    frames_timescale,
    frames_axes,
    frames_points,
    show_axes,
    show_points,
    add_axes!,
    add_point!,
    has_axes,
    has_point,
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
```julia-repl
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

        return new(from, to, order)
    end
end

ComputableAxesVector() = ComputableAxesVector(1, 2, 1)
function ComputableAxesVector(from, to, order::Int)
    return ComputableAxesVector(point_alias(from), point_alias(to), order)
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

function ComputableAxesProperties()
    return ComputableAxesProperties(ComputableAxesVector(), ComputableAxesVector())
end

# Frame Axes Function signatures definition, supporting up to the 2nd derivative 
# without allocations 
_TagAD1{T} = Autodiff.ForwardDiff.Tag{Autodiff.JSMDDiffTag, T}
_NodeFunAD1{T} = Autodiff.ForwardDiff.Dual{_TagAD1{T}, T, 1}

_TagAD2{T} = Autodiff.ForwardDiff.Tag{Autodiff.JSMDDiffTag, _NodeFunAD1{T}}
_NodeFunAD2{T} = Autodiff.ForwardDiff.Dual{_TagAD2{T}, _NodeFunAD1{T}, 1}

_FAxesFunIn{N, T, S} = Tuple{T, SVector{N,S}, SVector{N,S}}
_FAxesFunSig{O, T, S, N} = FunctionWrapper{Rotation{O, T}, _FAxesFunIn{N, T, S}}

_FAxesWrappers{O, T, N} = FunctionWrappersWrapper{Tuple{
    _FAxesFunSig{O, T, T, N},
    _FAxesFunSig{O, _NodeFunAD1{T}, T, N},
    _FAxesFunSig{O, _NodeFunAD1{T}, _NodeFunAD1{T}, N},
    _FAxesFunSig{O, _NodeFunAD2{T}, T, N}, 
    _FAxesFunSig{O, _NodeFunAD2{T}, _NodeFunAD2{T}, N}
}, true}


# This automatically generates the FunctionWrappersWrapper according to the above type 
# definitions for the given input function
function _FAxesWrappers{O, T, N}(fun::Function) where {O, T, N}
    
    typesT = (T, _NodeFunAD1{T}, _NodeFunAD1{T}, _NodeFunAD2{T}, _NodeFunAD2{T})
    typesS = (T, T, _NodeFunAD1{T}, T, _NodeFunAD2{T})

    inps = map((x,y)->_FAxesFunIn{N, x, y}, typesT, typesS)
    outs = map(x->Rotation{O, x}, typesT)

    fws = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end

    return _FAxesWrappers{O, T, N}(fws)

end

# Container to store frame axes update functions 
struct FrameAxesFunctions{T,O,N}
    fun::NTuple{O, _FAxesWrappers{O,T,N}}
end

Base.getindex(af::FrameAxesFunctions, i) = af.fun[i]

# Default rotation function for axes that do not require updates
_get_fixedrot(::T, ::SVector{N,T}, ::SVector{N,T}) where {N,T} = Rotation{N / 3}(T(1)I)

function _get_fixedrot(::T, ::SVector{N, S}, ::SVector{N, S}) where {N, S, T}
    Rotation{N/3}(promote_type(T, S)(1)I)
end

# Constructors for FrameAxesFunctions 
@generated function FrameAxesFunctions{T}(funs::Function...) where {T}
    O = length(funs)

    expr = :(tuple($([
        Expr(:call, 
            Expr(:curly, :_FAxesWrappers, :O, :T, Expr(:call, :(*), 3, :O)), 
            Expr(:ref, :funs, i)) for i in 1:O
        ]...))
    )

    return quote
        O = length(funs)
        @inbounds FrameAxesFunctions{T,$O,3 * $O}($(expr))
    end
end

# Constructor to filter out some of the specified functions!
@generated function FrameAxesFunctions{T,O}(funs::Function...) where {T,O}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = :(tuple($([
        Expr(:call, 
            Expr(:curly, :_FAxesWrappers, :O, :T, Expr(:call, :(*), 3, :O)), 
            Expr(:ref, :funs, i)) for i in 1:O
        ]...))
    )

    return quote
        @inbounds FrameAxesFunctions{T,O,3 * O}($(expr))
    end
end

# Default constructors for dummy axes function updates 
@generated function FrameAxesFunctions{T,O}() where {T,O}
    expr = :(tuple($([_FAxesWrappers{O, T, 3*O}(_get_fixedrot) for i in 1:O]...)))
    return quote
        FrameAxesFunctions{T,O,3 * O}($(expr))
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
- `R` -- rotation matrix for fixed relative axes 
- `f` -- `FrameAxesFunctions` container 
- `angles` -- vector storing the libration angles retrived from ephemerides
"""
struct FrameAxesNode{O,T,N} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    id::Int
    parentid::Int
    comp::ComputableAxesProperties
    R::Vector{Rotation{O, T}}
    epochs::Vector{T}
    nzo::Vector{Int}
    f::FrameAxesFunctions{T,O,N}
    angles::Vector{DiffCache{MVector{N,T}, Vector{T}}}
end         

get_node_id(ax::FrameAxesNode) = ax.id

function Base.show(io::IO, ax::FrameAxesNode{O,T}) where {O,T}
    pstr = "FrameAxesNode{$O, $T}(name=$(ax.name), class=$(ax.class), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

# -------------------------------------
# POINTS
# -------------------------------------

# Frame Point Function signatures definition, supporting up to the 2nd derivative 
# without allocations 
_FPointFunIn{N, T} = Tuple{MVector{N, T}, T}
_FPointFunSig{N, T} = FunctionWrapper{Nothing, _FPointFunIn{N, T}}

_FPointWrappers{N, T} = FunctionWrappersWrapper{Tuple{
    _FPointFunSig{N, T},
    _FPointFunSig{N, _NodeFunAD1}, 
    _FPointFunSig{N, _NodeFunAD2}
}, true}

# This automatically generates the FunctionWrappersWrapper according to the above type 
# definitions for the given input function
function _FPointWrappers{N, T}(fun::Function) where {N, T}

    types = (T, _NodeFunAD1{T}, _NodeFunAD2{T})
    
    inps = map(x->_FPointFunIn{N, x}, types)
    outs = map(x->Nothing, types)

    fws = map(inps, outs) do A, R 
        FunctionWrapper{R, A}(fun)
    end

    return _FPointWrappers{N, T}(fws)

end

# Container to store frame point update functions 
struct FramePointFunctions{T,O,N}
    fun::NTuple{O,_FPointWrappers{N, T}}
end

Base.getindex(pf::FramePointFunctions, i) = pf.fun[i]

# Default state-vector update function for points that do not require updates
_empty_stv_update!(::AbstractVector{T}, ::T) where {T} = nothing

# Constructors for FramePointFunctions 
@generated function FramePointFunctions{T}(funs::Function...) where {T}
    O = length(funs)

    expr = :(tuple($([
        Expr(:call, 
            Expr(:curly, :_FPointWrappers, Expr(:call, :(*), 3, :O), :T), 
            Expr(:ref, :funs, i)) for i in 1:O
        ]...))
    )

    return quote
        O = length(funs)
        @inbounds FramePointFunctions{T, $O, 3*$O}($(expr))
    end
end

# Constructor to filter out some of the specified functions!
@generated function FramePointFunctions{T,O}(funs::Function...) where {T,O}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))

    expr = :(tuple($([
        Expr(:call, 
            Expr(:curly, :_FPointWrappers, Expr(:call, :(*), 3, :O), :T), 
            Expr(:ref, :funs, i)) for i in 1:O
        ]...))
    )
    return quote
        @inbounds FramePointFunctions{T, O, 3*O}($(expr))
    end
end

# Default constructors for dummy point function updates 
@generated function FramePointFunctions{T,O}() where {T,O}
    expr = :(tuple($([_FPointWrappers{3*O, T}(_empty_stv_update!) for i in 1:O]...)))
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
struct FramePointNode{O,T,N} <: AbstractGraphNode
    name::Symbol
    class::Symbol
    axesid::Int
    parentid::Int
    NAIFId::Int
    stv::Vector{DiffCache{MVector{N,T}, Vector{T}}}
    epochs::DiffCache{Vector{T}, Vector{T}}
    nzo::Vector{MVector{2, Int}}
    f::FramePointFunctions{T,O,N}
end

get_node_id(p::FramePointNode) = p.NAIFId

function Base.show(io::IO, p::FramePointNode{O,T}) where {O,T}
    pstr = "FramePointNode{$O, $T}(name=$(p.name), class=$(p.class), NAIFId=$(p.NAIFId), axes=$(p.axesid)"
    p.parentid == p.NAIFId || (pstr *= ", parent=$(p.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

# -------------------------------------
# FRAMES
# -------------------------------------

struct FrameSystemProperties
    ebid::Vector{Int}  # ephemeris body ids
    eaid::Vector{Int}  # ephemeris axes ids 
end
FrameSystemProperties() = FrameSystemProperties(Int64[], Int64[])

@inline ephemeris_points(fsp::FrameSystemProperties) = fsp.ebid
@inline ephemeris_axes(fsp::FrameSystemProperties) = fsp.eaid

"""
    FrameSystem{O, T, S, E}

A `FrameSystem` instance manages a collection of user-defined `FramePointNode` and 
`FrameAxesNode` objects, enabling efficient computation of arbitrary transformations 
between them. It is created by specifying the maximum transformation order `O`, the outputs 
datatype `T` and an `AbstractTimeScale` instance `S`. Additionally, an `AbstractEphemerisProvider` 
instance `E` can be provided to compute transformations that require ephemeris data. 

The following transformation orders are accepted: 
- **1**: position 
- **2**: position and velocity 
- **3**: position, velocity and acceleration
- **4**: position, velocity, acceleration and jerk

By specifying the maximum transformation the `FrameSystem` memory usage and performance can 
be optimised and tailored to the user's needs.

--- 

    FrameSystem{O, T}()

Create a `FrameSystem` object of order `O` and datatype `T`. The `BarycentricDynamicalTime` 
is automatically assigned as the default time scale. The resulting object is constructed 
with a `NullEphemerisProvider`, which does not allow the computation of transformation that 
involve ephemeris files.

### Examples 
```julia-repl
julia> F = FrameSystem{2, Float64}();

julia> @axes ICRF 1 

julia> @axes ECLIPJ2000 17 

julia> add_axes_inertial!(F, ICRF)

julia> add_axes_eclipj2000!(F, ECLIPJ2000, ICRF)

julia> rotation6(F, ICRF, ECLIPJ2000, 0.0)
Rotation{2, Float64}
[...]

julia> rotation9(F, ICRF, ECLIPJ2000, 0.0)
ERROR: Insufficient frame system order: transformation requires at least order 3.

```

---

    FrameSystem{O, T, S}() 

Create a `FrameSystem` object of order `O`, datatype `T` and time scale `S`. The resulting 
object is constructed with a `NullEphemerisProvider`, which does not allow the computation 
of transformation that involve ephemeris files.

### Examples 

```julia-repl
julia> F = FrameSystem{1, Float64, TerrestrialTime}();

julia> @axes ICRF 1 

julia> @axes ECLIPJ2000 17 

julia> add_axes_inertial!(F, ICRF)

julia> add_axes_eclipj2000!(F, ECLIPJ2000, ICRF)

julia> ep_tt = Epoch("2023-02-10T12:00:00 TT")
2023-02-10T12:00:00.000 TT

julia> rotation3(F, ICRF, ECLIPJ2000, ep_tt)
Rotation{1,Float64}([...])

julia> ep_tdb = Epoch("2023-02-10T12:00:00 TDB")
2023-02-10T12:00:00.000 TDB

julia> rotation3(F, ICRF, ECLIPJ2000, ep_tdb)
ERROR: ArgumentError: Incompatible epoch timescale: expected TerrestrialTime, found BarycentricDynamicalTime.
[...]
```
---

    FrameSystem{O, T}(eph::AbstractEphemerisProvider)

Create a `FrameSystem` object of order `O` and datatype `T` by providing an instance of an 
`AbstractEphemerisProvider` subtype. The timescale is automatically set to the one associated 
to the ephemeris files loaded in `eph`. This constructor shall be used when the user desires 
to compute transformations that involve ephemeris data. 

!!! note 
    All the kernels that will be used must be loaded within `eph`. Once the `FrameSystem` 
    has been created, no additional kernel can be added nor removed.

### Examples 
```julia-repl
julia> using Ephemerides 

julia> eph = EphemerisProvider(DE440_KERNEL_PATH);

julia> F = FrameSystem{2, Float64}(eph)
FrameSystem{2, Float64, BarycentricDynamicalTime, EphemerisProvider}(
  eph: 1-kernel EphemerisProvider,
  points: EMPTY
  axes: EMPTY
)
```

### See also 
See also [`add_axes_inertial!`](@ref), [`add_point_root!`](@ref), [`vector3`](@ref) and [`rotation3`](@ref)
"""
struct FrameSystem{O,T<:Number,S<:AbstractTimeScale,E<:AbstractEphemerisProvider,N}
    eph::E
    prop::FrameSystemProperties
    points::MappedNodeGraph{FramePointNode{O,T,N},SimpleGraph{Int}}
    axes::MappedNodeGraph{FrameAxesNode{O,T,N},SimpleGraph{Int}}
end

function FrameSystem{O,T,S}(
    eph::E, points::Vector{Int}, axes::Vector{Int}
) where {O,T<:Number,S<:AbstractTimeScale,E<:AbstractEphemerisProvider}
    if O < 1 || O > 4
        throw(ArgumentError("FrameSystem order must be between 1 and 4."))
    end

    return FrameSystem{O,T,S,E,3 * O}(
        eph,
        FrameSystemProperties(points, axes),
        MappedGraph(FramePointNode{O,T,3O}),
        MappedGraph(FrameAxesNode{O,T,3O}),
    )
end

function FrameSystem{O,T}(eph::E) where {O,T,E}
    points = ephem_available_points(eph)
    axes = ephem_available_axes(eph)

    tsid = ephem_timescale(eph)

    if tsid == 1
        S = TDB
    elseif tsid == 2
        S = TCB
    else
        throw(EphemerisError("Unrecognised ephemeris timescale."))
    end

    return FrameSystem{O,T,typeof(S)}(eph, points, axes)
end

function FrameSystem{O,T,S}() where {O,T,S}
    return FrameSystem{O,T,S}(NullEphemerisProvider(), Int64[], Int64[])
end

FrameSystem{O,T}() where {O,T<:Number} = FrameSystem{O,T,BarycentricDynamicalTime}()

frames_order(::FrameSystem{O}) where {O} = O
frames_timescale(::FrameSystem{O,T,S}) where {O,T,S} = S

frames_points(fs::FrameSystem) = fs.points
frames_axes(fs::FrameSystem) = fs.axes

function add_point!(fs::FrameSystem{O,T}, p::FramePointNode{O,T}) where {O,T}
    return add_vertex!(fs.points, p)
end

function add_axes!(fs::FrameSystem{O,T}, ax::FrameAxesNode{O,T}) where {O,T}
    return add_vertex!(fs.axes, ax)
end

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
    return nothing
end

function _mappedgraph_tree!(s::String, g::MappedNodeGraph)
    if !isempty(g.nodes)
        s *= "\n$(g.nodes[1].name)\n"
        s = _mappedgraph_tree!(s, g, get_node_id(g.nodes[1]), 2, 1)
    end
    return s
end

function _mappedgraph_tree!(s::String, g::MappedNodeGraph, pid::Int, idx::Int, del::Int=1)
    @inbounds for i in idx:length(g.nodes)
        if g.nodes[i].parentid == pid
            s *= "$(" "^(del))├── $(g.nodes[i].name) \n"
            s = _mappedgraph_tree!(s, g, get_node_id(g.nodes[i]), i, del + 1)
        end
    end
    return s
end

function Base.summary(io::IO, ::FrameSystem{O,T,S,E}) where {O,T,S,E}
    return println(io, "FrameSystem{$O, $T, $S, $E}")
end

function Base.show(io::IO, fs::FrameSystem{O,T,S,E}) where {O,T,S,E}
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
    return println(io, ")")
end
