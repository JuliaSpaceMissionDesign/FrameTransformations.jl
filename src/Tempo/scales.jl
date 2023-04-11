export TIMESCALES, @timescale, add_timescale, TimeSystem

struct TimeScaleNode{T} <: AbstractGraphNode
    name::Symbol
    id::Int
    parentid::Int
    ffp::FunctionWrapper{T,Tuple{T}}
    ftp::FunctionWrapper{T,Tuple{T}}
end

get_node_id(s::TimeScaleNode) = s.id

function Base.show(io::IO, s::TimeScaleNode{T}) where {T}
    pstr = "TimeScaleNode{$T}(name=$(s.name), id=$(s.id)"
    s.parentid == s.id || (pstr *= ", parent=$(s.parentid)")
    pstr *= ")"
    return println(io, pstr)
end

struct TimeSystem{T<:Number}
    scales::MappedNodeGraph{TimeScaleNode{T},SimpleDiGraph{Int}}
    function TimeSystem{T}() where {T}
        return new(MappedDiGraph(TimeScaleNode{T}))
    end
end

timescale_alias(s::AbstractTimeScale) = timescale_id(s)
timescale_alias(s::Int) = s
timescale_id(::AbstractTimeScale) = nothing
timescale_name(::AbstractTimeScale) = nothing

"""
    @timescale(name, id, type)

Create an `AbstractTimeScale` instance to alias the given `id`.

### Inputs 
- `name` -- Acronym to denote the new time scale 
- `id` -- Integer whose alias is the new time scale 
- `type` -- Name of the type representing the new time scale 
"""
macro timescale(name::Symbol, id::Int, type::Symbol)
    type = Symbol(format_camelcase(Symbol, String(type)))
    type_str = String(type)
    name_split = join(split(type_str, r"(?=[A-Z])"), " ")
    name_str = String(name)

    scaleid_expr = :(@inline Tempo.timescale_id(::$type) = $id)
    name_expr = :(@inline Tempo.timescale_name(::$type) = Symbol($name_str))
    show_expr = :(Base.show(io::IO, ::$type) = print(io, "$($(name_str))"))

    return quote
        """
        $($type_str)

        A type representing the $($name_split) ($($name_str)) time scale.
        """
        struct $type <: AbstractTimeScale end

        """
            $($name_str)

        The singleton instance of the [`$($type_str)`](@ref) type representing
        the $($name_split) ($($name_str)) time scale.
        """
        const $(esc(name)) = $(esc(type))()

        $(esc(scaleid_expr))
        $(esc(name_expr))
        $(esc(show_expr))
        nothing
    end
end

function _zero_offset(seconds::T) where {T}
    @error "a zero-offset transformation has been applied in the TimeSystem"
    return T(0.0)
end

"""
    add_vertex!(s::TimeSystem, ts::TimeScaleNode)

Register a new node in the `TimeSystem`.
"""
function add_vertex!(s::TimeSystem{T}, ts::TimeScaleNode{T}) where {T}
    return add_vertex!(s.scales, ts)
end

@inline has_timescale(s::TimeSystem, sid::Int) = has_vertex(s.scales, sid)

function add_timescale(
    scales::TimeSystem{T}, ts::S, ffp::Function=_zero_offset; ftp=nothing, parent=nothing
) where {T,S<:AbstractTimeScale}
    name, id = Tempo.timescale_name(ts), Tempo.timescale_id(ts)
    pid = isnothing(parent) ? nothing : timescale_alias(parent)

    if has_timescale(scales, id)
        # Check if a set of timescale with the same ID is already registered within 
        # the given time system 
        throw(
            ErrorException(
                "[Tempo] TimeScale with id $id is already registered in the given TimeSystem",
            ),
        )
    end

    if name in map(x -> x.name, scales.scales.nodes)
        # Check if timescale with the same name also does not already exist
        throw(
            ErrorException(
                "[Tempo] TimeScale with name $name is already registered in the given TimeSystem",
            ),
        )
    end

    # if the scale has a parent
    if !isnothing(parent)
        # Check if the root axes is not present
        isempty(scales.scales) && throw(ErrorException("[Tempo] missing root timescale"))

        # Check if the parent scale are registered in system
        if !has_timescale(scales, pid)
            throw(
                ErrorException(
                    "[Tempo] the specified parent timescale with id $pid is not " *
                    "registered in the given TimeSystem",
                ),
            )
        end
    else
        if !isempty(scales.scales)
            throw(ErrorException("[Tempo] TimeSystem has already a root point!"))
        end
        pid = id
    end

    # Create a new node 
    tsnode = TimeScaleNode{T}(name, id, pid, ffp, isnothing(ftp) ? _zero_offset : ftp)

    # Insert the new timescale in the graph
    add_vertex!(scales, tsnode)

    # Connect
    if !isnothing(parent)
        # add transformation from parent to new timescale
        add_edge!(scales.scales, pid, id)
        if !isnothing(ftp)
            # add the transformation from new timescale to parent 
            add_edge!(scales.scales, id, pid)
        end
    end
end

function apply_offsets(
    scales::TimeSystem{N}, sec::Number, from::S1, to::S2
) where {N<:Number,S1<:AbstractTimeScale,S2<:AbstractTimeScale}
    return apply_offsets(
        scales, sec, get_path(scales.scales, timescale_alias(from), timescale_alias(to))
    )
end

@inline function apply_offsets(
    scales::TimeSystem{N}, sec::N, path::Vector{Int}
) where {N<:Number}
    # initilize 
    offsec = sec

    tsi = get_mappednode(scales.scales, path[1])
    tsip1 = get_mappednode(scales.scales, path[2])

    offsec += apply_offsets(offsec, tsi, tsip1)

    @inbounds for i in 2:(length(path) - 1)
        tsi = tsip1
        tsip1 = get_mappednode(scales.scales, path[i + 1])
        offsec += apply_offsets(offsec, tsi, tsip1)
    end
    return offsec
end

@inline function apply_offsets(
    sec::N, ts1::TimeScaleNode{N}, ts2::TimeScaleNode{N}
) where {N<:Number}
    if ts1.parentid == ts2.id
        # This is the case in which the inverse transformation (from child to parent)
        # is used 
        offset = ts1.ftp(sec)
    else # ts2.parentid == ts1.id
        # In this case the direct transfromation (from parent to child is used)
        offset = ts2.ffp(sec)
    end
    return offset
end

@inline function apply_offsets(
    ::TimeSystem{N}, sec::N, ::S, ::S
) where {N<:Number,S<:AbstractTimeScale}
    return sec
end

const TIMESCALES_NAMES = (
    :TerrestrialTime,
    :InternationalAtomicTime,
    :CoordinatedUniversalTime,
    :GeocentricCoordinateTime,
    :BarycentricCoordinateTime,
    :BarycentricDynamicalTime,
    :UniversalTime,
    :HighPrecisionBarycentricDynamicalTime,
    :GlobalPositioningSystemTime,
)

const TIMESCALES_ACRONYMS = (:TT, :TAI, :UTC, :TCG, :TCB, :TDB, :UT1, :TDBH, :GPS)

for i in eachindex(TIMESCALES_ACRONYMS)
    @eval begin
        @timescale $(TIMESCALES_ACRONYMS[i]) $i $(TIMESCALES_NAMES[i])
        export $(TIMESCALES_ACRONYMS[i]), $(TIMESCALES_NAMES[i])
    end
end

const TIMESCALES::TimeSystem{Float64} = TimeSystem{Float64}()

"""
    TIMESCALES

Time scales graph, containing at least: $(String.(TIMESCALES_ACRONYMS))

It can be easily extended using the `@timescale` to create new `TimeScaleNode` aliases and 
`add_timescale` method to define its relation with the other nodes in the graph. 

### Example

```@example
# define a new timescale type alias
@timescale NTS 100 NewTimeScale

# define offset to and from another timescale in the graph 
offset_ffp(seconds) = 1.0
offset_ftp(seconds) = -1.0

# connect to the graph, with the parent node (TDB in this example)
add_timescale(TIMESCALES, NTS, offset_ffp, parent=TDB, ftp=offset_ftp)
```
"""
TIMESCALES

# Create time scales graph
add_timescale(TIMESCALES, TT, _zero_offset)
add_timescale(TIMESCALES, TDB, offset_tt2tdb; parent=TT, ftp=offset_tdb2tt)
add_timescale(TIMESCALES, TAI, offset_tt2tai; parent=TT, ftp=offset_tai2tt)
add_timescale(TIMESCALES, TCG, offset_tt2tcg; parent=TT, ftp=offset_tcg2tt)
add_timescale(TIMESCALES, TCB, offset_tdb2tcb; parent=TDB, ftp=offset_tcb2tdb)
add_timescale(TIMESCALES, UTC, offset_tai2utc; parent=TAI, ftp=offset_utc2tai)
add_timescale(TIMESCALES, TDBH, offset_tt2tdbh; parent=TT)
add_timescale(TIMESCALES, GPS, offset_gps2tai; parent=TAI, ftp=offset_tai2gps)

# Note: the UT1 timescale is added in the Orient module to exploit the IERS_EOP data. 
