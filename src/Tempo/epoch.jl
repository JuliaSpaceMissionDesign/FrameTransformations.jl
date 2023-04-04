export Epoch, j2000, j2000s, j2000c, second, timescale, value

"""
    Epoch

A type to represent Epoch-like data. 
Epochs here are represented as seconds + fraction since a refence epoch, which 
is considered to be `2000-01-01T12:00:00`, i.e. [`J2000`](@ref).

### Fields 

- `scale` -- `TimeScale` to represent the date 
- `seconds` -- seconds since `origin`

### Constructors

- `Epoch{S}(seconds::T) where {S<:AbstractTimeScale, T<:AbstractFloat}` -- default
- `Epoch(seconds::T, ::S) where {S<:AbstractTimeScale, T<:AbstractFloat}` 
- `Epoch(dt::DateTime, ::S) where {S<:AbstractTimeScale}` -- construct from `DateTime`
"""
struct Epoch{S,T}
    scale::S
    seconds::T
    function Epoch{S}(seconds::T) where {S<:AbstractTimeScale,T<:AbstractFloat}
        return new{S,T}(S(), seconds)
    end
end

"""
    timescale(ep::Epoch)

Epoch timescale.
"""
timescale(ep::Epoch) = ep.scale

"""
    value(ep::Epoch)

Full `Epoch` value.
"""
value(ep::Epoch) = ep.seconds

function Epoch(sec::T, ::S) where {S<:AbstractTimeScale,T<:AbstractFloat}
    return Epoch{S}(sec)
end

function Epoch(dt::DateTime, ::S) where {S<:AbstractTimeScale}
    return Epoch{S}(j2000s(dt))
end

function Epoch(sec::T, ::Type{S}) where {S<:AbstractTimeScale,T<:AbstractFloat}
    return Epoch{S}(sec)
end

function Epoch(dt::DateTime, ::Type{S}) where {S<:AbstractTimeScale}
    return Epoch{S}(j2000s(dt))
end

"""
    j2000(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (days)
"""
j2000(e::Epoch) = e.seconds / DAY2SEC

"""
    j2000s(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (seconds)
"""
j2000s(e::Epoch) = value(e)

"""
    j2000c(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (centuries)
"""
j2000c(e::Epoch) = value(e) / CENTURY2SEC

function Base.show(io::IO, ep::Epoch)
    return print(io, DateTime(ep), " ", timescale(ep))
end

function DateTime(ep::Epoch)
    y, m, d, H, M, S = jd2calhms(DJ2000, ep.seconds / DAY2SEC)
    fint = floor(Int64, S)
    f = S - fint
    return DateTime(y, m, d, H, M, fint, f)
end

Epoch(e::Epoch) = e
Epoch{S,T}(e::Epoch{S,T}) where {S,T} = e

"""
    Epoch(s::AbstractString, scale::S) where {S<:AbstractTimeScale}

Construct an `Epoch` from a `str` in the [`AbstractTimeScale`](@ref) (`scale`). 

!!! note 
    This constructor requires that the `str` is in the format `yyyy-mm-ddTHH:MM:SS.sss`.
"""
function Epoch(s::AbstractString, scale::S) where {S<:AbstractTimeScale}
    y, m, d, H, M, sec, sf = parse_iso(s)
    _, jd2 = calhms2jd(y, m, d, H, M, sec + sf)
    return Epoch(jd2 * 86400, scale)
end

# This is the default timescale used by Epochs 
const DEFAULT_EPOCH_TIMESCALE = :TDB

"""
    Epoch(s::AbstractString)

Construct an `Epoch` from a `str`.

!!! note 
    This constructor requires that the `str` is in the format:
    - **ISO** -- `yyyy-mm-ddTHH:MM:SS.ffff` : assume J2000 as origin
    - **J2000** -- `DDDD.ffff` : parse Julian Dates since J2000 (days)
    - **JD** -- `JD DDDDDDDDD.ffffff` : parse Julian Date (days)
    - **MJD** -- `MJD DDDDDDDDD.ffffff` : parse Julian Dates since (days)
    In all cases the `TimeScale` should be added at the end of the string, 
    separated by a whitespace. If it is not declated, [`TDB`](@ref) will be 
    used as timescale. 

### Examples
```@example 
julia> # standard ISO string 
julia> Epoch("2050-01-01T12:35:15.0000 TT")
2050-01-01T12:35:14.9999 TT

julia> # standard ISO string (without scale)
julia> Epoch("2050-01-01T12:35:15.0000")
2050-01-01T12:35:14.9999 TDB

julia> # Parse Julian Dates 
julia> Epoch("JD 2400000.5")
1858-11-17T00:00:00.0000 TDB

julia> # Parse Modified Julian Dates
julia> Epoch("MJD 51544.5")
2000-01-01T12:00:00.0000 TDB

julia> # Parse Julian Dates since J2000 
julia> Epoch("12.0")
2000-01-13T12:00:00.0000 TDB

julia> # All Julian Date parsers allow timescales 
julia> Epoch("12.0 TT")
2000-01-13T12:00:00.0000 TT
```
"""
function Epoch(s::AbstractString)
    scale = eval(DEFAULT_EPOCH_TIMESCALE) # default timescale

    # ISO 
    m = match(r"\d{4}-", s)
    if !isnothing(m) && length(m.match) != 0
        sub = split(s, " ")
        if length(sub) == 2 # check for timescale
            scale = eval(Symbol(sub[2]))
        end
        return Epoch(sub[1], scale)
    end

    # JD
    m = match(r"JD", s)
    mjd = match(r"MJD", s)
    if !isnothing(m) && isnothing(mjd) && length(m.match) != 0
        sub = split(s, " ")
        if length(sub) == 3 # check for timescale
            scale = eval(Symbol(sub[3]))
        end
        days = parse(Float64, sub[2])
        sec = (days - DJ2000) * DAY2SEC
        return Epoch(sec, scale)
    end

    # MJD
    m = mjd
    if !isnothing(m) && length(m.match) != 0
        sub = split(s, " ")
        if length(sub) == 3 # check for timescale
            scale = eval(Symbol(sub[3]))
        end
        days = parse(Float64, sub[2])
        sec = (days - DMJD) * DAY2SEC
        return Epoch(sec, scale)
    end

    # J2000
    sub = split(s, " ")
    if length(sub) == 2 # check for timescale
        scale = eval(Symbol(sub[2]))
    end
    return Epoch(parse(Float64, sub[1]) * 86400.0, scale)
end

# ----
# Operations

function Base.:-(e1::Epoch{S}, e2::Epoch{S}) where {S}
    return e1.seconds - e2.seconds
end

function Base.:+(e::Epoch, x::N) where {N<:Number}
    return Epoch(e.seconds + x, timescale(e))
end
function Base.:-(e::Epoch, x::N) where {N<:Number}
    return Epoch(e.seconds - x, timescale(e))
end

function (::Base.Colon)(start::Epoch, step::T, stop::Epoch) where {T<:AbstractFloat}
    step = start < stop ? step : -step
    return StepRangeLen(start, step, floor(Int64, (stop - start) / step) + 1)
end

(::Base.Colon)(start::Epoch, stop::Epoch) = (:)(start, 86400.0, stop)

function Base.isless(e1::Epoch{S}, e2::Epoch{S}) where {S}
    return value(e1) < value(e2)
end

function Base.isapprox(e1::Epoch{S}, e2::Epoch{S}; kwargs...) where {S}
    return isapprox(value(e1), value(e2); kwargs...)
end

Base.convert(::Type{S}, e::Epoch{S}) where {S<:AbstractTimeScale} = e
Base.convert(::S, e::Epoch{S}) where {S<:AbstractTimeScale} = e

function Base.convert(
    to::S2, e::Epoch{S1}; system::TimeSystem=TIMESCALES
) where {S1<:AbstractTimeScale,S2<:AbstractTimeScale}
    return Epoch{S2}(apply_offsets(system, e.seconds, timescale(e), to))
end
