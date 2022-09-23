export Epoch, 
       j2000, j2000s, j2000c, 
       second, timescale, value

"""
    Epoch <: AbstractDateTimeEpoch

A type to represent Epoch-like data. 
Epochs here are represented as seconds + fraction since a refence epoch, which 
is considered to be `2000-01-01T12:00:00`, i.e. [`J2000`](@ref).

### Fields 

- `scale` -- `TimeScale` to represent the date 
- `second` -- seconds since `origin`
- `fraction` -- fraction of seconds

### Constructor 

`Epoch{S}(second::N, frac::T) where {N<:Number, S<:TimeScale, T<:AbstractFloat}`
"""
struct Epoch{S, N, T} <: AbstractDateTimeEpoch
    scale::S 
    second::N
    fraction::T
    function Epoch{S}(second::N, frac::T) where {N<:Integer, S<:TimeScale, 
                                                 T<:AbstractFloat}
        return new{S, N, T}(S(), second, frac)
    end
end

"""
    second(ep::Epoch)

Seconds from origin.
"""
second(ep::Epoch) = ep.second

"""
    timescale(ep::Epoch)

Epoch timescale.
"""
timescale(ep::Epoch) = ep.scale

"""
    value(ep::Epoch)

Full `Epoch` value.
"""
value(ep::Epoch) = ep.second + ep.fraction

""" 
    Epoch(sec::T, scale::S) where {S<:TimeScale, T<:AbstractFloat}
    Epoch(sec::N, frac::T, scale::S) where {N<:Integer, S<:TimeScale, 
        T<:AbstractFloat}

Construct an `Epoch` from seconds and fraction since J2000 in the scale `S`.
"""
function Epoch(sec::N, frac::T, scale::S) where {N<:Integer, S<:TimeScale, 
    T<:AbstractFloat}
    return Epoch{typeof(scale)}(sec, frac)
end

function Epoch(sec::T, scale::S) where {S<:TimeScale, T<:AbstractFloat}
    ssec = floor(Int64, sec)
    return Epoch(ssec, sec-ssec, scale)
end

"""
    Epoch{S}(js::T where {T<:Number, S<:TimeScale}

Construct an `Epoch` from seconds since J2000 in the scale `S`.
"""
function Epoch{S}(js::T) where {T<:Number, S<:TimeScale}
    sec = floor(Int, js)
    frac = js - sec
    return Epoch{S}(sec, frac)
end

"""
    Epoch(dt::DateTime, scale::S) where {S<:TimeScale}

Construct an `Epoch` from [`DateTime`](@ref) type and timescale `S` 
"""
function Epoch(dt::DateTime, scale::S) where {S<:TimeScale}
    Epoch(j2000s(dt), scale)
end

"""
    j2000(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (days)
"""
j2000(e::Epoch) = e.second/DAY2SEC + e.fraction/DAY2SEC

"""
    j2000s(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (seconds)
"""
j2000s(e::Epoch) = value(e)

"""
    j2000c(e::Epoch)

Convert `Epoch` in Julian Date since J2000 (centuries)
"""
j2000c(e::Epoch) = value(e)/CENTURY2SEC

function Base.show(io::IO, ep::Epoch) 
    print(io, DateTime(ep), " ", timescale(ep))
end

function DateTime(ep::Epoch) 
    # TODO: is there a more elegant way to do this? 
    # There are cases where there is a loss of precision
    days = ep.second ÷ DAY2SEC
    frac = ep.fraction + (ep.second - DAY2SEC*(days - 0.5))
    if frac > 0
        fint = frac ÷ DAY2SEC
        days += fint 
        frac -= fint * DAY2SEC
    elseif frac < 0
        fint = (frac-DAY2SEC) ÷ DAY2SEC
        days += fint 
        frac -= fint * DAY2SEC
    end
    y, m, d = jd2cal(DJ2000, days)
    H, M, S, f = fd2hmsf(frac/DAY2SEC)
    fint = floor(Int64, f)
    S += fint 
    f -= fint
    DateTime(y, m, d, H, M, S, f)
end

Epoch(e::Epoch) = e 
Epoch{S, N, T}(e::Epoch{S, N, T}) where {S, N, T} = e

"""
    Epoch(s::AbstractString, scale::S) where {S<:TimeScale}

Construct an `Epoch` from a `str` in the [`TimeScale`](@ref) (`scale`). 

!!! note 
    This constructor requires that the `str` is in the format `yyyy-mm-ddTHH:MM:SS.sss`.
"""
function Epoch(s::AbstractString, scale::S) where {S<:TimeScale}
    y, m, d, H, M, s, sf = parse_iso(s)
    _, jd2 = calhms2jd(y, m, d, H, M, s+sf)
    days = floor(Int64, jd2)
    fd = jd2 - days 
    h, m, s, frac = fd2hmsf(fd)
    sec = s + 60 * (m + 60*(h + 24*days))
    return Epoch(sec, frac, scale)
end

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

julia> # Parse Julian Dates since J200 
julia> Epoch("12.0")
2000-01-13T12:00:00.0000 TDB

julia> # All Julian Date parsers allow timescales 
julia> Epoch("12.0 TT")
2000-01-13T12:00:00.0000 TT
```
"""
function Epoch(s::AbstractString)
    scale = TDB # default timescale

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
    return Epoch(parse(Float64, sub[1])*86400.0, scale)
end

# ----
# Operations

function Base.:-(e1::Epoch{S}, e2::Epoch{S}) where S
    return (e1.second - e2.second) + (e1.fraction - e2.fraction)
end
function Base.:+(e::Epoch, x::N) where {N<:Number} 
    xint = floor(Int64, x)
    sec = e.second + xint
    frac = e.fraction + (x - xint)
    Epoch(sec, frac, timescale(e))
end
function Base.:-(e::Epoch, x::N) where {N<:Number} 
    xint = floor(Int64, x)
    sec = e.second - xint
    frac = e.fraction - (x - xint)
    Epoch(sec, frac, timescale(e))
end

function (::Base.Colon)(start::Epoch, step::T, stop::Epoch) where {T<:AbstractFloat}
    step = start < stop ? step : -step
    StepRangeLen(start, step, floor(Int64, (stop-start)/step)+1)
end

(::Base.Colon)(start::Epoch, stop::Epoch) = (:)(start, 86400.0, stop)

function Base.isless(e1::Epoch{S}, e2::Epoch{S}) where {S}
    return value(e1) < value(e2)
end

function Base.isapprox(e1::Epoch{S}, e2::Epoch{S}; kwargs...) where {S}
    return isapprox(value(e1), value(e2); kwargs...)
end