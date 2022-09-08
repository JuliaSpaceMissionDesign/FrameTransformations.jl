using ..Utils: two_sum, apply_offset

export Epoch

"""
    Epoch{S, N, T, O} <: AbstractDateTimeEpoch

A type to represent Epoch-like data. 

### Fields 

- `scale` -- `TimeScale` to represent the date 
- `second` -- seconds since `origin`
- `fraction` -- fraction of seconds
- `error` -- seconds error
- `origin` -- date origin, see `AbstractEpochOrigin` subtypes

### Constructor 

`Epoch{S}(second::N, frac::T, err::T=zero(T), origin::O=J2000) 
where {N<:Number, S<:TimeScale, T<:AbstractFloat, O<:AbstractEpochOrigin}`

"""
struct Epoch{S, N, T, O} <: AbstractDateTimeEpoch
    scale::S 
    second::N
    fraction::T
    error::T 
    origin::O
    function Epoch{S}(
        second::N, 
        frac::T, 
        err::T=zero(T), 
        origin::O=J2000) where {N<:Number, S<:TimeScale, 
                                T<:AbstractFloat, 
                                O<:AbstractEpochOrigin}
        return new{S, N, T, O}(S(), second, frac, err, origin)
    end
end

second(ep::Epoch) = ep.second
timescale(ep::Epoch) = ep.scale
value(ep::Epoch) = ep.second + ep.fraction

function offsets(js1::N1, js2::N2, 
    origin::O) where {N1, N2, O<:AbstractEpochOrigin}
    js1 -= offset(origin)
    sum, residual = two_sum(js1, js2)
    epoch = floor(Int, sum)
    off = (sum - epoch) + residual
    return epoch, off
end

"""
    Epoch{S}(js::T; origin::O=J2000) 
    where {T<:Number, S<:TimeScale, O<:AbstractEpochOrigin}

Construct an `Epoch` from seconds since `origin` in the scale `S`.
"""
function Epoch{S}(js::T; origin::O=J2000) where {T<:Number, S<:TimeScale, O<:AbstractEpochOrigin}
    ep, off = offsets(js, 0.0, origin)
    return Epoch{S}(ep, off, 0.0, origin)
end

"""
    Epoch(js::T, scale::S=J2000; origin::O=J2000) 
    where {T<:Number, S<:TimeScale, O<:AbstractEpochOrigin}

Construct an `Epoch` from seconds since `origin` in the scale `scale`.
"""
function Epoch(js::T, scale::S=J2000; origin::O=J2000) where {T<:Number, S<:TimeScale, O<:AbstractEpochOrigin}
    Epoch{typeof(scale)}(js; origin=origin)
end

"""
    Epoch(dt::DateTime, scale::S=J2000; origin::O=J2000) 
    where {S<:TimeScale, O<:AbstractEpochOrigin}

Construct an `Epoch` from [`DateTime`](@ref).
"""
function Epoch(dt::DateTime, scale::S=J2000; origin::O=J2000) where {S<:TimeScale, O<:AbstractEpochOrigin}
    Epoch(j2000seconds(dt), scale; origin=origin)
end 

Epoch(e::Epoch) = e 

function DateTime(ep::Epoch)
    return DateTime(value(ep))
end

"""
    Epoch(s::AbstractString, scale::S; origin::O=J2000) 
    where {S<:TimeScale, O<:AbstractEpochOrigin}

Construct an `Epoch` from a `str` in the [`TimeScale`](@ref) (`scale`). 

!!! note 
    This constructor requires that the `str` is in the format `yyyy-mm-ddTHH:MM:SS.sss`.
"""
function Epoch(s::AbstractString, scale::S; origin::O=J2000) where {S<:TimeScale, O<:AbstractEpochOrigin}
    Epoch(DateTime(s), scale; origin=origin)
end

"""
    Epoch(s::AbstractString; origin::O=J2000) where {O<:AbstractEpochOrigin}

Construct an `Epoch` from a `str`.

!!! note 
    This constructor requires that the `str` is in the format:
    - **ISO** -- `yyyy-mm-ddTHH:MM:SS.ssss ttt`
    - **JD** -- `JD DDDDDDDD.ddd`
    - **J2000** -- `DDDDDDDD.ddd`
"""
function Epoch(s::AbstractString; origin::O=J2000) where {O<:AbstractEpochOrigin}
    # ISO 
    m = match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{1,}", s)
    if length(m.match) != 0 
        sub = split(s, " ")
        if length(sub) == 1
            return Epoch(m.match, eval(Symbol(sub[2])); origin=origin)
        else 
            return Epoch(m.match, TDB; origin=origin)
        end
    end 
    # JD
    m = match(r"JD", s)
    if length(m.match) != 0 
        sub = split(s, " ")
        return Epoch(parse(Float64, sub[2])*86400.0, TDB; origin=eval(Symbol(sub[1])))
    end
    # J2000
    return Epoch(parse(Float64, s)*86400.0, TDB; origin=origin)
end

function Base.show(io::IO, ep::Epoch) 
    print(io, DateTime(ep), " ", timescale(ep))
end