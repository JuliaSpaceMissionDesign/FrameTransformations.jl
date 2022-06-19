using ..Utils: two_sum, apply_offset

export Epoch, second, timescale, value,
    ISO, J2000, JD, MJD, MJD2000, 
    offset

"""
    EpochFormat
All epoch representation are subtypes of this abstract type.
The following representations are defines:
* [`ISO`](@ref) — Standardized  epoch representation
* [`J2000`](@ref) — Julian Date 2000
* [`JD`](@ref) — Julian Date
* [`MJD`](@ref) — Modified Julian Date
* [`MJD2000`](@ref) — Modified Julian Date 2000
"""
abstract type EpochFormat end 

"""
    ISOEpoch
Standardized epoch representation in the form `YYYY-MM-DDThh:mm:ss.ffffffffff`
"""
struct ISOEpoch <: EpochFormat end 

"""
    JulianDate
Origin for Epoch representation, starts at `-4712-01-01T12:00`.
"""
struct JulianDate <: EpochFormat end

"""
    ModifiedJulianDate
Origin for Epoch representation, starts at `1858-11-17T00:00`.
"""
struct ModifiedJulianDate <: EpochFormat end

"""
    JulianDate2000
Origin for Epoch representation, starts at `2000-01-01T12:00`.
"""
struct JulianDate2000 <: EpochFormat end

"""
    ModifiedJulianDate2000
Origin for Epoch representation, starts at `2000-01-01T00:00`.
"""
struct ModifiedJulianDate2000 <: EpochFormat end

const ISO = ISOEpoch()
const JD = JulianDate()
const MJD = ModifiedJulianDate()
const J2000 = JulianDate2000()
const MJD2000 = ModifiedJulianDate2000()

const J2000_TO_JULIAN = 2.451545e6 * 86400
const J2000_T0_MJD = 51544.5 * 86400
const J2000_TO_MJD2000 = 0.5 * 86400

offset(::ISOEpoch)                  = 0.0
offset(::JulianDate2000)            = 0.0
offset(::ModifiedJulianDate2000)    = J2000_TO_MJD2000
offset(::JulianDate)                = J2000_TO_JULIAN
offset(::ModifiedJulianDate)        = J2000_T0_MJD

struct Epoch{S<:TimeScale, T} <: AbstractDateTimeEpoch
    scale::S 
    second::Int64
    fraction::T
    error::T 
    format::Symbol
    function Epoch{S}(second::Int64, frac::T, err::T=zero(T), fmt::Symbol=:ISO) where {S<:TimeScale, T<:AbstractFloat}
        return new{S, T}(S(), second, frac, err, fmt)
    end
end
Epoch{S,T}(ep::Epoch{S,T}) where {S,T} = ep

second(ep::Epoch) = ep.second
timescale(ep::Epoch) = ep.scale
value(ep::Epoch) = ep.second + ep.fraction

function offsets(js1::Number, js2::Number, origin::EpochFormat)
    js1 -= offset(origin)
    sum, residual = two_sum(js1, js2)
    epoch = floor(Int64, sum)
    off = (sum - epoch) + residual
    return epoch, off
end

function Epoch{S}(js::T; origin::Symbol=:J2000, format::Symbol=:ISO) where {S<:TimeScale, T<:Number}
    ep, off = offsets(js, 0.0, eval(origin))
    return Epoch{S}(ep, off, 0.0, format)
end

function Epoch(js::T, scale::TimeScale; origin::Symbol=:J2000, format::Symbol=:ISO) where {T<:Number}
    Epoch{typeof(scale)}(js; origin=origin, format=format)
end

function Epoch(dt::DateTime, scale::TimeScale; origin::Symbol=:J2000, format::Symbol=:ISO)
    Epoch(j2000seconds(dt), scale; origin=origin, format=format)
end 

function Epoch(s::AbstractString, scale::TimeScale=TDB; origin::Symbol=:J2000, format::Symbol=:ISO)
    # TODO: extend to other string formats
    Epoch(DateTime(s), scale; origin=origin, format=format)
end

function DateTime(ep::Epoch)
    return DateTime(value(ep))
end

# FIXME: not working if origin is JD - something weird with DateTime type to be fixed
function Base.show(io::IO, ep::Epoch) 
    if ep.format == :ISO
        Base.show(io, DateTime(ep))
    elseif ep.format == :DAY 
        print(io, "JD ", value(ep)/SECONDS_PER_DAY)
    elseif ep.format == :SEC 
        print(io, "JD (seconds) ", value(ep))
    else 
        throw(error("Format $(ep.format) Epochs cannot be show in the REPL."))
    end
end 