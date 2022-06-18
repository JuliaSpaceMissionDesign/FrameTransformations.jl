export Instant, TimeUnit,
    milliseconds, seconds, minutes, hours, days, weeks, months, quarters, years, centuries,
    -, *, /, +, value, unit,
    SECONDS_PER_MILLISEC,
    SECONDS_PER_MINUTE,
    SECONDS_PER_HOUR,
    SECONDS_PER_DAY,
    SECONDS_PER_WEEK,
    SECONDS_PER_MONTH, 
    SECONDS_PER_QUARTER,
    SECONDS_PER_YEAR,
    SECONDS_PER_CENTURY

using ..Utils: apply_offset

const SECONDS_PER_MILLISEC = 1.0/1000
const SECONDS_PER_MINUTE   = 60.0
const SECONDS_PER_HOUR     = 60.0 * 60.0
const SECONDS_PER_DAY      = 60.0 * 60.0 * 24.0
const SECONDS_PER_WEEK     = 60.0 * 60.0 * 24.0 * 7.0
const SECONDS_PER_MONTH    = 60.0 * 60.0 * 24.0 * 30.0
const SECONDS_PER_QUARTER  = 60.0 * 60.0 * 24.0 * 30.0 * 4.0
const SECONDS_PER_YEAR     = 60.0 * 60.0 * 24.0 * 365.25
const SECONDS_PER_CENTURY  = 60.0 * 60.0 * 24.0 * 365.25 * 100.0

abstract type TimeUnit end

Base.broadcastable(u::TimeUnit) = Ref(u)

struct Millisecond <: TimeUnit end 
struct Second <: TimeUnit end
struct Minute <: TimeUnit end
struct Hour <: TimeUnit end
struct Day <: TimeUnit end
struct Week <: TimeUnit end 
struct Month <: TimeUnit end 
struct Quarter <: TimeUnit end 
struct Year <: TimeUnit end
struct Century <: TimeUnit end

const milliseconds = Millisecond()
const seconds = Second()
const minutes = Minute()
const hours = Hour()
const days = Day()
const weeks = Week()
const months = Month()
const quarters = Quarter()
const years = Year()
const centuries = Century()

# Time units expressed in seconds
factor(::Millisecond) = SECONDS_PER_MILLISEC
factor(::Second) = 1.0
factor(::Minute) = SECONDS_PER_MINUTE
factor(::Hour) = SECONDS_PER_HOUR
factor(::Day) = SECONDS_PER_DAY
factor(::Week) = SECONDS_PER_WEEK
factor(::Month) = SECONDS_PER_MONTH
factor(::Quarter) = SECONDS_PER_QUARTER
factor(::Year) = SECONDS_PER_YEAR
factor(::Century) = SECONDS_PER_CENTURY

struct Instant{U<:TimeUnit, T}
    unit::U
    value::Int64
    fraction::T
    error::T
end

function Instant(unit, δt)
    isfinite(δt) || throw(error("Impossible to initialize a non-finite `Instant`."))

    values = δt * factor(unit)
    int_values = floor(Int64, values)
    fraction = values - int_values
    return Instant(unit, int_values, fraction, zero(fraction))
end

"""
    (u::TimeUnit)(inst::Instant)
Convert instant to new TimeUnit `u`
"""
(u::TimeUnit)(inst::Instant) = Instant(u, inst.value, inst.fraction, inst.error)

"""
    unit(inst::Instant)
Return the unit of the Instant `inst`.
"""
unit(inst::Instant) = inst.unit

"""
    value(inst::Instant)
Return the unitless value of the Instant `inst`.
"""
value(inst::Instant) = (inst.fraction + inst.value) / factor(unit(inst))

name(::Millisecond) = "milliseconds"
name(::Second) = "seconds"
name(::Minute) = "minutes"
name(::Hour) = "hours"
name(::Day) = "days"
name(::Week) = "weeks"
name(::Month) = "months"
name(::Quarter) = "quarters"
name(::Year) = "years"
name(::Century) = "centuries"

# Initialize empty Instants
Base.zero(inst::Instant) = Instant(unit(inst), zero(value(inst)))
Base.zero(inst::Type{<:Instant{U}}) where {U} = Instant(U(), 0.0)
Base.zero(inst::Type{<:Instant{U,T}}) where {U, T} = Instant(U(), zero(T))
Base.eltype(inst::Instant) = typeof(value(inst))
Base.eltype(inst::Type{<:Instant{U,T}}) where {U, T} = T

# Default representation of Instants
function Base.show(io::IO, inst::Instant)
    u = unit(inst)
    v = value(inst)
    print(io, v, " ", name(u))
end

# ----
# Operations

Base.:*(dt::Number, unit::TimeUnit) = Instant(unit, dt)
Base.:*(unit::TimeUnit, dt::Number) = Instant(unit, dt)
Base.:*(A::TimeUnit, B::AbstractArray) = broadcast(*, A, B)
Base.:*(A::AbstractArray, B::TimeUnit) = broadcast(*, A, B)

Base.:-(inst::Instant) = Instant(unit(inst), -inst.value, -inst.fraction, -inst.error)

function Base.:+(inst1::Instant{U}, inst2::Instant{U}) where U
    value, fraction, error = apply_offset(inst1.value, inst1.fraction, inst1.error, inst2.value, inst2.fraction, inst2.error)
    return Instant(U(), value, fraction, error)
end

Base.:-(inst1::Instant, inst2::Instant) = inst1 + (-inst2)
Base.:*(x, inst::Instant) = Instant(unit(inst), value(inst) * x)
Base.:*(inst::Instant, x) = Instant(unit(inst), value(inst) * x)
Base.:/(inst::Instant, x) = Instant(unit(inst), value(inst) / x)

Base.isless(inst1::Instant{U}, inst2::Instant{U}) where {U} = isless(value(inst1), value(inst2))
Base.:(==)(inst1::Instant{U}, inst2::Instant{U}) where {U} = value(inst1) == value(inst2)
function Base.isapprox(inst1::Instant{U}, inst2::Instant{U}; kwargs...) where {U}
    return isapprox(value(inst1), value(inst2); kwargs...)
end

(::Base.Colon)(start::Instant{U,T}, stop::Instant{U,T}) where {U,T} = (:)(start, one(T) * U(), stop)

function (::Base.Colon)(start::Instant{U}, step::Instant{U}, stop::Instant{U}) where {U}
    step = start < stop ? step : -step
    StepRangeLen(start, step, floor(Int, value(stop-start)/value(step))+1)
end

Instant{U,T}(p::Instant{U,T}) where {U,T} = p

Base.step(r::StepRangeLen{<:Instant}) = r.step