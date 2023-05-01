export Date,
    Time,
    year,
    month,
    day,
    isleapyear,
    find_dayinyear,
    cal2jd,
    j2000,
    j2000s,
    j2000c,
    hour,
    minute,
    second,
    fraction_of_second,
    fraction_of_day,
    second_in_day,
    DateTime

function lastj2000dayofyear(year::N) where {N<:Integer}
    return 365 * year + year ÷ 4 - year ÷ 100 + year ÷ 400 - 730120
end

"""
    find_year(d::N) where {N<:Integer}

Find year from j2000 day
"""
function find_year(d::N) where {N<:Integer}
    j2d = ifelse(d isa Int32, widen(d), d)
    year = (400 * j2d + 292194288) ÷ 146097
    # The previous estimate is one unit too high in some rare cases
    # (240 days in the 400 years gregorian cycle, about 0.16%)
    if j2d <= lastj2000dayofyear(year - 1)
        year -= 1
    end
    return year
end

"""
    find_month(dayinyear::N, isleap::Bool) where {N<:Integer}

Find month from the day in the year and using if the year is leap or not.
"""
function find_month(dayinyear::N, isleap::Bool) where {N<:Integer}
    offset = ifelse(isleap, 313, 323)
    return ifelse(dayinyear < 32, 1, (10 * dayinyear + offset) ÷ 306)
end

"""
    find_day(dayinyear::N, month::N, isleap::Bool) where {N<:Integer}

Find day from the day in the year, the month and using if the year is leap or not.
"""
function find_day(dayinyear::N, month::N, isleap::Bool) where {N<:Integer}
    (!isleap && dayinyear > 365) && throw(
        ArgumentError("[Tempo] day of year cannot greater than 366 for a non-leap year."),
    )
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP, PREVIOUS_MONTH_END_DAY)
    return dayinyear - previous_days[month]
end

########
# DATE #
########

"""
    Date{N<:Integer}

Type to represent a calendar date.

### Fields

- `year` -- year
- `month` -- month 
- `day` -- day

### Constructors 

- `Date(year::N, month::N, day::N)` -- is the default constructor.

- `Date(offset::N) where {N<:Integer}` -- initialize from **integer** day 
    offset from `2000-01-01`.

- `Date(d::Date, offset::N) where {N<:Integer}` -- day offset from `d`.

- `Date(year::N, dayinyear::N) where {N<:Integer}` -- initialize giving the 
    year and the day of the year.

- `Date(dt::DateTime)` -- extract date from [`DateTime`](@ref) objects.
"""
struct Date
    year::Int
    month::Int
    day::Int
end

"""
    year(d::Date)

Get year associated to `Date` type.
"""
year(d::Date) = d.year

"""
    month(d::Date)

Get month associated to `Date` type.
"""
month(d::Date) = d.month

"""
    day(d::Date)

Get day associated to `Date` type.
"""
day(d::Date) = d.day

"""
    isleapyear(d::Date)::Bool

Find if `Date` has a leap year.
"""
isleapyear(d::Date) = isleapyear(year(d))

"""
    find_dayinyear(d::Date)

Find day in the year.
"""
find_dayinyear(d::Date) = find_dayinyear(month(d), day(d), isleapyear(d))

"""
    cal2jd(d::Date)

Convert Gregorian calendar date to Julian Date.

### Outputs
- `j2000` -- J2000 zero point: always 2451545
- `d` -- J2000 Date for 12 hrs
"""
cal2jd(d::Date) = cal2jd(year(d), month(d), day(d))

"""
    j2000(d::Date)

Convert Gregorian calendar date Julian Date past J2000
"""
j2000(d::Date) = j2000(cal2jd(d)...)

function Date(offset::N) where {N<:Integer}
    year = find_year(offset)
    dayinyear = offset - lastj2000dayofyear(year - 1)
    ly = isleapyear(year)
    month = find_month(dayinyear, ly)
    day = find_day(dayinyear, month, ly)
    return Date(year, month, day)
end

function Date(year::N, dayinyear::N) where {N<:Integer}
    if dayinyear <= 0
        throw(error("[Tempo] day in year must me ≥ than 0! $dayinyear provided."))
    end
    ly = isleapyear(year)
    month = find_month(dayinyear, ly)
    day = find_day(dayinyear, month, ly)
    return Date(year, month, day)
end

Date(d::Date, offset::N) where {N<:Integer} = Date(convert(N, j2000(d)) + offset)

function Base.show(io::IO, d::Date)
    return print(io, year(d), "-", lpad(month(d), 2, '0'), "-", lpad(day(d), 2, '0'))
end

# Operations 
function Base.isapprox(a::Date, b::Date; kwargs...)
    return a.year == b.year && a.month == b.month && a.day == b.day
end
Base.:+(d::Date, x::N) where {N<:Integer} = Date(d, x)
Base.:-(d::Date, x::N) where {N<:Integer} = Date(d, -x)

########
# TIME #
########

"""
    Time{T<:AbstractFloat}

A type representing a time of the day.

### Fields

- `hour` -- hour 
- `minute` -- minute 
- `second` -- second 
- `fraction` -- fraction of seconds

### Constructors

- `Time{T}(hour::N, minute::N, second::N, 
    fraction::T) where {N<:Integer, T<:AbstractFloat}` -- default constructor. 

- `Time(hour::N, minute::N, second::T) where {N<:Integer, T<:AbstractFloat}` -- 
    initialize a `Time` type automatically computing seconds fractions 

- `Time(secondinday::Integer, fraction::T) where {T<:AbstractFloat}` -- type from 
    second in the day and fraction of seconds. 

- `Time(dt::DateTime)` -- extract time from [`DateTime`](@ref) objects.
"""
struct Time{T}
    hour::Int
    minute::Int
    second::Int
    fraction::T
    function Time(
        hour::N, minute::N, second::N, fraction::T
    ) where {N<:Integer,T<:AbstractFloat}
        if hour < 0 || hour > 23
            throw(ArgumentError("`hour` must be an integer between 0 and 23."))
        elseif minute < 0 || minute > 59
            throw(ArgumentError("`minute` must be an integer between 0 and 59."))
        elseif second < 0 || second >= 61
            throw(ArgumentError("`second` must be an integer between 0 and 61."))
        elseif fraction < 0 || fraction > 1
            throw(ArgumentError("`fraction` must be a number between 0 and 1."))
        end
        return new{T}(hour, minute, second, fraction)
    end
end

function Time(hour::N, minute::N, second::T) where {N<:Integer,T<:AbstractFloat}
    sec, frac = divrem(second, 1)
    return Time(hour, minute, convert(N, sec), frac)
end

function Time(secondinday::Integer, fraction::T) where {T<:AbstractFloat}
    if secondinday < 0 || secondinday > 86400
        throw(
            ArgumentError(
                "[Tempo] seconds are out of range. Must be between 0 and 86400, provided $secondinday.",
            ),
        )
    end
    hour = secondinday ÷ 3600
    secondinday -= 3600 * hour
    minute = secondinday ÷ 60
    secondinday -= 60 * minute
    return Time(hour, minute, secondinday, fraction)
end

"""
    hour(t::Time)::Integer

Get the current hour.
"""
hour(t::Time) = t.hour

"""
    minute(t::Time)::Integer

Get the current minute.
"""
minute(t::Time) = t.minute

"""
    second(::Type{<:AbstractFloat}, t::Time)::AbstractFloat
    second(::Type{<:Integer}, t::Time)::Integer 
    second(t::Time)::Int64 

Get the current second.
"""
second(::Type{<:AbstractFloat}, t::Time) = t.fraction + t.second
second(::Type{<:Integer}, t::Time) = t.second
second(t::Time) = second(Int64, t)

function subsecond(fraction, n, r)
    n % 3 == 0 || throw(ArgumentError("`n` must be divisible by 3."))
    factor = ifelse(Int === Int32, widen(10), 10)^n
    rounded = round(fraction, r; digits=n)
    return round(Int64, rounded * factor, r) % 1000
end

function subsecond(fraction, n)
    r = ifelse(subsecond(fraction, n + 3, RoundNearest) == 0, RoundNearest, RoundToZero)
    return subsecond(fraction, n, r)
end

subsecond(t::Time, n) = subsecond(t.fraction, n)

"""
    millisecond(t::Time)::Integer

Get the current millisecond.
"""
millisecond(t::Time) = subsecond(t.fraction, 3)

"""
    microsecond(t::Time)::Integer

Get the current microsecond.
"""
microsecond(t::Time) = subsecond(t.fraction, 6)

"""
    nanosecond(t::Time)::Integer

Get the current nanosecond.
"""
nanosecond(t::Time) = subsecond(t.fraction, 9)

hms2fd(t::Time) = hms2fd(t.hour, t.minute, t.second + t.fraction)

"""
    fraction_of_day(t::Time)::AbstractFloat
    hms2fd(t::Time)

Find fraction of day.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.fraction_of_day(t)
0.5213002592592593  # days
```
"""
fraction_of_day(t::Time) = hms2fd(t::Time)

"""
    fraction_of_second(t::Time)::AbstractFloat

Find fraction of seconds.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.fraction_of_second(t)
0.3423999999999978  # seconds
```
"""
fraction_of_second(t::Time) = t.fraction

"""
    second_in_day(t::Time)::AbstractFloat

Find second in the day.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.second_in_day(t)
45040.3424  # seconds
```
"""
second_in_day(t::Time) = t.fraction + t.second + 60 * t.minute + 3600 * t.hour

function Base.show(io::IO, t::Time)
    h = lpad(hour(t), 2, '0')
    m = lpad(minute(t), 2, '0')
    s = lpad(second(t), 2, '0')
    f = lpad(millisecond(t), 3, '0')
    return print(io, h, ":", m, ":", s, ".", f)
end

############
# DATETIME #
############

"""
    DateTime{N<:Integer, T<:AbstractFloat} <: AbstractDateTimeEpoch

A type wrapping a date and a time since a reference date.

### Fields

- `date` -- `Date` part of the type
- `time` -- `Time` part of the type

### Constructors 

- `DateTime{T}(date::Date{N}, time::Time{N, T})` -- default constructor.

- `DateTime(year::N, month::N, day::N, hour::N, min::N, 
    sec::N, frac::T=0.0) where {N<:Integer, T<:AbstractFloat}` -- full constructor

- `DateTime(s::AbstractString)` -- parse `DateTime` from ISO string using 
    [`parse_iso`](@ref).

- `DateTime(seconds::T) where {T<:AbstractFloat}` -- parse from seconds since J2000.

- `DateTime(d::Date, sec::T) where {T<:AbstractFloat}` -- parse as seconds since `d`.

- `DateTime{N, T}(dt::DateTime) where {N, T}` -- ghost constructor

- `DateTime(e::Epoch)` -- construct from `Epoch`
"""
struct DateTime{T<:AbstractFloat}
    date::Date
    time::Time{T}
end

function DateTime(
    year::N, month::N, day::N, hour::N, min::N, sec::N, frac::T=0.0
) where {N<:Integer,T<:AbstractFloat}
    return DateTime(Date(year, month, day), Time(hour, min, sec, frac))
end

function DateTime(s::AbstractString)
    length(split(s)) != 1 && throw(error("[Tempo] cannot parse $s as `DateTime`."))
    dy, dm, dd, th, tm, ts, tms = parse_iso(s)
    return DateTime(dy, dm, dd, th, tm, ts, tms)
end

function DateTime(seconds::T) where {T<:AbstractFloat}
    y, m, d, H, M, Sf = jd2calhms(DJ2000, seconds / DAY2SEC)
    s = floor(Int64, Sf)
    return DateTime(y, m, d, H, M, s, Sf - s)
end

function DateTime(d::Date, sec::T) where {T<:AbstractFloat}
    jd1 = j2000(d) + sec / DAY2SEC
    y, m, d, H, M, Sf = jd2calhms(DJ2000, jd1)
    s = floor(Int64, Sf)
    return DateTime(y, m, d, H, M, s, Sf - s)
end

Date(dt::DateTime) = dt.date
Time(dt::DateTime) = dt.time
DateTime{T}(dt::DateTime{T}) where {T} = dt

"""
    year(d::DateTime)

Get year associated to `DateTime` type.
"""
year(dt::DateTime) = year(Date(dt))

"""
    month(d::DateTime)

Get month associated to `DateTime` type.
"""
month(dt::DateTime) = month(Date(dt))

"""
    day(d::DateTime)

Get day associated to `DateTime` type.
"""
day(dt::DateTime) = day(Date(dt))

"""
    hour(d::DateTime)

Get hour associated to `DateTime` type.
"""
hour(dt::DateTime) = hour(Time(dt))

"""
    minute(d::DateTime)

Get minute associated to `DateTime` type.
"""
minute(dt::DateTime) = minute(Time(dt))

"""
    second(d::DateTime)

Get second associated to `DateTime` type.
"""
second(dt::DateTime) = second(Float64, Time(dt))

Base.show(io::IO, dt::DateTime) = print(io, Date(dt), "T", Time(dt))

"""
    j2000(dt::DateTime)

Convert `DateTime` in Julian days since J2000
"""
function j2000(dt::DateTime)
    jd1, jd2 = calhms2jd(year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt))
    return j2000(jd1, jd2)
end

"""
    j2000s(dt::DateTime)

Convert `DateTime` to seconds since J2000
"""
function j2000s(dt::DateTime)
    return j2000(dt::DateTime) * DAY2SEC
end

"""
    j2000c(dt::DateTime)

Convert `DateTime` in Julian Date since J2000 (centuries)
"""
function j2000c(dt::DateTime)
    return j2000(dt) / CENTURY2DAY
end

Base.isless(d1::DateTime, d2::DateTime) = j2000(d1) < j2000(d2)
Base.:(==)(d1::DateTime, d2::DateTime) = j2000(d1) == j2000(d2)

function Base.isapprox(d1::DateTime, d2::DateTime; kwargs...)
    return isapprox(j2000(d1), j2000(d2); kwargs...)
end

function Base.:+(d1::DateTime, δs::N) where {N<:Number}
    return DateTime(j2000s(d1) + δs)
end

function Base.:-(d1::DateTime, δs::N) where {N<:Number}
    return DateTime(j2000s(d1) - δs)
end
