export Date, Time, DateTime,
       year, day, calendar, month, day, find_dayinyear, is_leapyear,
       hour, minute, second, millisecond, microsecond, nanosecond, 
       find_fractionofday, find_fractionofsecond, find_secondinday,
       J2000_DATE, JULIAN_DATE, MODIFIED_JULIAN_DATE, FIFTIES_DATE, CCSDS_DATE,
       GALILEO_DATE, GPS_DATE, UNIX_DATE, MIDNIGHT, NOON


const PREVIOUS_MONTH_END_DAY_LEAP = (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
const PREVIOUS_MONTH_END_DAY      = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

function last_j2000_day_of_year(calendar::Symbol, year::N) where {N<:Integer}
    if calendar == :proleptic_julian
        return 365 * year + (year + 1) ÷ 4 - 730123
    elseif calendar == :julian
        return 365 * year + year ÷ 4 - 730122
    end

    return 365 * year + year ÷ 4 - year ÷ 100 + year ÷ 400 - 730120
end

function find_calendar(year::N, month::N, day::N)::Symbol where {N<:Integer}
    if year < 1583
        if year < 1
            return :proleptic_julian
        elseif year < 1582 || month < 10 || (month < 11 && day < 5)
            return :julian
        end
    end
    return :gregorian
end

function find_year(calendar::Symbol, j2000day::N) where {N<:Integer}
    j2kday = ifelse(j2000day isa Int32, widen(j2000day), j2000day)
    if calendar == :proleptic_julian
        return -((-4 * j2kday - 2920488) ÷ 1461)
    elseif calendar == :julian
        return -((-4 * j2kday - 2921948) ÷ 1461)
    end

    year = (400 * j2kday + 292194288) ÷ 146097

    # The previous estimate is one unit too high in some rare cases
    # (240 days in the 400 years gregorian cycle, about 0.16%)
    if j2kday <= last_j2000_day_of_year(:gregorian, year - 1)
        year -= 1
    end

    return year
end

function find_month(dayinyear::N, isleap::Bool) where {N<:Integer}
    offset = ifelse(isleap, 313, 323)
    return ifelse(dayinyear < 32, 1, (10 * dayinyear + offset) ÷ 306)
end

function find_day(dayinyear::N, month::N, isleap::Bool) where {N<:Integer}
    (!isleap && dayinyear > 365) &&
        throw(ArgumentError("Day of year cannot be 366 for a non-leap year."))
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP, PREVIOUS_MONTH_END_DAY)
    return dayinyear - previous_days[month]
end

function find_dayinyear(month::N, day::N, isleap::Bool) where {N<:Integer}
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP, PREVIOUS_MONTH_END_DAY)
    return day + previous_days[month]
end

function is_leapyear(calendar::Symbol, year::N) where {N<:Integer}
    if calendar in (:proleptic_julian, :julian)
        return year % 4 == 0
    end

    return year % 4 == 0 && (year % 400 == 0 || year % 100 != 0)
end

function j2000(calendar::Symbol, year::N, month::N, day::N) where {N<:Integer}
    d1 = last_j2000_day_of_year(calendar, year - 1)
    d2 = find_dayinyear(month, day, is_leapyear(calendar, year))
    return d1 + d2
end

function j2000(year::N, month::N, day::N) where {N<:Integer}
    calendar = find_calendar(year, month, day)
    return j2000(calendar, year, month, day)
end

# ----
# Date


"""
    Date{N<:Integer}

Type to represent a calendar date.

### Fields

- `year` -- year
- `month` -- month 
- `day` -- day
- `calendar` -- calendar used to represent the date. Can be `:gregorian`, 
                `:julian` or `:proleptic_julian`

### Constructors 

- `Date{N}(year::N, month::N, day::N, calendar::Symbol)` -- is the default constructor.

- `Date(year::N, month::N, day::N) where {N<:Integer}` -- computes automatically 
   the calendar.

- `Date(offset::N) where {N<:Integer}` -- initialize from **integer** day 
   offset from `2000-01-01`.

- `Date(d::Date, offset::N) where {N<:Integer}` -- day offset from `d`.

- `Date(year::N, dayinyear::N) where {N<:Integer}` -- initialize giving the 
   year and the day of the year.

- `Date(dt::DateTime)` -- extract date from [`DateTime`](@ref) objects.
"""
struct Date{N<:Integer}
    year::N 
    month::N 
    day::N 
    calendar::Symbol
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
    calendar(d::Date)

Get calendar associated to `Date` type.
"""
calendar(d::Date) = d.calendar

"""
    is_leapyear(d::Date)::Bool

Find if `Date` has a leap year.
"""
is_leapyear(d::Date) = is_leapyear(find_calendar(year(d), month(d), day(d)), year(d))

"""
    find_dayinyear(d::Date)

Find day in the year.
"""
find_dayinyear(d::Date) = find_dayinyear(month(d), day(d), is_leapyear(d))

"""
    j2000(d::Date)

Convert `Date` to days since J2000.

### Example 
```julia-repl
julia> d = Date(65)
2000-03-06
julia> Tempo.j2000(d)
65
```
"""
j2000(d::Date) = j2000(calendar(d), year(d), month(d), day(d))


function Date(year::N, month::N, day::N) where {N<:Integer}
    if month < 1 || month > 12
        throw(ArgumentError("Invalid month number: $month"))
    end

    check = Date(j2000(year, month, day))
    if check.year != year || check.month != month || check.day != day
        throw(ArgumentError("Invalid date."))
    end
    return Date(year, month, day, find_calendar(year, month, day))
end

function Date(offset::N) where {N<:Integer}
    calendar = :gregorian
    if offset < -152384
        if offset > -730122
            calendar = :julian
        else
            calendar = :proleptic_julian
        end
    end
    year = find_year(calendar, offset)
    dayinyear = offset - last_j2000_day_of_year(calendar, year - 1)
    month = find_month(dayinyear, is_leapyear(calendar, year))
    day = find_day(dayinyear, month, is_leapyear(calendar, year))
    return Date(year, month, day, calendar)
end

function Date(year::N, dayinyear::N) where {N<:Integer}
    calendar = :gregorian
    if year < 1583
        if year < 1
            calendar = :proleptic_julian
        else
            calendar = :julian
        end
    end
    if dayinyear <= 0
        throw(error("Day in year must me ≥ than 0! $dayinyear provided.")) 
    end
    leap = is_leapyear(calendar, year)
    month = find_month(dayinyear, leap)
    day = find_day(dayinyear, month, leap)

    return Date(year, month, day, calendar)
end

Date(d::Date, offset::N) where {N<:Integer} = Date(j2000(d) + offset)

function Base.show(io::IO, d::Date)
    return print(io, year(d), "-", lpad(month(d), 2, '0'), "-", lpad(day(d), 2, '0'))
end

function Base.isapprox(a::Date, b::Date; kwargs...)
    return a.year == b.year &&
           a.month == b.month &&
           a.day == b.day
end

Base.:+(d::Date, x::N) where {N<:Integer} = Date(d, x)
Base.:-(d::Date, x::N) where {N<:Integer} = Date(d, -x)

# ----
# Time 
"""
    Time{N<:Integer, T<:AbstractFloat}

A type representing a time of the day.

### Fields

- `hour` -- hour 
- `minute` -- minute 
- `second` -- second 
- `fraction` -- fraction of seconds

### Constructors

- `Time(hour::N, minute::N, second::N, 
    fraction::T) where {N<:Integer, T<:AbstractFloat}` -- default constructor. 

- `Time(hour::N, minute::N, second::T) where {N<:Integer, T<:AbstractFloat}` -- 
    initialize a `Time` type automatically computing seconds fractions 

- `Time(secondinday::Integer, fraction::T) where {T<:AbstractFloat}` -- type from 
    second in the day and fraction of seconds. 

- `Time(dt::DateTime)` -- extract time from [`DateTime`](@ref) objects.
"""
struct Time{N, T}
    hour::N
    minute::N
    second::N 
    fraction::T 
    function Time(hour::N, minute::N, 
        second::N, fraction::T) where {N<:Integer, T<:AbstractFloat}
        if hour < 0 || hour > 23
            throw(ArgumentError("`hour` must be an integer between 0 and 23."))
        elseif minute < 0 || minute > 59
            throw(ArgumentError("`minute` must be an integer between 0 and 59."))
        elseif second < 0 || second >= 61
            throw(ArgumentError("`second` must be an integer between 0 and 61."))
        elseif fraction < 0 || fraction > 1
            throw(ArgumentError("`fraction` must be a number between 0 and 1."))
        end
        return new{N, T}(hour, minute, second, fraction)
    end
end

function Time(hour::N, minute::N, second::T) where {N<:Integer, T<:AbstractFloat}
    sec, frac = divrem(second, 1)
    return Time(hour, minute, convert(N, sec), frac)
end

function Time(secondinday::Integer, fraction::T) where {T<:AbstractFloat}
    if secondinday < 0 || secondinday > 86400
        throw(ArgumentError("Seconds are out of range. Must be between 0 and 86400, provided $secondinday."))
    end
    hour = secondinday ÷ 3600
    secondinday -= 3600 * hour
    minute = secondinday ÷ 60
    secondinday -= 60 * minute
    return Time(hour, minute, secondinday, fraction)
end

function subsecond(fraction, n, r)
    n % 3 == 0 || throw(ArgumentError("`n` must be divisible by 3."))
    factor = ifelse(Int === Int32, widen(10), 10)^n
    rounded = round(fraction, r; digits=n)
    return round(Int64, rounded * factor, r) % 1000
end

function subsecond(fraction, n)
    r = ifelse(subsecond(fraction, n+3, RoundNearest) == 0, RoundNearest, RoundToZero)
    return subsecond(fraction, n, r)
end


subsecond(t::Time, n) = subsecond(t.fraction, n)

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

"""
    find_fractionofday(t::Time)::AbstractFloat

Find fraction of day.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.find_fractionofday(t)
0.5213002592592593  # days
```
"""
find_fractionofday(t::Time) = (t.fraction + t.second) / 86400 + t.minute / 1440 + t.hour / 24

"""
    find_fractionofsecond(t::Time)::AbstractFloat

Find fraction of seconds.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.find_fractionofsecond(t)
0.3423999999999978  # seconds
```
"""
find_fractionofsecond(t::Time) = t.fraction

"""
    find_secondinday(t::Time)::AbstractFloat

Find second in the day.

### Example

```julia-repl
julia> t = Time(12, 30, 40.3424)
12:30:40.3423
julia> Tempo.find_secondinday(t)
45040.3424  # seconds
```
"""
find_secondinday(t::Time) = t.fraction + t.second + 60 * t.minute + 3600 * t.hour

function Base.show(io::IO, t::Time)
    h = lpad(hour(t), 2, '0')
    m = lpad(minute(t), 2, '0')
    s = lpad(second(t), 2, '0')
    f = lpad(find_fractionofsecond(t), 9, '0')
    return print(io, h, ":", m, ":", s, ".", f[3:6])
end


"""
    DateTime{N<:Integer, T<:AbstractFloat} <: AbstractDateTimeEpoch

A type wrapping a date and a time since a reference date.

### Fields

- `date` -- `Date` part of the type
- `time` -- `Time` part of the type

### Constructors 

- `DateTime{N, T}(date::Date{N}, time::Time{N, T})` -- default constructor.

- `DateTime(year::N, month::N, day::N, hour::N, min::N, 
    sec::N, frac::T=0.0) where {N<:Integer, T<:AbstractFloat}` -- full constructor

- `DateTime(s::AbstractString)` -- parse `DateTime` from ISO string using 
    [`parse_iso`](@ref).

- `DateTime(seconds::T) where {T<:AbstractFloat}` -- parse from seconds since J2000.

- `DateTime(d::Date{N}, sec::T) where {T<:AbstractFloat, N<:Integer}` -- parse 
    as seconds since `d` at midnight.

- `DateTime{N, T}(dt::DateTime) where {N, T}` -- ghost constructor
"""
struct DateTime{N<:Integer, T<:AbstractFloat} <: AbstractDateTimeEpoch
    date::Date{N}
    time::Time{N, T}
end

function DateTime(year::N, month::N, day::N, hour::N, min::N, 
    sec::N, frac::T=0.0) where {N<:Integer, T<:AbstractFloat}
    return DateTime(Date(year, month, day), Time(hour, min, sec, frac))
end

function DateTime(s::AbstractString)
    length(split(s)) != 1  && throw(error("Cannot parse $s as `DateTime`."))
    dy, dm, dd, th, tm, ts, tms = parse_iso(s)
    return DateTime(dy, dm, dd, th, tm, ts, tms)
end

function sec2hms(s::T) where {T<:AbstractFloat}
    mins, sec_ = divrem(s, 60)
    sec = floor(Int, sec_)
    frac = sec_ - sec
    hr, min = divrem(mins, 60)
    return abs(floor(Int, hr)), floor(Int, min), sec, frac
end

function DateTime(seconds::T) where {T<:AbstractFloat}
    if seconds < 0.0
        throw(error("[Tempo/DateTime] CANNOT PARSE DateTime WITH NEGATIVE SECONDS!"))
    end

    nday, remsec = divrem(seconds+86400.0/2, 86400.0)
    hrs, min, sec, frac = sec2hms(remsec)
    DateTime(Date(floor(Int, nday)), Time(hrs, min, sec, frac))
end 

function DateTime(d::Date{N}, sec::T) where {T<:AbstractFloat, N<:Integer}
    nday, remsec = divrem(sec, 86400.0)
    h, m, s, frac = sec2hms(remsec)
    DateTime(Date(d, floor(N, nday)), Time(h, m, s, frac))
end

Date(dt::DateTime) = dt.date
Time(dt::DateTime) = dt.time
DateTime{N, T}(dt::DateTime) where {N, T} = dt

Base.show(io::IO, dt::DateTime) = print(io, Date(dt), "T", Time(dt))

"""
    j2000(dt::DateTime)

Find days since J2000 (i.e. `2000-01-01T12:00:00.0000`).
"""
j2000(dt::DateTime) = j2000(Date(dt)) + find_fractionofday(Time(dt)) - 0.5

"""
    j2000seconds(dt::DateTime)

Find seconds since J2000 (i.e. `2000-01-01T12:00:00.0000`).
"""
j2000seconds(dt::DateTime) = j2000(dt) * 86400


function Base.isless(d1::DateTime, d2::DateTime)
    return j2000(d1) < j2000(d2)
end

function Base.:(==)(d1::DateTime, d2::DateTime)
    return j2000(d1) == j2000(d2)
end

function Base.isapprox(d1::DateTime, d2::DateTime; kwargs...) where {U}
    return isapprox(j2000(d1), j2000(d2); kwargs...)
end

function Base.:+(d1::DateTime, x::N) where {N<:Number}
    Δt = find_secondinday(Time(d1))
    return DateTime(Date(d1), x + Δt)
end 

function Base.:-(d1::DateTime, x::N) where {N<:Number}
    Δt = find_secondinday(Time(d1))
    return DateTime(Date(d1), Δt - x)
end 

const JULIAN_DATE = Date(-4712, 1, 1)
const MODIFIED_JULIAN_DATE = Date(1858, 11, 17)
const FIFTIES_DATE = Date(1950, 1, 1)
const CCSDS_DATE = Date(1958, 1, 1)
const GALILEO_DATE = Date(1999, 8, 22)
const GPS_DATE = Date(1980, 1, 6)
const J2000_DATE = Date(2000, 1, 1)
const UNIX_DATE = Date(1970, 1, 1)

const MIDNIGHT = Time(0, 0, 0, 0.0)
const NOON = Time(12, 0, 0, 0.0)