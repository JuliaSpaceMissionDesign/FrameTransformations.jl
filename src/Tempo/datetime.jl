export Date, Time, DateTime,
    year, month, day, calendar, isleapyear,
    hour, minute, second, millisecond, microsecond, nanosecond,
    finddayinyear, findfractionofday, findfractionofsecond, findsecondinday,
    j2000, j2000seconds, offset, offsets,
    J2000_DATE, JULIAN_DATE, MODIFIED_JULIAN_DATE, FIFTIES_DATE, CCSDS_DATE,
    GALILEO_DATE, GPS_DATE, UNIX_DATE, 
    MIDNIGHT, NOON

function last_j2000_dayofyear(calendar::Symbol, year::Int)
    if calendar == :proleptic_julian
        return 365 * year + (year + 1) ÷ 4 - 730123
    elseif calendar == :julian
        return 365 * year + year ÷ 4 - 730122
    end

    return 365 * year + year ÷ 4 - year ÷ 100 + year ÷ 400 - 730120
end

function findyear(calendar::Symbol, j2000day::Int)
    j2kday = ifelse(j2000day isa Int32, widen(j2000day), j2000day)
    if calendar == :proleptic_julian
        return -((-4 * j2kday - 2920488) ÷ 1461)
    elseif calendar == :julian
        return -((-4 * j2kday - 2921948) ÷ 1461)
    end

    year = (400 * j2kday + 292194288) ÷ 146097

    # The previous estimate is one unit too high in some rare cases
    # (240 days in the 400 years gregorian cycle, about 0.16%)
    if j2kday <= last_j2000_dayofyear(:gregorian, year - 1)
        year -= 1
    end

    return year
end

function isleapyear(calendar::Symbol, year::Int)
    if calendar in (:proleptic_julian, :julian)
        return year % 4 == 0
    end

    return year % 4 == 0 && (year % 400 == 0 || year % 100 != 0)
end

const PREVIOUS_MONTH_END_DAY_LEAP = (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
const PREVIOUS_MONTH_END_DAY      = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

function findmonth(dayinyear::Int, isleap::Bool)
    offset = ifelse(isleap, 313, 323)
    return ifelse(dayinyear < 32, 1, (10 * dayinyear + offset) ÷ 306)
end

function findday(dayinyear::Int, month::Int, isleap::Bool)
    (!isleap && dayinyear > 365) &&
        throw(ArgumentError("Day of year cannot be 366 for a non-leap year."))
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP, PREVIOUS_MONTH_END_DAY)
    return dayinyear - previous_days[month]
end

function finddayinyear(month::Int, day::Int, isleap::Bool)
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP, PREVIOUS_MONTH_END_DAY)
    return day + previous_days[month]
end

function findcalendar(year::Int64, month::Int64, day::Int64)
    if year < 1583
        if year < 1
            return :proleptic_julian
        elseif year < 1582 || month < 10 || (month < 11 && day < 5)
            return :julian
        end
    end
    return :gregorian
end

function j2000(calendar::Symbol, year::Int, month::Int, day::Int)
    d1 = last_j2000_dayofyear(calendar, year - 1)
    d2 = finddayinyear(month, day, isleapyear(calendar, year))
    return d1 + d2
end

function j2000(year::Int, month::Int, day::Int)
    calendar = findcalendar(year, month, day)
    return j2000(calendar, year, month, day)
end

struct Date
    year::Int 
    month::Int 
    day::Int 
    calendar::Symbol
end 

year(d::Date) = d.year
month(d::Date) = d.month
day(d::Date) = d.day
calendar(d::Date) = d.calendar
isleapyear(d::Date) = isleapyear(findcalendar(year(d), month(d), day(d)), year(d))
finddayinyear(d::Date) = finddayinyear(month(d), day(d), isleapyear(d))
j2000(d::Date) = j2000(calendar(d), year(d), month(d), day(d))

function Date(offset::Integer)
    calendar = :gregorian
    if offset < -152384
        if offset > -730122
            calendar = :julian
        else
            calendar = :proleptic_julian
        end
    end

    year = findyear(calendar, offset)
    dayinyear = offset - last_j2000_dayofyear(calendar, year - 1)

    month = findmonth(dayinyear, isleapyear(calendar, year))
    day = findday(dayinyear, month, isleapyear(calendar, year))

    return Date(year, month, day, calendar)
end

Date(d::Date, offset::Integer) = Date(j2000(d) + offset)

function Date(year::Integer, month::Integer, day::Integer)
    if month < 1 || month > 12
        throw(ArgumentError("Invalid month number: $month"))
    end

    check = Date(j2000(year, month, day))
    if check.year != year || check.month != month || check.day != day
        throw(ArgumentError("Invalid date."))
    end
    return Date(year, month, day, findcalendar(year, month, day))
end

function Date(year::Integer, dayinyear::Integer)
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
    leap = isleapyear(calendar, year)
    month = findmonth(dayinyear, leap)
    day = findday(dayinyear, month, leap)

    return Date(year, month, day, calendar)
end

function Base.show(io::IO, d::Date)
    return print(io, year(d), "-", lpad(month(d), 2, '0'), "-", lpad(day(d), 2, '0'))
end

function Base.isapprox(a::Date, b::Date; kwargs...)
    return a.year == b.year &&
           a.month == b.month &&
           a.day == b.day
end

Base.:-(d::Date) = -j2000(d)
Base.:-(d1::Date, d2::Date) = (j2000(d1) + (-d2))days
Base.:+(d::Date, x::Int) = Date(d, x)
Base.:-(d::Date, x::Int) = Date(d, -x)
function Base.:+(d::Date, i::Instant) 
    ref = days(i)
    ref.fraction != 0 || mod(value(ref), 1) != 0 && throw(
        error("It is not possible add a non-integer (day) `Instant` to `Date`! Use `DateTime` type."))
    return d + floor(Int64, value(ref))
end
function Base.:-(d::Date, i::Instant) 
    ref = days(i)
    ref.fraction != 0 || mod(value(ref), 1) != 0 && throw(
        error("It is not possible add a non-integer (day) `Instant` to `Date`! Use `DateTime` type."))
    return d - floor(Int64, value(ref))
end

struct Time{T}
    hour::Int
    minute::Int
    second::Int 
    fraction::T 
    function Time(hour, minute, second, fraction::T) where {T}
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

function Time(hour::Int, minute::Int, second::Number)
    sec, frac = divrem(second, 1)
    return Time(hour, minute, sec, frac)
end

function Time(secondinday::Int, fraction::Number)
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
hour(t::Time) = t.hour
minute(t::Time) = t.minute
second(::Type{Float64}, t::Time) = t.fraction + t.second
second(::Type{Int}, t::Time) = t.second
second(t::Time) = second(Int, t)
millisecond(t::Time) = subsecond(t.fraction, 3)
microsecond(t::Time) = subsecond(t.fraction, 6)
nanosecond(t::Time) = subsecond(t.fraction, 9)

findfractionofday(t::Time) = (t.fraction + t.second) / 86400 + t.minute / 1440 + t.hour / 24
findfractionofsecond(t::Time) = t.fraction
findsecondinday(t::Time) = t.fraction + t.second + 60 * t.minute + 3600 * t.hour

function Base.show(io::IO, t::Time)
    h = lpad(hour(t), 2, '0')
    m = lpad(minute(t), 2, '0')
    s = lpad(second(t), 2, '0')
    f = lpad(findfractionofsecond(t), 9, '0')
    return print(io, h, ":", m, ":", s, ".", f[3:6])
end

struct DateTime{T} <: AbstractDateTimeEpoch
    date::Date
    time::Time{T}
end

function DateTime(year::Int, month::Int=1, day::Int=1, hour::Int=0, min::Int=0, sec::Int=0, frac::T=0.0) where T
    return DateTime(Date(year, month, day), Time(hour, min, sec, frac))
end

function DateTime(s::AbstractString)
    length(split(s)) != 1  && throw(error("Cannot parse $s as `DateTime` as it has a `TimeScale`! Please use `Epoch` instead."))
    dy, dm, dd, th, tm, ts, tms = parse_iso(s)
    return DateTime(dy, dm, dd, th, tm, ts, tms)
end

function sec2hms(s::Float64)
    mins, sec_ = divrem(s, 60)
    sec = floor(Int64, sec_)
    frac = sec_ - sec
    hr, min = divrem(mins, 60)
    return abs(floor(Int64, hr)), floor(Int64, min), sec, frac
end

# !!!! seconds since `d` at midnight
function DateTime(d::Date, seconds::Float64)
    nday, remsec = divrem(seconds, SECONDS_PER_DAY)
    hrs, min, sec, frac = sec2hms(remsec)
    DateTime(Date(d, floor(Int64, nday)), Time(hrs, min, sec, frac))
end

# !!!! seconds since 2000-01-01T12:00:00.0000
function DateTime(seconds::Float64)
    nday, remsec = divrem(seconds+SECONDS_PER_DAY/2, SECONDS_PER_DAY)
    hrs, min, sec, frac = sec2hms(remsec)
    DateTime(Date(floor(Int, nday)), Time(hrs, min, sec, frac))
end 

# !!!! referred true J2000 Epoch (2000-01-01T12:00:00.0000)
j2000(dt::DateTime) = j2000(Date(dt)) + findfractionofday(Time(dt)) - 0.5
j2000seconds(dt::DateTime) = j2000(dt) * SECONDS_PER_DAY

function Base.isless(d1::DateTime, d2::DateTime)
    return j2000(d1) < j2000(d2)
end

function Base.:(==)(d1::DateTime, d2::DateTime)
    return j2000(d1) == j2000(d2)
end

function Base.isapprox(d1::DateTime, d2::DateTime; kwargs...) where {U}
    return isapprox(j2000(d1), j2000(d2); kwargs...)
end

(::Base.Colon)(start::DateTime, stop::DateTime) = (:)(start, 1days, stop)
function (::Base.Colon)(start::DateTime, step::Instant{U}, stop::DateTime) where {U}
    step = start < stop ? step : -step
    StepRangeLen(start, step, floor(Int, value(seconds(stop-start))/value(seconds(step)))+1)
end

DateTime{T}(dt::DateTime) where {T} = dt

function Base.:-(d1::DateTime, d2::DateTime)
    Δday = (j2000(Date(d1)) + (-Date(d2)))days 
    Δt = (findsecondinday(Time(d1)) - findsecondinday(Time(d2)))seconds 
    return seconds(Δday) + Δt
end

function Base.:+(d1::DateTime, x::Number)
    Δt = findsecondinday(Time(d1))
    return DateTime(Date(d1), x + Δt)
end 

function Base.:-(d1::DateTime, x::Number)
    Δt = findsecondinday(Time(d1))
    return DateTime(Date(d1), Δt - x)
end 

function Base.:-(d1::DateTime, i::Instant)
    Δt = findsecondinday(Time(d1))
    return DateTime(Date(d1), Δt - value(seconds(i)))
end 

function Base.:+(d1::DateTime, i::Instant)
    Δt = findsecondinday(Time(d1))
    return DateTime(Date(d1), Δt + value(seconds(i)))
end 

Date(dt::DateTime) = dt.date
Time(dt::DateTime) = dt.time

Base.show(io::IO, dt::DateTime) = print(io, Date(dt), "T", Time(dt))

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

offset(ref::Date, target::Date) = (j2000(ref) - j2000(target))days # results in days
offset(ref::Time, target::Time) = (findsecondinday(ref) - findsecondinday(target))seconds
function offsets(ref::DateTime, target::DateTime)
    offset(ref.date, target.ref), offset(ref.time, target.time)
end
