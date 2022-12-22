const MTAB = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const PREVIOUS_MONTH_END_DAY_LEAP = (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
const PREVIOUS_MONTH_END_DAY      = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

"""
    isleapyear(year::Integer)

Check if a Gregorian year is leap. Return a `Bool`.
"""
function isleapyear(year::Integer)
    return year % 4 == 0 && (year % 400 == 0 || year % 100 != 0)
end

function find_dayinyear(month::N, day::N, isleap::Bool) where {N<:Integer}
    previous_days = ifelse(
        isleap, 
        PREVIOUS_MONTH_END_DAY_LEAP[month], 
        PREVIOUS_MONTH_END_DAY[month]
    )
    return day + previous_days
end

"""
    hms2fd(h::N, m::N, s::T) where {N <: Integer, T <: AbstractFloat}

Convert hours, minutes, seconds to day fraction. The day fraction is returned 
converted in type `T`.
"""
function hms2fd(h::N, m::N, s::T) where {N <: Integer, T <: AbstractFloat}
    # Validate arguments
    if h < 0 || h > 23 
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid hour provided, must be between 0 and 23"
        ))
    elseif m < 0 || m > 59 
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid minutes provided, must be between 0 and 59"
        ))
    elseif s < 0.0 || s >= 60.0 
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid seconds provided, must be between 0.0 and 59.99"
         ))
    end 
    return T(((60.0*(60.0*h + m))+s)/86400.0)
end

"""
    fd2hms(fd::T) where {T<:AbstractFloat}

Convert day fraction to hour, minute, second.
"""
function fd2hms(fd::T) where {T<:AbstractFloat}
    secinday = fd * 86400.0
    if secinday < 0 || secinday > 86400
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "seconds are out of range: must be between 0 and 86400, provided $secinday"
        ))
    end
    hours = Integer(secinday ÷ 3600)
    secinday -= 3600 * hours
    mins = Integer(secinday ÷ 60)
    secinday -= 60 * mins
    return hours, mins, secinday
end

"""
    fd2hmsf(fd::T) where {T<:AbstractFloat}

Convert day fraction to hour, minute, second, fraction of seconds.
"""
function fd2hmsf(fd::T) where {T<:AbstractFloat}
    h, m, sid = fd2hms(fd)
    sec = Integer(sid ÷ 1)
    fsec = sid - sec 
    return h, m, sec, fsec 
end 

"""
    cal2jd(Y::N, M::N, D::N) where {N<:Integer}

Convert Gregorian Calendar to Julian Date.

### Inputs
- `Y, M, D` -- year, month and day in Gregorian calendar

### Outputs
- `j2000` -- J2000 zero point: always 2451545
- `d` -- J2000 Date for 12 hrs

### References
    
- Explanatory Supplement to the Astronomical Almanac,
  P. Kenneth Seidelmann (ed), University Science Books (1992),
  Section 12.92 (p604).

- Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
  Computing, 76, 279-293 (2006), Section 3.

- [ERFA software library](https://github.com/liberfa/erfa/blob/master/src/cal2jd.c)
"""
function cal2jd(Y::N, M::N, D::N) where {N<:Integer}
    # Validate year and month
    if Y < 1583
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid year provided, must be greater than 1583"
        ))
        
    elseif M < 1 || M > 12
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid month provided, must be between 1 and 12"
        ))
        
    end

    # If February in a leap year, 1, otherwise 0
    ly = isleapyear(Y)
    
    # Validate day, taking into account leap years
    if (D < 1) || (D > (MTAB[M] + ly))
        throw(
            EpochConversionError(
                String(Symbol(@__MODULE__)),
                "invalid day provided, shall be between 1 and $(MTAB[M]+ly)"
        ))
    
    end

    Y = Y - 1
    # find j2000 day of the year 
    d1 = 365 * Y + Y ÷ 4 - Y ÷ 100 + Y ÷ 400 - 730120
    # find day in the year
    d2 = find_dayinyear(M, D, ly)
    # compute days since 01-01-2000
    d = d1 + d2
    return convert(N, DJ2000), convert(N, d)
end

"""
    calhms2jd(Y::N, M::N, id::N, h::N, m::N, sec::T) where {N<:Integer, T<:AbstractFloat}

Convert Gregorian Calendar date and time to Julian Date.

### Inputs
- `Y, M, D` -- year, month and day in Gregorian calendar
- `h, m, sec` -- hour, minute and second

### Outputs
- `jd1` -- J2000 zero point: always 2451545.0
- `jd2` -- J2000 Date for 12 hrs
"""
function calhms2jd(Y::N, M::N, id::N, h::N, m::N, 
    sec::T) where {N<:Integer, T<:AbstractFloat}
    jd1, jd2 = cal2jd(Y, M, id)
    fd = hms2fd(h, m, sec)
    return jd1, jd2+fd
end