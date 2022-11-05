"""
    MTAB

Month days table.
"""
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
    previous_days = ifelse(isleap, PREVIOUS_MONTH_END_DAY_LEAP[month], PREVIOUS_MONTH_END_DAY[month])
    return day + previous_days
end

"""
    hms2fd(ihour::N, imin::N, sec::T) where {N <: Integer, T <: AbstractFloat}

Convert hours, minutes, seconds to day fraction. The day fraction is returned 
converted in type `T`.
"""
function hms2fd(ihour::N, imin::N, sec::T) where {N <: Integer, T <: AbstractFloat}
    # Validate arguments and return status.
    if ihour < 0 || ihour > 23 
        throw(error("[Tempo] invalid hour provided, must be between 0 and 23"))
    elseif imin < 0 || imin > 59 
        throw(error("[Tempo] invalid minutes provided, must be between 0 and 59"))
    elseif sec < 0.0 || sec >= 60.0 
        throw(error("[Tempo] invalid seconds provided, must be between 0.0 and 60.0"))
    end 
    return T(((60.0*(60.0*ihour + imin))+sec)/86400.0)
end


"""
    fd2hms(fd::T) where {T<:AbstractFloat}

Convert day fraction to hour, minute, second.
"""
function fd2hms(fd::T) where {T<:AbstractFloat}
    secinday = fd * 86400.0
    if secinday < 0 || secinday > 86400
        throw(ArgumentError("[Tempo] Seconds are out of range. Must be between 0 and 86400, provided $secinday."))
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
    cal2jd(iyear::N, imonth::N, iday::N) where {N<:Integer}

Convert Gregorian Calendar to Julian Date.

### Inputs
- `iyear, imonth, iday` -- year, month and day in Gregorian calendar

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
function cal2jd(iyear::N, imonth::N, iday::N) where {N<:Integer}
    # Validate year and month
    if iyear < 1583
        throw(error("[Tempo] invalid year provided, must be greater than 1583"))
    elseif imonth < 1 || imonth > 12
        throw(error("[Tempo] invalid month provided, must be between 1 and 12"))
    end

    # If February in a leap year, 1, otherwise 0
    ly = isleapyear(iyear)
    
    # Validate day, taking into account leap years
    if (iday < 1) || (iday > (MTAB[imonth] + ly))
        throw(error("[Tempo] invalid day provided, shall be between 1 and $(MTAB[imonth]+ly)"))
    end

    iyear = iyear - 1
    # find j2000 day of the year 
    d1 = 365 * iyear + iyear ÷ 4 - iyear ÷ 100 + iyear ÷ 400 - 730120
    # find day in the year
    d2 = find_dayinyear(imonth, iday, ly)
    # compute days since 01-01-2000 at 12:00
    d = d1 + d2
    return convert(N, DJ2000), convert(N, d)
end

"""
    calhms2jd(iyear::N, imonth::N, id::N, ihour::N, imin::N, sec::T) where {N<:Integer, T<:AbstractFloat}

Convert Gregorian Calendar date and time to Julian Date.

### Inputs
- `iyear, imonth, iday` -- year, month and day in Gregorian calendar
- `ihour, imin, sec` -- hour, minute and second

### Outputs
- `j2000` -- J2000 zero point: always 2451545
- `d` -- J2000 Date for 12 hrs
"""
function calhms2jd(iyear::N, imonth::N, id::N, ihour::N, imin::N, 
    sec::T) where {N<:Integer, T<:AbstractFloat}
    jd1, jd2 = cal2jd(iyear, imonth, id)
    fd = hms2fd(ihour, imin, sec)
    return jd1, jd2+fd-0.5
end

""" 
    jd2cal(dj1::T, dj2::T) where {T<:AbstractFloat}

Julian Date to Gregorian year, month, day, and fraction of a day.

### Inputs

-  `dj1,dj2` -- Julian Date (Notes 1, 2)

### Outputs 

- `iyear::Integer` -- year
- `imonth::Integer` -- month 
- `id::Integer` -- day 
- `fd::AbstractFloat` -- fraction of day 

### Notes 

1. The earliest valid date is 0 (-4713 Jan 1). The largest value accepted is 1e9.

2. The Julian Date is apportioned in any convenient way between the arguments 
   `dj1` and `dj2`. For example, JD=2450123.7 could be expressed in any of these 
   ways, among others:

   | dj1       	| dj2     	|                      	|
   |-----------	|---------	|----------------------	|
   | 2450123.7 	| 0.0     	| (JD method)          	|
   | 2451545.0 	| -1421.3 	| (J2000 method)       	|
   | 2400000.5 	| 50123.2 	| (MJD method)         	|
   | 2450123.5 	| 0.2     	| (date & time method) 	|

### References
    
- Explanatory Supplement to the Astronomical Almanac,
  P. Kenneth Seidelmann (ed), University Science Books (1992),
  Section 12.92 (p604).

- Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
  Computing, 76, 279-293 (2006), Section 3.

- [ERFA software library](https://github.com/liberfa/erfa/blob/master/src/jd2cal.c)
"""
function jd2cal(dj1::N1, dj2::N2) where {N1 <: Real, N2 <: Real}
    dj = dj1 + dj2
    if dj < -68569.5 || dj > 1e9
        throw(error("[Tempo] Invalid date provided, shall be between -68569.5 and 1e9."))
    end

    # Copy the date, big then small, and re-align to midnight
    if abs(dj1) ≥ abs(dj2) 
        d1 = dj1 
        d2 = dj2 
    else
        d1 = dj2 
        d2 = dj1 
    end
    d2 -= 0.5 

    #  Separate day and fraction
    f1 = mod(d1, 1.0)
    f2 = mod(d2, 1.0)
    fd = mod(f1+f2, 1.0)
    if fd < 0.0
        fd += 1.0 
    end
    d = round(Int, d1-f1) + round(Int, d2-f2) + round(Int, f1+f2-fd)
    jd = round(Int, d) + 1

    # Express day in Gregorian calendar
    f = jd + 1401 + (((4*jd + 274277)÷146097)*3)÷4 - 38
    e = 4*f + 3
    g = mod(e, 1461)÷4
    h = 5 * g + 2
    D = mod(h, 153)÷5 + 1
    M = mod(h÷153 + 2, 12) + 1
    Y = e÷1461 - 4716 + (12+2-M)÷12
    return Y, M, D, fd
end

"""
    jd2calhms(dj1::N1, dj2::N2) where {N1 <: Real, N2 <: Real}

Julian Date to Gregorian year, month, day, hour, minute, seconds.

### Inputs

-  `dj1,dj2` -- Julian Date (Notes 1, 2)

### Outputs 

- `iyear` -- year
- `imonth` -- month 
- `id` -- day 
- `hour`
- `minute`
- `second`
"""
function jd2calhms(dj1::N1, dj2::N2) where {N1 <: Real, N2 <: Real}
    y, m, d, fd = jd2cal(dj1, dj2)
    h, min, sec = fd2hms(fd)
    return y, m, d, h, min, sec 
end


const LEAP_TABLE = (
    ( 1972,  1, 10.0       ),
    ( 1972,  7, 11.0       ),
    ( 1973,  1, 12.0       ),
    ( 1974,  1, 13.0       ),
    ( 1975,  1, 14.0       ),
    ( 1976,  1, 15.0       ),
    ( 1977,  1, 16.0       ),
    ( 1978,  1, 17.0       ),
    ( 1979,  1, 18.0       ),
    ( 1980,  1, 19.0       ),
    ( 1981,  7, 20.0       ),
    ( 1982,  7, 21.0       ),
    ( 1983,  7, 22.0       ),
    ( 1985,  7, 23.0       ),
    ( 1988,  1, 24.0       ),
    ( 1990,  1, 25.0       ),
    ( 1991,  1, 26.0       ),
    ( 1992,  7, 27.0       ),
    ( 1993,  7, 28.0       ),
    ( 1994,  7, 29.0       ),
    ( 1996,  1, 30.0       ),
    ( 1997,  7, 31.0       ),
    ( 1999,  1, 32.0       ),
    ( 2006,  1, 33.0       ),
    ( 2009,  1, 34.0       ),
    ( 2012,  7, 35.0       ),
    ( 2015,  7, 36.0       ),
    ( 2017,  1, 37.0       )
)
const LEAP_RELEASE = 2021;


"""
    leapseconds(iyear::N, imonth::N) where {N<:Integer}

For a given UTC date, calculate Delta(AT) = TAI-UTC.

!!!! warning 
    A new version of the tables called in this function must be produced 
    whenever a new leap second is announced. 

### Inputs 

- `iyear` -- UTC year 
- `imonth` -- UTC month

### Output 

- `Δt` -- TAI - UTC in seconds

### References 

- [ERFA software library](https://github.com/liberfa/erfa/blob/master/src/dat.c)
"""
function leapseconds(iyear::N, imonth::N) where {N<:Integer}

    # If pre-UTC year, set warning status and return 0.0
    if iyear < LEAP_TABLE[1][1] 
        @warn "[Tempo] UTC not available for year $iyear, 0 is returned"
        return 0.0
    end

    # If suspiciously late year, proceed
    if iyear > LEAP_RELEASE + 5
        @debug "[Tempo] current year is 5 years ahead of the leapsecond release."
    end

    # Combine year and month to form a date-ordered integer
    m = 12*iyear + imonth
    # ...and use it to find the preceding table entry.
    @inbounds for i in length(LEAP_TABLE):-1:1
        if m >= 12*LEAP_TABLE[i][1] + LEAP_TABLE[i][2]
            return LEAP_TABLE[i][3]
        end
    end
end

#######
# UTC #
#######

"""
    utc2tai(utc1, utc2)

Time scale transformation:  Coordinated Universal Time, [`UTC`](@ref), to 
International Atomic Time, [`TAI`](@ref).

### Input 

- `utc1, utc2` --  UTC as a 2-part (quasi) Julian Date

### Output

- `tai1, tai2` -- TAI as a 2-part Julian Date


### Notes 

1. utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

2. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the JD day represents UTC days whether the
    length is 86399, 86400 or 86401 SI seconds.  

### References
    
- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 12.92 (p604).

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- [ERFA software library](https://github.com/liberfa/erfa/blob/master/src/utctai.c)
"""
function utc2tai(utc1, utc2)

    # Put the two parts of the UTC into big-first order
    big1 = abs(utc1) >= abs(utc2)

    if big1 
        u1 = utc1
        u2 = utc2
    else
        u1 = utc2
        u2 = utc1
    end
    
    # Get TAI-UTC at 0h today
    iy, im, id, fd = jd2cal(u1, u2)
    Δt0 = leapseconds(iy, im)
    
    # Get TAI-UTC at 0h tomorrow (to detect jumps)
    iyt, imt, _, _ = jd2cal(u1+1.5, u2-fd)
    Δt24 = leapseconds(iyt, imt)

    # Detect any jump
    # Remove any scaling applied to spread leap into preceding day
    fd = fd*(86400.0+Δt24-Δt0)/86400.0

    # Today's calendar date to 2-part JD
    z1, z2 = cal2jd(iy, im, id)
    
    # Assemble the TAI result, preserving the UTC split and order
    a2 = z2 - 0.5 # translate to midnight
    a2 += fd + Δt0/86400.0

    if big1 
        tai1 = u1 
        tai2 = a2 
    else 
        tai1 = a2 
        tai2 = u1 
    end
    return tai1, tai2 

end


"""
    tai2utc(tai1, tai2)

Time scale transformation:  International Atomic Time, [`TAI`](@ref) to 
Coordinated Universal Time, [`UTC`](@ref).
.

### Input 

- `tai1, tai2` -- TAI as a 2-part Julian Date

### Output

- `utc1, utc2` --  UTC as a 2-part (quasi) Julian Date


### Notes 

1. tai1+tai2 is Julian Date, apportioned in any convenient way
    between the two arguments, for example where tai1 is the Julian
    Day Number and tai2 is the fraction of a day.  The returned utc1 
    and utc2 form an analogous pair.

2. JD cannot unambiguously represent UTC during a leap second unless
    special measures are taken.  The convention in the present
    function is that the JD day represents UTC days whether the
    length is 86399, 86400 or 86401 SI seconds.  

### References
    
- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992),
    Section 12.92 (p604).

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- [ERFA software library](https://github.com/liberfa/erfa/blob/master/src/taiutc.c)
"""
function tai2utc(tai1, tai2)
    # Put the two parts of the UTC into big-first order
    big1 = abs(tai1) >= abs(tai2)
    if big1 
        a1 = tai1
        a2 = tai2
    else
        a1 = tai2
        a2 = tai1
    end

    # Initial guess for UTC
    u1 = a1
    u2 = a2
    #  Iterate (though in most cases just once is enough)
    for _ in 1:2
        g1, g2 = utc2tai(u1, u2)
        u2 += a1 - g1 
        u2 += a2 - g2 
    end

    # Return the UTC result, preserving the TAI order
    if big1 
        utc1 = u1 
        utc2 = u2 
    else 
        utc1 = u2 
        utc2 = u1 
    end
    return utc1, utc2

end

const DAY2SEC = 86400.0
const YEAR2SEC = 60.0 * 60.0 * 24.0 * 365.25
const CENTURY2SEC = 60.0 * 60.0 * 24.0 * 365.25 * 100.0
const CENTURY2DAY = 36525.0

"""
    j2000(jd)
    j2000(jd1, jd2)

Convert Julian Date (in days) to J2000 Julian Date
"""
@inline j2000(jd) = jd - DJ2000
@inline j2000(jd1, jd2) = abs(jd1) > abs(jd2) ? (jd1-DJ2000)+jd2 : (jd2-DJ2000)+jd1

"""
    j2000s(jd)
    j2000s(jd1, jd2)

Convert Julian Date (in days) in seconds past J2000 
"""
@inline j2000s(jd) = j2000(jd) * DAY2SEC
@inline j2000s(jd1, jd2) = j2000(jd1, jd2) * DAY2SEC

"""
    j2000c(jd)
    j2000c(jd1, jd2)

Convert Julian Date (in days) to Julian Centuries
"""
@inline j2000c(jd) = j2000(jd)/CENTURY2DAY
@inline j2000c(jd1, jd2) = j2000(jd1, jd2)/CENTURY2DAY
