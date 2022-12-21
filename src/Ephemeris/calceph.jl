export CalcephProvider, 
       ephem_timespan, 
       ephem_timescale, 
       ephem_position_records,
       ephem_orient_records,
       ephem_compute_order!

using CALCEPH: Ephem as CalcephEphemHandler, 
               prefetch, 
               timespan, 
               timeScale,
               positionRecords,
               orientationRecords,
               unsafe_compute!, 
               useNaifId, unitKM, unitSec

using Basic: AstronautGenericError, EphemerisError

"""
    CalcephProvider

Calceph-based ephemeris handler.
"""
struct CalcephProvider <: AbstractEphemerisProvider
    ptr::CalcephEphemHandler
    function CalcephProvider(files::Vector{<:AbstractString})
        ptr = CalcephEphemHandler(unique(files))
        prefetch(ptr)
        new(ptr)
    end
end
CalcephProvider(file::AbstractString) = CalcephProvider([file])

function ephem_load(::Type{CalcephProvider}, files::Vector{<:AbstractString})
    return CalcephProvider(files)
end

"""
    ephem_position_records(eph::CalcephProvider)

Get ephemeris an array of `CALCEPH.PositionRecord`s, providing detailed 
informations on the content of the ephemeris file.
"""
ephem_position_records(eph::CalcephProvider) = positionRecords(eph.ptr)


"""
    ephem_orient_records(eph::CalcephProvider)

Get ephemeris an array of `CALCEPH.OrientationRecord`s, providing detailed 
informations on the content of the ephemeris file.
"""
ephem_orient_records(eph::CalcephProvider) = orientationRecords(eph.ptr)

"""
    ephem_timespan(eph::CalcephProvider)

Returns the first and last time available in the ephemeris file associated to 
an ephemeris file.

### Input/s:

- `eph` : ephemeris

### Output/s:

Returns a tuple containing:

- `firsttime` -- Julian date of the first time
- `lasttime` -- Julian date of the last time
- `continuous` -- information about the availability of the quantities over the 
                  time span
    
    It returns the following value in the parameter continuous :

        1. if the quantities of all bodies are available for any time between 
            the first and last time.
        2. if the quantities of some bodies are available on discontinuous time 
            intervals between the first and last time.
        3. if the quantities of each body are available on a continuous time 
            interval between the first and last time, but not available for any 
            time between the first and last time.

"""
ephem_timespan(eph::CalcephProvider) = timespan(eph.ptr)

"""
    ephem_timescale(eph::CalcephProvider) 

Retrieve `Basic` timescale associated with ephemeris handler `eph`.
"""
function ephem_timescale(eph::CalcephProvider) 
    tsid = timeScale(eph.ptr)
    if tsid == 1
        return TDB
    elseif tsid == 2
        return TCB
    else
        throw(
            AstronautGenericError(
                String(Symbol(@__MODULE__)),
                "unknown time scale identifier: $tsid"
            )
        )
    end
end

"""
    ephem_compute_order!(res, eph, jd0, time, target, center, order)

### Inputs
- `res` -- results
- `eph` -- ephemeris provider
- `jd0` and `time` -- Two-parts Julian date
"""
function ephem_compute_order!(res, eph::AbstractEphemerisProvider, jd0, time, target, 
    center, order) end

# TODO: jd0+time da scrivere come datetime!
function ephem_compute_order!(res, eph::CalcephProvider, jd0::Float64, time::Float64, 
    target::Int, center::Int, order::Int)
    stat = unsafe_compute!(res, eph.ptr, jd0, time, target, center, useNaifId+unitKM+unitSec, order)
    stat == 1 && throw(EphemerisError(String(Symbol(@__MODULE__)), "ephemeris data for "*
                    "point with NAIFId $target with respect to point $center is not available "*
                    "at JD $(jd0+time)"))
    nothing
end

# function ephem_compute_order!(res, eph::SpiceProvider, jd0::Float64, time::Float64, 
#     target::Int, center::Int, order::Int)

#     et = ((jd0 + time) - DJ2000) * 86400

#     unsafe_compute!(res, eph.ptr, jd0, time, target, center, useNaifId+unitKM+unitSec, order)
#     nothing
# end