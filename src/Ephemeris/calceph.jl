export CalcephProvider, 
       ephem_timespan, 
       ephem_timescale, 
       ephem_position_records,
       ephem_orient_records,
       ephem_available_points, 
       ephem_available_axes,
       ephem_compute_order!,
       ephem_orient_order!,
       ephem_load

using CALCEPH: Ephem as CalcephEphemHandler, 
               prefetch, 
               timespan, 
               timeScale,
               positionRecords,
               orientationRecords,
               unsafe_compute!, unsafe_orient!,
               useNaifId, unitKM, unitSec, unitRad

using Basic: AstronautGenericError, EphemerisError

"""
    CalcephProvider(file::String)
    CalcephProvider(files::Vector{String})
    
Create a `CalcephProvider` instance by loading a single or multiples ephemeris kernel 
files specified by `files`.

!!! note 
    Since this object is immutable, kernels cannot be added nor removed from the generated
    `CalcephProvider` instance.  

### Example 
```jldoctest
julia> eph1 = CalcephProvider("PATH_TO_KERNEL")
CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} [...]))

julia> eph2 = CalcephProvider(["PATH_TO_KERNEL_1", "PATH_TO_KERNEL_2"])
CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} [...]))
```
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

Get an array of `CALCEPH.PositionRecord`, providing detailed informations on the content of 
the ephemeris file.
"""
ephem_position_records(eph::CalcephProvider) = positionRecords(eph.ptr)


"""
    ephem_available_points(eph::CalcephProvider)

Return a list of NAIFIds representing bodies with available ephemeris data. 
"""
function ephem_available_points(eph::CalcephProvider)
    try 
        rec = ephem_position_records(eph)
        tids = map(x->x.target, rec)
        cids = map(x->x.center, rec)

        return unique([tids..., cids...])
    catch 
        return Int64[]
    end
end


"""
    ephem_orient_records(eph::CalcephProvider)

Get ephemeris an array of `CALCEPH.OrientationRecord`s, providing detailed 
informations on the content of the ephemeris file.
"""
ephem_orient_records(eph::CalcephProvider) = orientationRecords(eph.ptr)


"""
    ephem_available_points(eph::CalcephProvider)

Return a list of Frame IDs representing axes with available orientation data. 
"""
function ephem_available_axes(eph::CalcephProvider)
    try 
        rec = ephem_orient_records(eph)
        tids = map(x->x.target, rec)
        cids = map(x->x.frame, rec)

        return unique([tids..., cids...])
    catch 
        return Int64[]
    end
end

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

Interpolate the position and/or its derivatives up to `order` for one body `target` relative 
to another `center` at the time `jd0` + `time`, expressed as a Julian Date. This function reads 
the ephemeris files associated to `eph` and stores the results to `res`.

The returned array `res` must be large enough to store the results. The size of this array 
must be equal to 3*order: 

- res[1:3] contain the position (x, y, z) and is always valid 
- res[4:6] contain the velocity (dx/dt, dy/dt, dz/dt) for order ≥ 1 
- res[7:9] contain the acceleration (d²x/dt², d²y/dt², d²z/dt²) for order ≥ 2
- res[10:12] contain the jerk (d³x/dt³, d³y/dt³, d³z/dt³) for order ≥ 3

The values stores in `res` are always returned in km, km/s, km/s², km/s³

### See also 
See also [`ephem_orient_order!`](@ref)
"""
function ephem_compute_order!(res, eph::CalcephProvider, jd0::Float64, time::Float64, 
            target::Int, center::Int, order::Int)
    stat = unsafe_compute!(res, eph.ptr, jd0, time, target, center, useNaifId+unitKM+unitSec, order)
    stat == 0 && throw(EphemerisError(String(Symbol(@__MODULE__)), "ephemeris data for "*
                    "point with NAIFId $target with respect to point $center is not available "*
                    "at JD $(jd0+time)"))
    nothing
end


"""
    ephem_orient_order!(res, eph, jd0, time, target, order)

Interpolate the orientation and its derivatives up to `order` for the `target` body at the 
time `jd0` + `time`, expressed as a Julian Date. This function reads the ephemeris files 
associated to `eph` and stores the results to `res`.

The returned array `res` must be large enough to store the results. The size of this array 
must be equal to 3*order: 

- res[1:3] contain the angles 
- res[4:6] contain the 1st derivative  for order ≥ 1 
- res[7:9] contain the 2nd derivative for order ≥ 2
- res[10:12] contain the 3rd derivative for order ≥ 3

The values stores in `res` are always returned in rad, rad/s, rad/s², rad/s³

### See also 
See also [`ephem_compute_order!`](@ref)
"""
function ephem_orient_order!(res, eph::CalcephProvider, jd0::Float64, time::Float64, 
            target::Int, order::Int)

    stat = unsafe_orient!(res, eph.ptr, jd0, time, target, useNaifId+unitRad+unitSec, order)
    stat == 0 && throw(EphemerisError(String(Symbol(@__MODULE__)), "ephemeris data for "*
                "frame with NAIFId $target is not available at JD $(jd0+time)"))
    nothing
end
