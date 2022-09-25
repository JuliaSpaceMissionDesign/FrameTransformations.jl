using Calceph: Ephem as CalcephEphemHandler, 
               unsafe_compute!, useNaifId, unitKM, unitSec,
               getpositionrecords, gettimespan, gettimescale

using NodeGraphs: NodeGraph
import Basic: register!, connect!
using Basic.Bodies: NAIFId

export CalcephProvider, 
       ephem_timespan, ephem_timescale, ephem_position_records

"""
    CalcephProvider

Calceph-based ephemeris reading.
"""
struct CalcephProvider <: AbstractEphemerisProvider
    ptr::CalcephEphemHandler
    function CalcephProvider(files::Vector{<:AbstractString})
        new(CalcephEphemHandler(files))
    end
end
CalcephProvider(file::AbstractString) = CalcephProvider([file])

"""
    ephem_position_records(eph::CalcephProvider)

Get ephemeris an array of `Calceph.PositionRecord`s, providing detailed 
informations on the content of the ephemeris file.
"""
ephem_position_records(eph::CalcephProvider) = getpositionrecords(eph.ptr)

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
ephem_timespan(eph::CalcephProvider) = gettimespan(eph.ptr)

"""
    ephem_timescale(eph::CalcephProvider) 

Retrieve `Basic` timescale associated with ephemeris handler `eph`.
"""
function ephem_timescale(eph::CalcephProvider) 
    tsid = gettimescale(eph.ptr)
    if tsid == 1
        return TDB
    elseif tsid == 2
        return TCB
    else
        throw(error("[Ephemeris] unknown time scale identifier: $tsid"))
    end
end

"""
    register!(g::NodeGraph, eph::CalcephProvider) 

Register and connect the bodies present in the ephemeris file in a body graph.
"""
function register!(g::NodeGraph{NAIFId, N, G, N}, 
    eph::CalcephProvider) where {N<:Integer, G}
    pos_records = ephem_position_records(eph)
    for p in pos_records
        connect!(g, NAIFId(p.center), NAIFId(p.target))
    end
end

