"""
    AbstractEphemerisProvider

Abstract type to represent ephemeris types.
"""
abstract type AbstractEphemerisProvider end

"""
    ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}

Load a generic ephemeris file.
"""
function ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}
    return E(file)
end

# TODO: add NotImplementedErrors to the abstract implementation of those functions.

function ephem_position_records(::AbstractEphemerisProvider) end

function ephem_available_points(::AbstractEphemerisProvider) end 

function ephem_orint_records(::AbstractEphemerisProvider) end

function ephem_available_axes(::AbstractEphemerisProvider) end

function ephem_timespan(::AbstractEphemerisProvider) end

function ephem_timescale(::AbstractEphemerisProvider) end

function ephem_compute_order!(res, eph::AbstractEphemerisProvider, 
    jd::Number, time::Number, target::Int, center::Int, order::Int) end

function ephem_orient_order!(res, eph::AbstractEphemerisProvider, 
    jd::Number, time::Number, target::Int, order::Int) end
