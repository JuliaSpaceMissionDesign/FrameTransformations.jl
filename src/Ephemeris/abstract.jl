"""
    AbstractEphemerisProvider

Abstract type to represent ephemeris types.
"""
abstract type AbstractEphemerisProvider end

"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: AbstractEphemerisProvider end 

"""
    ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}

Load a generic ephemeris file.
"""
function ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}
    return E(file)
end

