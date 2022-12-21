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

