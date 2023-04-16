
"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: AbstractEphemerisProvider end

ephem_timescale(::NullEphemerisProvider) = TDB

# Precompilation 

precompile(ephem_timescale, (NullEphemerisProvider,))