
"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: AbstractEphemerisProvider end 

# TODO: l'ultimo valore cosa era?
ephem_timespan(::NullEphemerisProvider) = (0, 0, 1)
ephem_timescale(::NullEphemerisProvider) = TDB