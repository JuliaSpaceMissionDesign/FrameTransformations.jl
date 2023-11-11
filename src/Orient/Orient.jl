module Orient

using DelimitedFiles
using LinearAlgebra
using Logging

using ReferenceFrameRotations
using StaticArrays

using JSMDInterfaces.Ephemeris:
    AbstractEphemerisProvider, ephem_available_axes, ephem_orient!

using JSMDInterfaces.Math: interpolate
using JSMDUtils.Math: InterpAkima, arcsec2rad

using Tempo
using FrameTransformations.Utils: skew

# Earth
include("Earth/Earth.jl")

# Moon 
include("moon.jl")

# Planets
include("planets.jl")

# Ecliptic 
include("ecliptic.jl")

function __init__()
    if !Tempo.has_timescale(TIMESCALES, Tempo.timescale_id(UT1))
        Tempo.add_timescale!(TIMESCALES, UT1, offset_utc2ut1; parent=UTC)
    end
end

end
