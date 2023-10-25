module Orient

using DelimitedFiles
using LinearAlgebra
using Logging

using ReferenceFrameRotations
using RemoteFiles: @RemoteFile, download, path
using StaticArrays

using JSMDInterfaces.Ephemeris:
    AbstractEphemerisProvider, ephem_available_axes, ephem_orient!

using JSMDInterfaces
using JSMDInterfaces.FilesIO
using JSMDInterfaces.Math: interpolate

using JSMDUtils.Math: InterpAkima, arcsec2rad, skew

using Tempo

# Topocentric 
include("geodesy.jl")

# Earth
include("Earth/Earth.jl")

# Moon 
include("moon.jl")

# Planets
include("tpc.jl")
include("planets.jl")

# Ecliptic 
include("ecliptic.jl")
include("legacy.jl")

function __init__()
    if !Tempo.has_timescale(TIMESCALES, Tempo.timescale_id(UT1))
        Tempo.add_timescale!(TIMESCALES, UT1, offset_utc2ut1; parent=UTC)
    end
end

end
