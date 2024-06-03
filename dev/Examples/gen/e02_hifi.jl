using FrameTransformations
using Ephemerides
using LinearAlgebra
using ReferenceFrameRotations

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp";

EPH = EphemerisProvider([download(url_spk), download(url_pck)])

FRAMES = FrameSystem{3, Float64}()

add_axes_icrf!(FRAMES)

add_point_root!(FRAMES, :SSB, 0, 1)
add_point_ephemeris!(FRAMES, EPH, :EMB, 3)
add_point_ephemeris!(FRAMES, EPH, :Sun, 10)
add_point_ephemeris!(FRAMES, EPH, :Earth, 399)
add_point_ephemeris!(FRAMES, EPH, :Moon, 301)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
