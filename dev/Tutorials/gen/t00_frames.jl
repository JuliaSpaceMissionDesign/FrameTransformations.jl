using FrameTransformations
using Tempo

F = FrameSystem{2, Float64}()

F = FrameSystem{2, Float64, InternationalAtomicTime}()

has_point(F, 1)

has_axes(F, 1)

has_direction(F, :Root)

axes(F)

points(F)

order(F)

FrameTransformations.timescale(F)

using Ephemerides, Downloads

url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
E = EphemerisProvider(Downloads.download(url));

F = FrameSystem{2, Float64}()

add_axes_icrf!(F)
add_point_root!(F, :SSB, 0, 1)

add_point_ephemeris!(F, E, :Sun, 10)
add_point_ephemeris!(F, E, :EMB, 3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
