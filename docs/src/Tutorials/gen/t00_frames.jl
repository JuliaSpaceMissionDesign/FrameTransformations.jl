using FrameTransformations #hide
using Tempo #hide

F = FrameSystem{2, Float64}()

F = FrameSystem{2, Float64, InternationalAtomicTime}()

using Ephemerides, Downloads

url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(url));

F = FrameSystem{2, Float64}()

add_axes_icrf!(F)
add_point_root!(F, :SSB, 0, 1)

add_point_ephemeris!(F, eph, :Sun, 10)
add_point_ephemeris!(F, eph, :EMB, 3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
