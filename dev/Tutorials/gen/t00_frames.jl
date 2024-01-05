using FrameTransformations #hide

F = FrameSystem{2, Float64}()

F = FrameSystem{2, Float64, InternationalAtomicTime}()

using Ephemerides, Downloads

url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(url));

F = FrameSystem{2, Float64}(eph)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
