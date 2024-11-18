using FrameTransformations
using Ephemerides

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

const EPH = EphemerisProvider([download(url_spk), download(url_pck)])
const F = FrameSystem{3,Float64}()

add_axes!(F, :ICRF, AXESID_ICRF)

using ReferenceFrameRotations
using LinearAlgebra

fun(t) = DCM(1.0I)
add_axes_projected!(F, :GCRF, AXESID_GCRF, :ICRF, fun)

R = rotation6(F, AXESID_ICRF, AXESID_GCRF, 1.0)

R[1]

R[2]

rot = angle_to_dcm(π / 4, :Z)

add_axes_fixedoffset!(F, :FOX, 2, AXESID_ICRF, rot)

R = rotation6(F, :ICRF, :FOX, 86400)

R[1]

R[2]

fun(t) = angle_to_dcm(-t, :Z)

add_axes_rotating!(F, :ROX, 3, :ICRF, fun)

R = rotation6(F, 2, 3, π / 4)

R[1]

R2 = rotation6(F, 1, 3, π / 4)

R2[2]

using JSMDUtils.Math

fun(t) = angle_to_dcm(-t, :Z)
dfun(t) = (angle_to_dcm(-t, :Z), Math.angle_to_δdcm([-t, -1], :Z))

add_axes_rotating!(F, :ROX2, 4, :ICRF, fun, dfun)

R2 = rotation6(F, 1, 3, π / 4)

R2[2]

FrameTransformations.add_axes_ephemeris!(F, EPH, :MOONPA, 31006, :ZXZ)

R = rotation6(F, :ICRF, :MOONPA, 86400.0)

R[1]

R[2]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
