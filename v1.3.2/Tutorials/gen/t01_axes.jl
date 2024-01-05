using FrameTransformations
using ReferenceFrameRotations

dcm  = angle_to_dcm(π/3, :Z)
δdcm = DCM(0I)

R = Rotation(dcm, δdcm)

R[1]

R[2]

v = [1., -6., 3., 0., 5., 0]
R*v

inv(R)

using Ephemerides

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

eph = EphemerisProvider([download(url_spk), download(url_pck)])

G = FrameSystem{3, Float64}(eph)

@axes ICRF 1 InternationalCelestialReferenceFrame

add_axes_inertial!(G, ICRF)

@axes ECLIPJ2000 17

add_axes_inertial!(G, ECLIPJ2000; parent=ICRF, dcm=DCM_ICRF_TO_ECLIPJ2000)

R = rotation6(G, ICRF, ECLIPJ2000, 10.0)

R[1]

R[2]

@axes FO1 2

rot = angle_to_dcm(π/4, :Z)

add_axes_fixedoffset!(G, FO1, ICRF, rot)

R = rotation6(G, ICRF, FO1, 86400)

R[1]

R[2]

@axes RotAx 3

fun(t) = angle_to_dcm(-t, :Z)

add_axes_rotating!(G, RotAx, FO1, fun)

R1 = rotation6(G, ICRF, RotAx, π/4)

R1[1]

R2 = rotation6(G, ICRF, RotAx, π/2)

R2[2]

using JSMDUtils.Math

@axes RotAx2 4

fun(t) = angle_to_dcm(-t, :Z)
dfun(t) = (angle_to_dcm(-t, :Z), Math.angle_to_δdcm([-t, -1], :Z))

add_axes_rotating!(G, RotAx2, FO1, fun, dfun)

R2 = rotation6(G, ICRF, RotAx2, π/2)

R2[2]

@axes ProjAx 500
@axes RotAx3 501

fun(t) = angle_to_dcm(-t, :Z)

add_axes_rotating!(G, RotAx3, ICRF, fun)
add_axes_projected!(G, ProjAx, ICRF, fun)

R1 = rotation6(G, ICRF, RotAx3, 50.0)
R2 = rotation6(G, ICRF, ProjAx, 50.0)

R1[1] - R2[1]

R1[2]

R2[2]

@axes SunFrame 600

@point SSB 0 SolarSystemBarycenter
@point Sun 10 SunPoint

add_point_root!(G, SSB, ICRF)
add_point_ephemeris!(G, Sun)

v1 = ComputableAxesVector(Sun, SSB, 1)
v2 = ComputableAxesVector(Sun, SSB, 2)

add_axes_computable!(G, SunFrame, ICRF, v1, v2, :XY)

R = rotation6(G, ICRF, SunFrame, 0.0)

@axes MoonPA 31006

add_axes_ephemeris!(G, MoonPA, :ZXZ)

R = rotation9(G, ICRF, MoonPA, 86400.0)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
