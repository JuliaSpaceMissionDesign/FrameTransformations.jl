using FrameTransformations

G = FrameSystem{2, Float64}()

@axes SATF 1 SatelliteFrame

add_axes_inertial!(G, SATF)

@point SC -10000 SpacecraftCenter

add_point_root!(G, SC, SATF)

G

@point SACL -10101 SolarArrayCenterLeft
@point SACR -10102 SolarArrayCenterRight
@point Antenna -10001

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixed!(G, SACL, SC, SATF, sa_offset_left)
add_point_fixed!(G, SACR, SC, SATF, sa_offset_right)
add_point_fixed!(G, Antenna, SC, SATF, an_offset)

vector3(G, SACL, SC, SATF, 0.0)

vector6(G, Antenna, SACR, SATF, 10.0)

@point TimeDependantAppendage -10003

fun(t) = [cos(t), sin(t), 0]

add_point_dynamical!(G, TimeDependantAppendage, SACL, SATF, fun)

vector6(G, TimeDependantAppendage, SC, SATF, π/3)

@point TimeDependantAppendage2 -10004

fun(t) = [cos(t), sin(t), 0]
dfun(t) = [cos(t), sin(t), 0, -sin(t), cos(t), 0]

add_point_dynamical!(G, TimeDependantAppendage2, SACL, SATF, fun, dfun)

vector6(G, TimeDependantAppendage2, SC, SATF, π/3)

@point UA -10002 UpdatableAppendage

add_point_updatable!(G, UA, SC, SATF)

ua_pos = [0.0, -1.0, 0.0]
update_point!(G, UA, ua_pos, 0.0)
vector3(G, Antenna, UA, SATF, 0.0)

using Ephemerides

spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(download(spk))

F = FrameSystem{2, Float64}(eph)

@axes ICRF 1

@point SSB 0
@point Sun 10
@point EMB 3
@point Earth 399

add_axes_inertial!(F, ICRF)
add_point_root!(F, SSB, ICRF)
add_point_ephemeris!(F, Sun)
add_point_ephemeris!(F, EMB)
add_point_ephemeris!(F, Earth)

F

vector6(F, EMB, SSB, ICRF, 1000.0)

vector3(F, Earth, SSB, ICRF, 0.0)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
