using BenchmarkTools 
using ForwardDiff
using Basic 

# CR3BP scheme 
include("frames4/Frames.jl")

@axes ICRF 1 InertialReferenceFrame
@axes MEME2000 2 MeanEarthMeanEquinoxJ2000
@axes ECLIPJ2000 4 EclipticJ2000
@axes IAU_EARTH 5 IauEarth
@axes IAU_MARS 13 IauMars

eph = CalcephProvider("/home/michele/spice/kernels/spk/de440.bsp")
FRAMES = FrameSystem{Float64}(eph);

add_axes_inertial!(FRAMES, ICRF)
add_axes_inertial!(FRAMES, MEME2000; parent=1, dcm=angle_to_dcm(π/3, :Z))
add_axes_inertial!(FRAMES, ECLIPJ2000; parent=MEME2000, dcm=angle_to_dcm(5π/7, :X))
add_axes_fixedoffset!(FRAMES, IAU_EARTH, ECLIPJ2000, angle_to_dcm(3π/7, :X))
add_axes_fixedoffset!(FRAMES, IAU_MARS, ICRF, angle_to_dcm(2π/3, :X))

get_rotation3(FRAMES, ICRF, MEME2000, 0.)
get_rotation6(FRAMES, ICRF, MEME2000, 0.)
get_rotation9(FRAMES, ICRF, MEME2000, 0.)

@benchmark get_rotation3($FRAMES, $ICRF, $MEME2000, 0.)
@benchmark get_rotation6($FRAMES, $ICRF, $MEME2000, 0.)
@benchmark get_rotation9($FRAMES, $ICRF, $MEME2000, 0.)

@benchmark get_rotation3($FRAMES, $ICRF, $ECLIPJ2000, 0.)
@benchmark get_rotation6($FRAMES, $ICRF, $ECLIPJ2000, 0.)
@benchmark get_rotation9($FRAMES, $ICRF, $ECLIPJ2000, 0.)

R3 = get_rotation3(FRAMES, MEME2000, ICRF, 0.)
R6 = get_rotation6(FRAMES, ECLIPJ2000, ICRF, 0.)
R3 = get_rotation3(FRAMES, IAU_EARTH, IAU_MARS, 0.)

@benchmark get_rotation3($FRAMES, $IAU_EARTH, $IAU_MARS, 0.)

@point SSB 0 SolarSystemBaricenter
@point Sun 10 SunPoint 
@point EMB 3 EarthMoonBarycenterPoint
@point Earth 399 EarthPoint 
@point Moon 301 MoonPoint


jd0 = ephem_timespan(eph)[1] + 1000

add_point_root!(FRAMES, SSB, 1)
add_point_ephemeris!(FRAMES, Sun, SSB)
add_point_ephemeris!(FRAMES, EMB, Sun)
add_point_ephemeris!(FRAMES, Earth, EMB)
add_point_ephemeris!(FRAMES, Moon, Earth)

get_vector3(FRAMES, SSB, Sun, IAU_MARS, jd0)
get_vector6(FRAMES, SSB, Sun, IAU_MARS, jd0)
get_vector9(FRAMES, SSB, Sun, IAU_MARS, jd0)

@benchmark get_vector3($FRAMES, $SSB, $Earth, $ICRF, $jd0)
@benchmark get_vector6($FRAMES, $SSB, $Earth, $ICRF, $jd0)
@benchmark get_vector9($FRAMES, $SSB, $Earth, $ICRF, $jd0)


path = get_path(frames_points(FRAMES), point_alias(SSB), point_alias(Earth))
@code_warntype get_vector3(FRAMES, SSB, Earth, ICRF, jd0)
@code_warntype _compute_vector3(FRAMES, jd0, axes_alias(ICRF), path)
@code_warntype _get_vector3_forward(FRAMES, jd0, path)
@code_warntype _get_vector3_backwards(FRAMES, jd0, path)
@code_warntype _compute_vector3(FRAMES.points.nodes[1], FRAMES.points.nodes[2], jd0)
@code_warntype _compute_vector3(FRAMES.points.nodes[2], jd0)


node = FRAMES.points.nodes[2]
@benchmark _compute_vector3($node, $jd0)
@benchmark ephem_compute_order()
y = @MVector zeros(3)
@benchmark ephem_compute_order!($y, $eph, $jd0, 0., 10, 0, 0)