using BenchmarkTools 

include("frames/Frames.jl")

path = "/home/michele/spice/kernels/"; 
f1 = path*"spk/de440.bsp";

eph = EphemerisKernels([f1, f2])
ep = get_timespan(eph)[1] + 1000;

FRAMES = FrameSystem{Float64}(eph);

icrf2meme = angle_to_dcm(pi*2, :X)
meme2eclip = angle_to_dcm(pi*2, pi/3, :XY)

offset = SA[1., 2., 3.]

@inertial_axes FRAMES ICRF 1 type=InternationalCelestialReferenceFrame
@inertial_axes(FRAMES, MEME2000, 2, type=MeanEquatorMeanEquinoxJ2000, 
                parent=ICRF, dcm=icrf2meme)

# Si possono fare sia come inerziali che come fissi
@fixed_axes FRAMES ECLIPJ2000 3 MEME2000 meme2eclip type=EclipticEquinoxJ2000

@root_point FRAMES ssb 0 1 type=SolarSystemBarycenterPoint

@ephemeris_point FRAMES sun 10 ssb type=SunPoint
@ephemeris_point FRAMES emb 3 ssb type=EarthMoonBarycenterPoint
@ephemeris_point FRAMES earth 399 emb type=EarthPoint
@ephemeris_point FRAMES moon 301 emb type=MoonPoint
@ephemeris_point FRAMES venus 299 ssb type=VenusPoint
@fixed_point FRAMES sc 302 earth ECLIPJ2000 offset type=EarthSpacecraftPoint


ys = MVector{6, Float64}(undef);

unsafe_compute_order!(ys, eph, ep, 0., 3, 0, cph_use_naifid+cph_km+cph_sec, 0);
ys[1:3] - get_vector3(FRAMES, ssb, emb, ICRF, ep)


unsafe_compute_order!(ys, eph, ep, 0., 399, 0, cph_use_naifid+cph_km+cph_sec, 0);
ys[1:3] - get_vector3(FRAMES, ssb, earth, ICRF, ep)


# ys[1:3]
# get_vector3(FRAMES, ssb, moon, ICRF, ep)

# @benchmark get_vector3($FRAMES, $ssb, $moon, $ICRF, $ep)
# @benchmark get_vector3($FRAMES, $venus, $sc, $ECLIPJ2000, $ep)

# @benchmark unsafe_compute_order!($ys, $eph, $ep, $0., $301, $0, $cph_use_naifid+$cph_km+$cph_sec, $0)
# @benchmark unsafe_compute_order!($ys, $eph, $ep, $0., $301, $0, $cph_use_naifid+$cph_km+$cph_sec, $1)
