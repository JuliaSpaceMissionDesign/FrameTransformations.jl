using BenchmarkTools 

include("frames/Frames.jl")

path = "/home/michele/spice/kernels/"; 
f1 = path*"spk/de440.bsp";

eph = EphemerisKernels([f1])
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
@ephemeris_point FRAMES moon 301 earth type=MoonPoint
@ephemeris_point FRAMES venus 299 ssb type=VenusPoint
@fixed_point FRAMES sc 302 earth ECLIPJ2000 offset type=EarthSpacecraftPoint

ys = MVector{6, Float64}(undef);

unsafe_compute_order!(ys, eph, ep, 0., 3, 0, cph_use_naifid+cph_km+cph_sec, 0);
ys[1:3] - get_vector3(FRAMES, ssb, emb, ICRF, ep)


unsafe_compute_order!(ys, eph, ep, 0., 399, 0, cph_use_naifid+cph_km+cph_sec, 0);
ys[1:3] - get_vector3(FRAMES, ssb, earth, ICRF, ep)


nruns = 100
y = Vector{SVector{3, Float64}}(undef, nruns);
tid = Vector{Int}(undef, nruns)
Threads.@threads for i = 1:nruns 
    y[i] = get_vector3(FRAMES, emb, earth, ECLIPJ2000, ep)
    tid[i] = Threads.threadid()
    get_vector3(FRAMES, sc, ssb, ICRF, ep+2)
end

@benchmark get_vector3($FRAMES, $ssb, $earth, $ICRF, $ep+2)

# ys[1:3]
# get_vector3(FRAMES, ssb, moon, ICRF, ep)

# @benchmark get_vector3($FRAMES, $ssb, $moon, $ICRF, $ep)
# @benchmark get_vector3($FRAMES, $venus, $sc, $ECLIPJ2000, $ep)

# @benchmark unsafe_compute_order!($ys, $eph, $ep, $0., $301, $0, $cph_use_naifid+$cph_km+$cph_sec, $0)
# @benchmark unsafe_compute_order!($ys, $eph, $ep, $0., $301, $0, $cph_use_naifid+$cph_km+$cph_sec, $1)

tid = Vector{Int}(undef, nruns);
y = Vector{MVector{3, Float64}}(undef, nruns);
Threads.@threads for i = 1:nruns 
    y[i] = MVector{3, Float64}(zeros(3)...)
    unsafe_compute_order!(y[i], eph, ep, 0., 399, 0, cph_use_naifid+cph_km+cph_sec, 0);
    tid[i] = Threads.threadid()
end

yex = zeros(3);
unsafe_compute_order!(yex, eph, ep, 0., 399, 3, cph_use_naifid+cph_km+cph_sec, 0);
err = zeros(nruns)
for i = 1:nruns 
    err[i] = norm(yex - y[i])
end

maximum(err)


nruns = 100
y = Vector{SVector{3, Float64}}(undef, nruns);
tid = Vector{Int}(undef, nruns);
eph = Ephem(f1);
CALCEPH.prefetch(eph)
Threads.@threads for i = 1:nruns 
    CALCEPH.compute(eph, ep, 0., 399, 0, cph_km+cph_sec+cph_use_naifid)
    tid[i] = Threads.threadid()
end
