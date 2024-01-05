using FrameTransformations
using Ephemerides
using LinearAlgebra
using ReferenceFrameRotations

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp";

eph = EphemerisProvider([download(url_spk), download(url_pck)])

FRAMES = FrameSystem{3, Float64}(eph)

@axes GCRF 1 GeocentricCelestialReferenceFrame

add_axes_inertial!(FRAMES, GCRF)

@point SSB 0
@point EMB 3
@point Sun 10
@point Earth 399
@point Moon 301

add_point_root!(FRAMES, SSB, GCRF)
add_point_ephemeris!(FRAMES, EMB)
add_point_ephemeris!(FRAMES, Earth)
add_point_ephemeris!(FRAMES, Moon)
add_point_ephemeris!(FRAMES, Sun)

tpc = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc"
iau = load(TPC(download(tpc)));

@axes IAU_EARTH 3
@axes IAU_MOON 4

add_axes_bcrtod!(FRAMES, IAU_EARTH, Earth, iau)
add_axes_bcrtod!(FRAMES, IAU_MOON, Moon, iau)

@axes LME2000 5

add_axes_bci2000!(FRAMES, LME2000, Moon, iau);

url_eop = "https://datacenter.iers.org/data/csv/finals2000A.data.csv"
eopfile = "iau2000a"

Orient.prepare_eop(download(url_eop), eopfile)
Orient.init_eop(eopfile * ".eop.dat")

@axes ITRF 6
@axes MOONPA_DE421 31006

add_axes_itrf!(FRAMES, ITRF, GCRF)
add_axes_pa421!(FRAMES, MOONPA_DE421)

FRAMES

@point SC -1_900_000

add_point_updatable!(FRAMES, SC, Moon, LME2000)

e = Epoch("2020-01-01T12:45:30.0 TDB");
x = [2274.0, 0.0, 0.0, 0.0, sqrt(4904.87/2274.0), 0.0];

update_point!(FRAMES, SC, x, e)

vector6(FRAMES, Moon, SC, LME2000, e)

vector6(FRAMES, Earth, SC, GCRF, e)

vector6(FRAMES, Earth, SC, IAU_EARTH, e)

vector6(FRAMES, Earth, SC, ITRF, e)

vector3(FRAMES, Moon, SC, IAU_MOON, e)

vector3(FRAMES, Moon, SC, MOONPA_DE421, e)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
