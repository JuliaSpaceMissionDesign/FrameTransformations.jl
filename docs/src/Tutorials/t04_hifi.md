# [Use Case: High Fidelity](@id tutorial_04_hifi)

Once the general structure of the [`FrameSystem`](@ref) is understood, we can pass to a use case in which we want to build and exploit our frame system to perform computations in a high-fidelity environment. 

### Frame system setup 

In this example, we plan on using ephemeris data to retrieve accurate positions of the planets and the orientation of certain reference frames. Therefore, we create an ephemeris provider object leveraging our own [`Ephemerides.jl`](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package and use it to generate a frame system instance:

```@setup init 
using FrameTransformations
using Downloads
```

```@repl init
using Ephemerides, Downloads

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

eph = EphemerisProvider([Downloads.download(url_spk), Downloads.download(url_pck)])

FRAMES = FrameSystem{3, Float64}(eph)
```

Once the graph is created, we assign the `GCRF` (Geocentric Celestial Reference Frame) as our set of inertial root-axes:

```@setup postInitHiFi
using FrameTransformations

FRAMES = FrameSystem{3, Float64}()
```

```@repl postInitHiFi
@axes GCRF 1 GeocentricCelestialReferenceFrame
add_axes_inertial!(FRAMES, GCRF)
```

These axes are practically coincident with the `ICRF`. To do so, we can leverage the possibility to define multiple aliases that are associated to the same integer ID. Thus let's create an alias of the `GCRF` and call it `ICRF`: 

```@repl  
@axes ICRF 1 InternationalCelestialReferenceFrame
```

In this scenario, we will be working within the Cislunar environment, therefore we will need the major bodies that influence this dynamic regime, i..e, the Earth, the Moon and the Sun. To do so, we also define the Solar System Barycenter (SSB) and the Earth-Moon Barycenter (EMB) as the ephemeris data of the remaining bodies is expressed with respect to those. 

For this example, we will assume the SSB is our root point:

```@setup postICRFHiFi
using FrameTransformations
using Ephemerides, Downloads

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

eph = EphemerisProvider([Downloads.download(url_spk), Downloads.download(url_pck)])

FRAMES = FrameSystem{3, Float64}(eph)

@axes GCRF 1 GeocentricCelestialReferenceFrame
@axes ICRF 1 InternationalCelestialReferenceFrame
add_axes_inertial!(FRAMES, GCRF)
```

```@repl postICRFHiFi
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
```

### Registering body-fixed axes

At this point, we want to be able to transform back and forth from body-fixed axes of the Earth 
and the Moon to the `GCRF`. The Earth and the Moon cases are actually _special cases_ because 
they have both high-accuracy and low-precision body-fixed rotation models. 

Let's start with the low precision ones: `IAU_EARTH` and `IAU_MOON`. To register such frames, 
we first need to parse a [`TPC`](@ref) file with the required constants:

```@repl init
tpc = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc";
iau = load(TPC(Downloads.download(tpc)));
```

At this point, we leverage the high-level interface provided by this package to automatically 
generates ad-hoc functions from the IAU Rotational Elements via the [`add_axes_bcrtod!`](@ref) function.
In particular, these axes are Body-Centric Rotating, True-of-Date axes, or BCR-TOD (which 
are the `FrameTransformations` way of representing SPICE `IAU_XXX` axes).

```@setup iauAxes
using FrameTransformations, Downloads
tpc = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc"
iau = load(TPC(Downloads.download(tpc)));

@axes GCRF 1 
@axes ICRF 1
@point Earth 399 
@point Moon 301

FRAMES = FrameSystem{2, Float64}()

add_axes_inertial!(FRAMES, ICRF)

add_point_root!(FRAMES, Earth, ICRF)
add_point_fixed!(FRAMES, Moon, Earth, ICRF, [0.0, 0.0, 0.0])
```

```@repl iauAxes
@axes IAU_EARTH 3 
@axes IAU_MOON 4

add_axes_bcrtod!(FRAMES, IAU_EARTH, Earth, iau);
add_axes_bcrtod!(FRAMES, IAU_MOON, Moon, iau);
```

The registration of this type of axes requires either an alias or the ID of the point associated to the body.

For convenience, let us also insert a set of inertial axes for the Moon. This can be done using the IAU-based Body-Centric Inertial J2000 axes available through the [`add_axes_bci2000!`](@ref) function, which defines a rotation from the `ICRF` to the local equatorial plane of the body at the J2000 epoch.

```@repl iauAxes
@axes LME2000 5

add_axes_bci2000!(FRAMES, LME2000, Moon, iau);
```

Finally, we complete the axes definition by inserting the high-precision Earth and Lunar body-fixed rotation models. 
For this purpose, `FrameTransformations` provides two high-level functions that can be used to ease these definitions: 
[`add_axes_itrf!`](@ref) and [`add_axes_pa421!`](@ref).

```@repl postICRFHiFi
@axes ITRF 6 
@axes MOONPA_DE421 31006 

add_axes_itrf!(FRAMES, ITRF, GCRF)
add_axes_pa421!(FRAMES, MOONPA_DE421)
```

The default ITRF model is the [`iau2006b`](@ref), but other approximations are also [available](@ref iers_models). If one was interested in the Moon's PA440 axes, a similar function named [`add_axes_pa440!`](@ref) is available.

!!! note 
    In order to use `IERS` associated reference frames, EOP must be loaded. 
    See also [`Orient.prepare_eop`](@ref), [`Orient.init_eop`](@ref).

    ```julia 
    eop_url = "https://datacenter.iers.org/data/csv/finals2000A.data.csv"
    # Download eop file 
    download(eop_url, "data.csv")

    # EOP data shall be transformed to FrameTransformations compatible format:
    Orient.prepare_eop("data.csv", "iau2000a")

    # Once transformed, before creating any IERS-associated frame the previously created 
    # file must be loaded:
    Orient.init_eop("iau2000a.eop.dat")
    ```

!!! note 
    To correctly use the [`add_axes_pa421!`](@ref) function, the frame system must contain an ephemeris provider that has loaded the necessary PCK kernels with the DE421 data.

### Using the `FrameSystem`

Now let's assume we have a spacecraft in orbit about the Moon in the `LME2000` axes, and we want to compute its position in the `GCRF` and with respect to both the Moon and the Earth surface. We then start by registering our spacecraft as an [updatable point](@ref updatable_points):

```@setup usingHiFi1
using FrameTransformations 
@point Moon 301 
@axes LME2000 5

FRAMES = FrameSystem{2, Float64}()

add_axes_inertial!(FRAMES, LME2000)
add_point_root!(FRAMES, Moon, LME2000)
```

```@repl usingHiFi1
@point SC -1_900_000
add_point_updatable!(FRAMES, SC, Moon, LME2000)

```

Therefore our final `FrameSystem` results in:

```@setup finalFramesHiFi
using FrameTransformations
using Ephemerides, Downloads

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
url_eop = "https://datacenter.iers.org/data/csv/finals2000A.data.csv"

let
    eopfile = "iau2000a"
    Orient.prepare_eop(Downloads.download(url_eop), eopfile)
    Orient.init_eop(eopfile * ".eop.dat")
end;

eph = EphemerisProvider([Downloads.download(url_spk), Downloads.download(url_pck)])

tpc = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc"
iau = load(TPC(Downloads.download(tpc)));

@axes GCRF 1 GeocentricCelestialReferenceFrame
@axes ICRF 1 InternationalCelestialReferenceFrame

@axes IAU_EARTH 3 
@axes IAU_MOON 4
@axes LME2000 5
@axes ITRF 6 
@axes MOONPA_DE421 31006 

@point SSB 0 
@point EMB 3
@point Sun 10 
@point Earth 399 
@point Moon 301 

@point SC -1_900_000

FRAMES = FrameSystem{3, Float64}(eph)

add_axes_inertial!(FRAMES, GCRF)

add_point_root!(FRAMES, SSB, GCRF)
add_point_ephemeris!(FRAMES, EMB)
add_point_ephemeris!(FRAMES, Earth)
add_point_ephemeris!(FRAMES, Moon)
add_point_ephemeris!(FRAMES, Sun)

add_axes_bcrtod!(FRAMES, IAU_EARTH, Earth, iau);
add_axes_bcrtod!(FRAMES, IAU_MOON, Moon, iau);
add_axes_bci2000!(FRAMES, LME2000, Moon, iau);

add_axes_itrf!(FRAMES, ITRF, GCRF)
add_axes_pa421!(FRAMES, MOONPA_DE421)

add_point_updatable!(FRAMES, SC, Moon, LME2000)

e = Epoch("2020-01-01T12:45:30.0 TDB");
x = [2274.0, 0.0, 0.0, 0.0, sqrt(4904.87/2274.0), 0.0];
update_point!(FRAMES, SC, x, e)
```

```@repl finalFramesHiFi 
FRAMES 
```

To begin exploiting our newly completed frame system, we assumed that our `SC` is on a circular equatorial orbit around the Moon at an altitude of 500 km, so that at the epoch `e`, the spacecraft state vector is updated as:

```@repl finalFramesHiFi
e = Epoch("2020-01-01T12:45:30.0 TDB");
x = [2274.0, 0.0, 0.0, 0.0, sqrt(4904.87/2274.0), 0.0];

update_point!(FRAMES, SC, x, e)
```

!!! note 
    The timescale used for the `Epoch` must be the same used in the `FrameSystem`.

At this point, we are completely free to compute the state of the spacecraft with respect to any other point registered in the frame system with respect to any known set of axes at the epoch `e`. For example, the state in the `LME2000` is: 

```@repl finalFramesHiFi
vector6(FRAMES, Moon, SC, LME2000, e)
```

The spacecraft state with respect to the Earth in the `GCRF`, `IAU_EARTH` and `ITRF` is instead:

```@repl finalFramesHiFi
vector6(FRAMES, Earth, SC, GCRF, e)
vector6(FRAMES, Earth, SC, IAU_EARTH, e)
vector6(FRAMES, Earth, SC, ITRF, e)
```

while the position with respect to the Moon in the `IAU_MOON` and `PA421` axes:

```@repl finalFramesHiFi
vector3(FRAMES, Moon, SC, IAU_MOON, e)
vector3(FRAMES, Moon, SC, MOONPA_DE421, e)
```

These last examples are intended to show how easily the state of a spacecraft with respect to any other body can be retrieved by properly leveraging the [`FrameSystem`](@ref) and the high-level routines provided by this package. 