```@meta
EditURL = "../e02_hifi.jl"
```

# [Use Case: High Fidelity](@id example_02_hifi)
_This example was generated on 2024-05-30T14:30:06.242._

Once the general structure of the [`FrameSystem`](@ref) is understood, we can pass to
a use case in which we want to build and exploit our frame system to perform computations
in a high-fidelity environment.

### Frame system setup

In this example, we plan on using ephemeris data to retrieve accurate positions
of the planets and the orientation of certain reference frames. Therefore, we
create an ephemeris provider object leveraging our own
[`Ephemerides.jl`](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package
and use it to generate a frame system instance:

````@example e02_hifi
using FrameTransformations
using Ephemerides
using LinearAlgebra
using ReferenceFrameRotations

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp";

EPH = EphemerisProvider([download(url_spk), download(url_pck)])

FRAMES = FrameSystem{3, Float64}()
````

Once the graph is created, we assign the `ICRF` (International Celestial Reference
Frame) as our set of inertial root-axes:

````@example e02_hifi
add_axes_icrf!(FRAMES)
````

These axes are practically coincident with the `GCRF`.

In this scenario, we will be working within the Cislunar environment, therefore
we will need the major bodies that influence this dynamic regime, i..e, the Earth,
the Moon and the Sun. To do so, we also define the Solar System Barycenter (SSB) and
the Earth-Moon Barycenter (EMB) as the ephemeris data of the remaining bodies is
expressed with respect to those.

For this example, we will assume the SSB is our root point:

````@example e02_hifi
add_point_root!(FRAMES, :SSB, 0, 1)
add_point_ephemeris!(FRAMES, EPH, :EMB, 3)
add_point_ephemeris!(FRAMES, EPH, :Sun, 10)
add_point_ephemeris!(FRAMES, EPH, :Earth, 399)
add_point_ephemeris!(FRAMES, EPH, :Moon, 301)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

