# [Points Creation and Translations](@id tutorial_02_points)

Similarly to [axes](@ref tutorial_01_axes), `FrameTransformations` also provides the capability to define custom and standard reference points (e.g., the Solar System Barycenter) and arbitrarily connect them through the [`FrameSystem`](@ref) In turn, this allows the computation of the relative position and its derivatives (up to order 3) between any two registered points and express it in any known set of axes. 

At the time being, the following types of points are supported:
- **Root point**: it is the root of the point graph.
- **Fixed points**: are those whose positions have a constant offset with respect their parent point in a given set of axes.
- **Dynamical points**: the position of these points depends only on time and is computed through custom user-defined functions.
- **Ephemeris points**: are those whose state-vector is retrieved from binary SPK kernels (e.g., DE440) that are loaded within the [`FrameSystem`](@ref).
- **Updatable points**: differently from all the other classes, the state vector for updatable points must be manually updated at a given epoch before it can be used in any transformation at the same epoch.

!!! note
    This package provides a dedicated function to register each type of supported points.

## Graph Initialisation

In this section we will display how to create a frame system to compute generic points transformation. Differently from the axes graph, each register point is also associated to a set of axes. Hence, this tutorial assumes the reader to already be familiar with the different types of axes and their definitions. First of all, we need to load both this package and an ephemeris reader. The latter will be used to retrieve the positions of the planets from the binary SPK kernels. In this example, we will use our own [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package and download the kernels from NAIF's website.

```@setup init 
using FrameTransformations
```

```@repl init
G = FrameSystem{2, Float64}()
```

To initialise the point graph, we first need to define a root point. This, in turn, must be associated to an arbitrary set of axes. Therefore, we begin by definining a generic `SatelliteFrame`, here considered as inertial, and then register a root point, called `SpacecraftCenter` in our graph. Similarly, to axes, the [`@point`](@ref) macro is used to define an acronym, an ID and a name of each point that we wish to register in the system. If a name is not provided, a default one is used. 

```@setup rootPoint 
using FrameTransformations 

G = FrameSystem{2, Float64}()

```
 
A root point can be registered using the [`add_point_root!`](@ref) function: 

```@repl rootPoint
@axes SATF 1 SatelliteFrame 

add_axes_inertial!(G, SATF)

@point SC -10000 SpacecraftCenter

add_point_root!(G, SC, SATF)
```

!!! note 
    For standard applications, it is good practice that the points's IDs are as in agreement with NAIF's numbering system. This becomes mandatory to properly read JPL's SPK kernels.

!!! note 
    The frame system uses an integer system based on the user-defined IDs to compute the transformations between axes and points. The name and acronym of the point are only used as aliases to provide a user-friendly interface to the transformations and do not have any other meaning. 

We can now see that our axes and point graphs are populating themselves:

```@setup pastRootPoint
using FrameTransformations 

G = FrameSystem{2, Float64}()
@axes SATF 1 SatelliteFrame 
@point SC -10000 SpacecraftCenter

add_axes_inertial!(G, SATF)
add_point_root!(G, SC, SATF)
```

```@repl pastRootPoint
G
```

## Fixed Points

Fixed points have a constant relative position vector with respect to their parent points in a given set of axes. Similarly to fixed-offset axes, these points are fixed w.r.t. their parents but might be moving with respect to others. 

In this example, we use the [`add_point_fixed!`](@ref) function to register the location of an antenna and two solar panels, which are generally fixed in the satellite body-fixed frame. To do so, we define a position offset in the form of a 3-elements vector with respect to the `SpacecraftCenter`.

```@repl pastRootPoint

@point SACL -10101 SolarArrayCenterLeft
@point SACR -10102 SolarArrayCenterRight
@point Antenna -10001

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixed!(G, SACL, SC, SATF, sa_offset_left)
add_point_fixed!(G, SACR, SC, SATF, sa_offset_right)
add_point_fixed!(G, Antenna, SC, SATF, an_offset)
```

As a result the graph is now populated with the new points and we can finally compute their relative positions and velocities with the proper transformation functions: 

```@setup pastFixedPoints
using FrameTransformations 

G = FrameSystem{2, Float64}()

@axes SATF 1 SatelliteFrame 
@point SC -10000 SpacecraftCenter
@point SACL -10101 SolarArrayCenterLeft
@point SACR -10102 SolarArrayCenterRight
@point Antenna -10001

add_axes_inertial!(G, SATF)
add_point_root!(G, SC, SATF)

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixed!(G, SACL, SC, SATF, sa_offset_left)
add_point_fixed!(G, SACR, SC, SATF, sa_offset_right)
add_point_fixed!(G, Antenna, SC, SATF, an_offset)
```

```@repl pastFixedPoints
vector3(G, SACL, SC, SATF, 0.0)
vector6(G, Antenna, SACR, SATF, 10.0)
```

As expected, since these points are fixed, the relative velocity vector is null.

## Dynamical Points
Dynamical points are generic time-dependent points whose position vector (and optionally its derivatives)
are only function of time. However, differently from ephemeris points, their position is computed through 
user-defined functions.

```@repl pastFixedPoints
@point TimeDependantAppendage -10003

fun(t) = [cos(t), sin(t), 0]

add_point_dynamical!(G, TimeDependantAppendage, SACL, SATF, fun)

vector6(G, TimeDependantAppendage, SC, SATF, π/3)
```

!!! note 
    To avoid allocations, `fun` should return a static array.

Similarly to rotating-axes, if the user only provides the function to compute the relative position, the remaining derivatives are automatically retrievied via automatic differentiation of `fun`. On the other hand, if those functions are specified, they must return a single vector that stacks all the components. For instance, for the second order derivative of `fun`, the function should return a 9-elements vector containing the relative position, velocity and acceleration. For example: 

```@repl pastFixedPoints
@point TimeDependantAppendage2 -10004

fun(t) = [cos(t), sin(t), 0]
dfun(t) = [cos(t), sin(t), 0, -sin(t), cos(t), 0]

add_point_dynamical!(G, TimeDependantAppendage2, SACL, SATF, fun, dfun)

vector6(G, TimeDependantAppendage2, SC, SATF, π/3)
```

We can again see that the results are in agreement with the previous example. For more details, consult the [`add_point_dynamical!`](@ref) documentation.

## [Updatable Points](@id updatable_points)

Updatable points are a class of point whose states at a given epoch must be manually updated 
before any other computation at the same epoch can occur. They can be inserted in the 
computational graphs as follows:

```@repl pastFixedPoints
@point UA -10002 UpdatableAppendage

add_point_updatable!(G, UA, SC, SATF)
```

If we now call a transformation involving this point, an error will be thrown because we have 
not registered any state for this point. To do so, we use the [`update_point!`](@ref) function and then 
evaluate the relative position: 

```@setup pastUpdatable
using FrameTransformations 

G = FrameSystem{2, Float64}()

@axes SATF 1 SatelliteFrame 
@point SC -10000 SpacecraftCenter
@point SACL -10101 SolarArrayCenterLeft
@point SACR -10102 SolarArrayCenterRight
@point Antenna -10001
@point UA -10002 UpdatableAppendage

add_axes_inertial!(G, SATF)
add_point_root!(G, SC, SATF)

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixed!(G, SACL, SC, SATF, sa_offset_left)
add_point_fixed!(G, SACR, SC, SATF, sa_offset_right)
add_point_fixed!(G, Antenna, SC, SATF, an_offset)
add_point_updatable!(G, UA, SC, SATF)
```

```@repl pastUpdatable 
ua_pos = [0.0, -1.0, 0.0]
update_point!(G, UA, ua_pos, 0.0)
vector3(G, Antenna, UA, SATF, 0.0)
```

Note that in the previous example, only the position has been updated but the current frame
system is of order two! Therefore, in this case, calling [`vector6`](@ref) will give an error since 
the computational graph is order-sentitive.

To correct that behaviour, also the higher order shall be updated:

```@repl pastUpdatable
update_point!(G, UA, [1.0, -1.0, 0.0, 0.0, 0.0, 0.0], .0)
vector6(G, Antenna, UA, SATF, 0.0)
```

!!! note  
    Updatable points do not store any history of the updated states. Meaning that each time 
    a state is updated at a different epoch, the information of the older epochs is completely lost.

## Ephemeris Points

Ephemeris points are a type of time-dependent points whose position and higher-order derivatives are retrieved from a binary SPK ephemeris kernel. However, differently from all other points, in this case the set of axes is automatically inferred from those contained in the ephemeris kernels. In case such set is not yet registered in the frame system, an error will be thrown. 

!!! note 
    To properly compute the position of these points, the [`FrameSystem`](@ref) object must contain an ephemeris provider that has loaded the necessary kernels. Additionally, in this case the ID of the registered points must match the ID contained in the SPK kernels. 

In this example, we define a new frame system `F` and give it an ephemeris provider that has loaded the DE421 SPK kernel, containing the position of the major planets and/or their barycenters of the Solar System.

```@repl init
using Ephemerides, Downloads
spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(spk))
F = FrameSystem{2, Float64}(eph)
```

```@setup ephemPoint 
using FrameTransformations
using Ephemerides, Downloads
spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(spk))
F = FrameSystem{2, Float64}(eph)
```

We now register some points that are stored in the kernels using the [`add_point_ephemeris!`](@ref) function: 

```@repl ephemPoint
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
```

Notice that this function does not requires the parent point. Indeed, the parent is automatically set to those contained in the descriptors of the ephemeris kernels. For instance, in the DE421, the Earth-Moon Barycenter (EMB) is defined with respect to the SSB, which the frame system automatically uses as parent for the EMB. Similarly, the EMB is the default parent point for the Earth.

!!! warning 
    If a parent point is not specified and the point in the kernels has not yet been registered, an error is thrown. 

Finally, we can retrieve the transformation data as usual: 

```@setup pastEphemPoint
using FrameTransformations
using Ephemerides, Downloads
spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(spk))
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
```

```@repl pastEphemPoint
vector6(F, EMB, SSB, ICRF, 1000.0)

vector3(F, Earth, SSB, ICRF, 0.0)
```