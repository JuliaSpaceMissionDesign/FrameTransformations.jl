# [Axes Creation and Rotations](@id tutorial_01_axes)

To compute relative orientations, `FrameTransformations` provides the capability to define custom and standard reference axes (e.g., the ITRF) and arbitrarily connect them through the [`FrameSystem`](@ref) In turn, this allows the computation of the relative orientation and its derivatives (up to order 3) between any two registered axes. 

At the time being, the following types of axes are supported:
- **Inertial axes**: these are the only ones which can be used as root axes to initialise the axes graph. 
- **Fixed offset axes**: they have a constant orientation with respect to their parent axes.
- **Rotating axes**: the orientation of these axes depends only on time and is computed through the custom functions provided by the user
- **Computable axes**: they are computed through two vectors that are defined within the frame system itself. Computable axes are the equivalent of SPICE's parameterized two-vector frames.
- **Projected axes**: the orientation of these axes depends only on time and is computed through the custom functions provided by the user. Projected axes are similar to rotating axis, except that all the positions, velocity, etc ... are rotated by the 0-order rotation (i.e. the derivatives of the rotation matrix are null, despite the rotation depends on time).
- **Ephemeris axes**: these are constructed by extracting the Euler rotation angles and their derivatives from the binary PCK kernels that are loaded within the [`FrameSystem`](@ref).

!!! note
    This package provides a dedicated function to register each type of supported axes. Additionally, higher-level functions to automatically register standard astronomical reference axes are also provided, e.g., [`add_axes_eclipj2000!`](@ref).

## Rotations

Before diving into the creation of the axes graph, it is worth highlighting that transformations that express the relative orientation or its time-derivatives between two generic set of axes are represented by a [`Rotation`](@ref) object, which stores a Direction Cosine Matrix (DCM) and its derivatives. This package leverages the already available [ReferenceFrameRotations.jl](https://github.com/JuliaSpace/ReferenceFrameRotations.jl) to define the DCM objects. 

A time-fixed rotation between two axes and its derivative can then be expressed as follows: 

```@setup initAxes
using FrameTransformations 
```

```@repl initAxes
using ReferenceFrameRotations

dcm  = angle_to_dcm(π/3, :Z)
δdcm = DCM(0I)

R = Rotation(dcm, δdcm)

R[1]

R[2]
```

A rotation object is returned by all the rotation functions that are applied to the frame system. It provide overloads to the basic algebraic operations so that multiplication and inversions can be efficiently computed leveraging the properties of rotation matrixes. For example, to rotate a generic vector `v`, we can simply do: 


```@setup rotation 
using FrameTransformations
using ReferenceFrameRotations

dcm  = angle_to_dcm(π/3, :Z)
δdcm = DCM(0I)

R = Rotation(dcm, δdcm)
```

```@repl rotation 
v = [1., -6., 3., 0., 5., 0]
R*v
```

The inverse can instead be taken as: 

```@repl rotation 
inv(R)
```

See the [Rotation API](@ref rotation_api) for more information on this object.


## Graph Initialisation

In this section we will display how to create a frame system to compute generic axes rotation. First of all, we need to load both this package and an ephemeris reader. The latter will be used to compute the orientation of the Moon's Principal Axes (PA) 440, whose Euler angles are defined in binary PCK kernels and to retrieve the positions of the planets. In this example, we will use our own [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package and download the kernels from NAIF's website.

```@setup init 
using FrameTransformations
```

```@repl init
using Ephemerides, Downloads

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

eph = EphemerisProvider([Downloads.download(url_spk), Downloads.download(url_pck)])

G = FrameSystem{3, Float64}(eph)
```

To initialise the axes computational graph, a set of root axes must be initially registered. These will serve as the uppermost node of the graph and have no parents, meaning their orientation is not specified. Only inertial axes can be used as root axes of the [`FrameSystem`](@ref). 

Before registering the axes, the [`@axes`](@ref) macro is used to define an acronym, an ID and a name of each set of axes that we wish to register in the system. If a name is not provided, a default one is used. In this example, we will use the `ICRF` as our base root inertial axes.

```@setup axesRoot 
using FrameTransformations
G = FrameSystem{3, Float64}()
```

```@repl axesRoot 
@axes ICRF 1 InternationalCelestialReferenceFrame

add_axes_inertial!(G, ICRF)
```

Once a set of root axes has been registered, any other type of axes can be added to the system.

!!! note 
    For standard applications, it is good practice that the axes's IDs are as in agreement with NAIF's numbering system. A list of IDs for the most common axes is provided in the Orient submodule.

!!! note 
    The frame system uses an integer system based on the user-defined IDs to compute the transformations between axes and points. The name and acronym of the axes are only used as aliases to provide a user-friendly interface to the transformations and do not have any other meaning. For instance, one could register a set of rotating axes named `ICRF`.


## Inertial Axes

Inertial axes are those that are fixed with respect to the star background. They are the only ones that can be used as root axes in the frame system but can also be defined through a relative orientation with respect to another set of inertial axis. 

In this example, we register the `ECLIPJ2000` as a set of inertial axes with respect to the `ICRF`. Since the set of root axes has already been defined, all the future usages of the [`add_axes_inertial!`](@ref) function require a parent set of axes and a DCM with the relative orientation.

```@setup axesInertial
using FrameTransformations 
G = FrameSystem{3, Float64}()
@axes ICRF 1 InternationalCelestialReferenceFrame
add_axes_inertial!(G, ICRF)
```

```@repl axesInertial

@axes ECLIPJ2000 17

add_axes_inertial!(G, ECLIPJ2000; parent=ICRF, dcm=DCM_ICRF_TO_ECLIPJ2000)

R = rotation6(G, ICRF, ECLIPJ2000, 10.0)
R[1]
R[2]
```

Notice that we leveraged the default DCM provided by Orient's to express the relative orientation between the `ICRF` and the `ECLIPJ2000`. For a set of default DCM exported by Orient, check the [API documentation](@ref orient_dcms). Additionally, since it is an inertial frame, the time derivative of the rotation is null. 


## Fixed-offset Axes

Fixed-offset axes have a constant orientation with respect to their parent axes in time. We previously saw that inertial axes can also be used to define axes with a fixed orientation with respect to their parents. However, while inertial axes do not rotate with respect to the star background, fixed offset axes are only constant with respect to their parent axes, but might be rotating with respect to some other inertial axes.

In this example, we register `FO1` as a set of axes with a fixed rotation of `π/4` around the Z-axis with respect to the `ICRF`.

```@setup axesFixedOffset
using FrameTransformations
using ReferenceFrameRotations
G = FrameSystem{3, Float64}()

@axes ICRF 1 InternationalCelestialReferenceFrame
add_axes_inertial!(G, ICRF)
```

```@repl axesFixedOffset 
@axes FO1 2

rot = angle_to_dcm(π/4, :Z)

add_axes_fixedoffset!(G, FO1, ICRF, rot)
```

```@setup axesFixedOffset2
using FrameTransformations
using ReferenceFrameRotations
G = FrameSystem{3, Float64}()

@axes ICRF 1 InternationalCelestialReferenceFrame
add_axes_inertial!(G, ICRF)

@axes FO1 2

rot = angle_to_dcm(π/4, :Z)

add_axes_fixedoffset!(G, FO1, ICRF, rot)
```

The state rotation matrix can then be obtained as: 

```@repl axesFixedOffset2
R = rotation6(G, ICRF, FO1, 86400)

R[1]

R[2]
```

Since `FO1` has a constant orientation with respect to the `ICRF`, the time derivative of the rotation matrix `R[2]` is infact null. For further information see the [`add_axes_fixedoffset!`](@ref) documentation.


## [Rotating Axes](@id rot_axes)

Rotating axes are generic, time-dependant, non-inertial axes. In order to register this kind of axes, a function (and optionally its derivatives) that expresses the relative orientation of this axes must be defined. This function shall return a Direction Cosine Matrix (DCM), available from the [ReferenceFrameRotations.jl](https://github.com/JuliaSpace/ReferenceFrameRotations.jl) package.

```@repl axesFixedOffset2
@axes RotAx 3 

fun(t) = angle_to_dcm(-t, :Z)

add_axes_rotating!(G, RotAx, FO1, fun)
```

```@setup axesRotating
using FrameTransformations
using ReferenceFrameRotations
G = FrameSystem{3, Float64}()

@axes ICRF 1 InternationalCelestialReferenceFrame
@axes FO1 2
@axes RotAx 3 

rot = angle_to_dcm(π/4, :Z)
fun(t) = angle_to_dcm(-t, :Z)

add_axes_inertial!(G, ICRF)
add_axes_fixedoffset!(G, FO1, ICRF, rot)
add_axes_rotating!(G, RotAx, FO1, fun)
```

If we now compute the orientation between the `ICRF` and `RotAx` at `π/4` we obtain an identity
rotation, since the orientation of `RotAx` is directed in the opposite direction of `FO1`.

```@repl axesRotating
R1 = rotation6(G, ICRF, RotAx, π/4)

R1[1]

R2 = rotation6(G, ICRF, RotAx, π/2)

R2[2]
```

Notice that, although we only provided a function that expresses the relative orientation, the frame system has automatically computed its time-derivative via Automatic Differentiation (AD) of `fun`. This becomes particularly useful for rapid prototyping or when the manual differentiation requires a lot of time. The functions for higher-order derivatives, must return the original DCM and its derivatives up to their orders. For example: 

```@repl axesFixedOffset2 
using JSMDUtils.Math

@axes RotAx2 4

fun(t) = angle_to_dcm(-t, :Z)
dfun(t) = (angle_to_dcm(-t, :Z), Math.angle_to_δdcm([-t, -1], :Z))

add_axes_rotating!(G, RotAx2, FO1, fun, dfun)

R2 = rotation6(G, ICRF, RotAx2, π/2)

R2[2]
```

We can see the results are in agreement with the previous example. For more details, consult the [`add_axes_rotating!`](@ref) documentation.

## Projected Axes

Projected axes are a particular type of inertial axes. In this case the rotation is built by
means of a time dependant function `f(t)`. However, all the derivatives of `f(t)` are assumed
to be zero. This axes type is usually used to build True-of-Date (TOD) axes sets. 

In this example, we illustrate this difference by registering two new set of axes with the same relative orientation with respect to the `ICRF`, one rotating and one projected. 

```@repl axesFixedOffset
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
```

As you can see, while the relative orientation `R[1]` is equal, the time-derivative of the projected-axes orientation is null. For further information see the [`add_axes_projected!`](@ref) documentation.

## Computable Axes

Computable axes are a kind of _time-dependant axes_. In this case, differently from the rotating axes, the axes and their derivatives are computed through two time-dependant vectors which are expressed using any type of point that is registered in the system. These axes are the equivalent of SPICE's two-vector frames.

In this example, we will register two ephemeris points, the Solar system barycenter and the Sun. For more information on how this operation is performed, see the [points tutorial](@ref tutorial_02_points). The two vectors that generate the set of computable axes are defined with the [`ComputableAxesVector`](@ref) object, by specifing the vector center and target point and its order, i.e., whether we are interested in the position, velocity or acceleration of that vector. A symbol is used to specify which direction the vectors have to align with. 

In this example, the axes are constructed with the X-axis parallel to the instantaneous SSB to Sun direction, whereas the secondary vector is chosen parallel to the SSB to Sun velocity vector (order 2). Then, the component of this vector orthogonal to the X-axis is used to create the Y-axis. 

```@setup ephemAxes
using FrameTransformations
using Ephemerides, Downloads

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

eph = EphemerisProvider([Downloads.download(url_spk), Downloads.download(url_pck)])

G = FrameSystem{3, Float64}(eph)

@axes ICRF 1 InternationalCelestialReferenceFrame
add_axes_inertial!(G, ICRF)
```

```@repl ephemAxes
@axes SunFrame 600

@point SSB 0 SolarSystemBarycenter 
@point Sun 10 SunPoint 

add_point_root!(G, SSB, ICRF)
add_point_ephemeris!(G, Sun)

v1 = ComputableAxesVector(Sun, SSB, 1)
v2 = ComputableAxesVector(Sun, SSB, 2)

add_axes_computable!(G, SunFrame, ICRF, v1, v2, :XY)

R = rotation6(G, ICRF, SunFrame, 0.0)
```

For further information see the [`add_axes_computable!`](@ref) documentation.

!!! note 
    The center and target point can either be specified with their integer IDs or their name aliases. 

## Ephemeris Axes

Ephemeris axes a are a type of time-dependent axes which are build by means of Euler angles contained within a binary PCK ephemeris kernel. For example, in practice these are used to express the orientation of high-accuracy Lunar body-fixed frames (i.e., the Principal Axes) or the Earth's ITRF.

!!! note 
    To properly compute the orientation of these axes, the [`FrameSystem`](@ref) object must contain an ephemeris provider that has loaded the necessary PCK kernels. Additionally, in this case the ID of the registered axes must match the ID contained in the PCK kernels. 

In this example, the ephemeris provider `eph` in the frame system `G` has loaded the DE421 PCK kernel containing the orientation of the Moon's Principal Axes (PA421). NAIF's system has assigned to such set of axes the ID `31006`. If a different ID was assigned to the `MoonPA`, the function would have thrown an error. A set of default axes IDs is also defined within the [Orient](@ref orient_axesid)'s submodule for ease of use.

The function also requires the user to specify the rotation sequence to convert the Euler angles to a proper rotation matrix.

```@repl ephemAxes
@axes MoonPA 31006

add_axes_ephemeris!(G, MoonPA, :ZXZ)

R = rotation9(G, ICRF, MoonPA, 86400.0)
```

For further information see the [`add_axes_ephemeris!`](@ref) documentation.
