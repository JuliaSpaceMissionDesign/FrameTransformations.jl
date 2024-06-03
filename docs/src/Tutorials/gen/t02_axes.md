```@meta
EditURL = "../t02_axes.jl"
```

# [Axes](@id tutorial_01_axes)
_This example was generated on 2024-05-31T11:33:00.405._

To compute relative orientations, `FrameTransformations` provides the capability to define
custom and standard reference axes (e.g., the ITRF) and arbitrarily connect them through
the [`FrameSystem`](@ref) In turn, this allows the computation of the relative orientation
and its derivatives (up to order 3) between any two registered axes.

At the time being, the following types of axes are supported:
- **Inertial axes**: these are the only ones which can be used as root axes to initialise
  the axes graph.
- **Fixed offset axes**: they have a constant orientation with respect to their parent axes.
- **Rotating axes**: the orientation of these axes depends only on time and is computed t
  through the custom functions provided by the user.
- **Ephemeris axes**: these are constructed by extracting the Euler rotation angles and their
  derivatives from the binary PCK kernels that are loaded within the [`FrameSystem`](@ref).

!!! note
    This package provides a dedicated function to register each type of supported axes.
    Additionally, higher-level functions to automatically register standard astronomical
    reference axes are also provided, e.g., [`add_axes_ecl2000!`](@ref).

## Graph Initialisation

In this section we will display how to create a frame system to compute generic axes rotation.
First of all, we need to load both this package and an ephemeris reader.
The latter will be used to compute the orientation of the Moon's Principal Axes (PA) 440,
 whose Euler angles are defined in binary PCK kernels and to retrieve the positions of the
planets. In this example, [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl)
package and download the kernels from NAIF's website.

````@example t02_axes
using FrameTransformations
using Ephemerides

url_pck = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc";
url_spk = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";

const EPH = EphemerisProvider([download(url_spk), download(url_pck)])
const F = FrameSystem{3, Float64}()
````

To initialise the axes graph, a set of root axes must be initially registered.
These will serve as the uppermost node of the graph and have no parents, meaning their
orientation is not specified. Only inertial axes can be used as root axes of the
[`FrameSystem`](@ref).

 In this example, we will use the `ICRF` as our base root inertial axes.

````@example t02_axes
add_axes_root!(F, :ICRF, AXESID_ICRF)
````

Once a set of root axes has been registered, any other type of axes can be added to the system.

!!! note
    For standard applications, it is good practice that the axes's IDs are as in agreement
    with NAIF's numbering system. A list of IDs for the most common axes is provided in
    this package.

!!! note
    The frame system uses an integer system based on the user-defined IDs to compute
    the transformations between axes and points.

Inertial axes are those that are fixed with respect to the star background.
They are the only ones that can be used as root axes in the frame system but can also be
defined through a relative orientation with respect to another set of inertial axis.

## [Inertial Axes](@id ine_axes)

In this example, we register the `GCRF` as a set of inertial axes with respect to
the `ICRF`. We assume that the two frames are equivalent, thus:

````@example t02_axes
using ReferenceFrameRotations
using LinearAlgebra

fun(t) = DCM(1.0I)
add_axes_inertial!(F, :GCRF, AXESID_GCRF, AXESID_ICRF, fun)

R = rotation6(F, AXESID_ICRF, AXESID_GCRF, 1.0)
````

````@example t02_axes
R[1]
````

````@example t02_axes
R[2]
````

Since it is an inertial frame, the time derivative of the rotation is null.

## [Fixed-offset Axes](@id fox_axes)

Fixed-offset axes have a constant orientation with respect to their parent axes in time.
We previously saw that inertial axes can also be used to define axes with a fixed orientation
 with respect to their parents. However, while inertial axes do not rotate with respect to
the star background, fixed offset axes are only constant with respect to their parent axes,
but might be rotating with respect to some other inertial axes.

In this example, we register `FOX` as a set of axes with a fixed rotation of `π/4` around
the Z-axis with respect to the `ICRF`.

````@example t02_axes
rot = angle_to_dcm(π/4, :Z)

add_axes_fixedoffset!(F, :FOX, 2, AXESID_ICRF, rot)
````

The state rotation matrix can then be obtained as:

````@example t02_axes
R = rotation6(F, 1, 2, 86400)
````

````@example t02_axes
R[1]
````

````@example t02_axes
R[2]
````

Since `FOX` has a constant orientation with respect to the `ICRF`, the time derivative of
the rotation matrix `R[2]` is, in fact, null. For further information see the
[`add_axes_fixedoffset!`](@ref) documentation.

## [Rotating Axes](@id rot_axes)

Rotating axes are generic, time-dependant, non-inertial axes. In order to register this
kind of axes, a function (and optionally its derivatives) that expresses the relative
orientation of this axes must be defined. This function shall return a Direction Cosine
Matrix (DCM), available from [ReferenceFrameRotations.jl](https://github.com/JuliaSpace/ReferenceFrameRotations.jl).

````@example t02_axes
fun(t) = angle_to_dcm(-t, :Z)

add_axes_rotating!(F, :ROX, 3, AXESID_ICRF, fun)
````

If we now compute the orientation between the `FOX` and `ROX` at `π/4` we obtain an identity
rotation, since the orientation of `ROX` is directed in the opposite direction of `FOX`.

````@example t02_axes
R = rotation6(F, 2, 3, π/4)
````

````@example t02_axes
R[1]
````

Notice that, although we only provided a function that expresses the relative orientation,
the frame system has automatically computed its time-derivative via Automatic Differentiation
(AD) of `fun`.

````@example t02_axes
R2 = rotation6(F, 1, 3, π/4)
````

````@example t02_axes
R2[2]
````

This becomes particularly useful for rapid prototyping or when the manual differentiation
requires a lot of time. The functions for higher-order derivatives, must return the original
DCM and its derivatives up to their orders. For example:

````@example t02_axes
using JSMDUtils.Math

fun(t) = angle_to_dcm(-t, :Z)
dfun(t) = (angle_to_dcm(-t, :Z), Math.angle_to_δdcm([-t, -1], :Z))

add_axes_rotating!(F, :ROX2, 4, AXESID_ICRF, fun, dfun)

R2 = rotation6(F, 1, 3, π/4)
````

````@example t02_axes
R2[2]
````

We can see the results are in agreement with the previous example.
For more details, see [`add_axes_rotating!`](@ref) documentation.

## Ephemeris Axes

Ephemeris axes a are a type of time-dependent axes which are build by means of Euler angles
contained within a binary PCK ephemeris kernel. For example, in practice these are used
to express the orientation of high-accuracy Lunar body-fixed frames (i.e., the Principal
Axes) or the Earth's ITRF.

!!! note
    To properly compute the orientation of these axes, the ephemeris provider used
    must contain the necessary PCK kernels.
    Additionally, in this case the ID of the registered axes must match the ID
    contained in the PCK kernels.

In this example, the ephemeris provider `EPH` has loaded the DE421
PCK kernel containing the orientation of the Moon's Principal Axes (PA421). NAIF's system
has assigned to such set of axes the ID `31006`. If a different ID was assigned to the
`MoonPA`, the function would have thrown an error.

The function also requires the user to specify the rotation sequence to convert the Euler
angles to a proper rotation matrix.

````@example t02_axes
add_axes_ephemeris!(F, EPH, :MOONPA, 31006, :ZXZ)

R = rotation6(F, 1, 31006, 86400.0)
````

````@example t02_axes
R[1]
````

````@example t02_axes
R[2]
````

For further information see the [`add_axes_ephemeris!`](@ref) documentation.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

