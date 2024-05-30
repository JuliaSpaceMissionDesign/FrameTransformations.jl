```@meta
EditURL = "../t00_frames.jl"
```

# [Frame System Overview](@id tutorial_00_frames)
_This example was generated on 2024-05-30T14:30:06.151._

The core object of `FrameTransformations` is the [`FrameSystem`](@ref), which provides
the capability to compute relative position, orientation and their time derivatives up to
order 3 (jerk), between standard and user-defined point and axes. It works by creating two
separate graphs that silently store and manage all the parent-child relationships between
the user-registered axes and points, in the form of `FramePointNode` and `FrameAxesNode`.

These two objects define two precise entities:
- **Axes**: defines an orientation in space. These are related each other by means of a
  [`Rotation`](@ref) transformation which relate one axes to a parent axes in
  a certain time interval.
- **Points**: defines a location in space. These are related each other by
  means of a `Translation` transformation which relate one point to a parent point in a
  particular axes in a certain time interval.

!!! note
    A single [`FrameSystem`](@ref) instance simultaneously handles both the axes and
    point graphs, regardless of what the user registers in it. For instance, if no
    points are added, the point graph will remain empty.

Additionally, any node can have several childs, each with different transformations with
respect to the parent node. However, they shall be **registered** within the
[`FrameSystem`](@ref) before being used in a transformation or as parents of other nodes.

````@example t00_frames
using FrameTransformations #hide
using Tempo #hide
````

## Basic Constructors
The creation of a generic [`FrameSystem`](@ref) requires the definition of the maximum
desired transformation order and of its `DataType`, which in most applications is a `Float64`.
The transformation order is always one greater than the maximum desired time derivative.
For instance, if the user only desires to compute position and velocity components (i.e.,
order 1 time-derivative), the transformation order to be used is 2. Thus, the maximum
allowed transformation order is 4.

In this example, we highlight the most basic way to initialise a [`FrameSystem`](@ref):

````@example t00_frames
F = FrameSystem{2, Float64}()
````

From this example, you can see that within the frame system there are both point and axes
graphs. However, at the moment they are completely empty since the graph was just created.

Each [`FrameSystem`](@ref) object is assigned a reference timescale that is used to perform
computations with epochs and to parse ephemeris files. The default timescale is the
`BarycentricDynamicalTime`, however, the user is free to select the most suited timescale
for his applications. In this example, we set the `InternationalAtomicTime` as the reference scale.

````@example t00_frames
F = FrameSystem{2, Float64, InternationalAtomicTime}()
````

## Ephemerides Support

In certain scenarios, the transformations require usage of binary ephemeris kernels, e.g.,
the JPL's DE440 files. To support this applications, this package has an interface relying
on [JSMDInterfaces.jl](https://github.com/JuliaSpaceMissionDesign/JSMDInterfaces.jl)
`AbstractEphemerisProvider`s. Currently, this package is shipped with extension for the
following two ephemeris readers:
* [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl)
* [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl)

Once the desired ephemeris provider is created, it can be used to register points or axes.
In this example we begin loading an old DE421 kernels to pass to the ephemeris reader.

````@example t00_frames
using Ephemerides, Downloads

url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(url));

F = FrameSystem{2, Float64}()
````

Before registering any node, a set of root axes and a root node shall be anyway registered.

````@example t00_frames
add_axes_icrf!(F)
add_point_root!(F, :SSB, 0, 1)
````

Points from the `EphemerisProvider` can be now registered.

````@example t00_frames
add_point_ephemeris!(F, eph, :Sun, 10)
add_point_ephemeris!(F, eph, :EMB, 3)
````

Here the parent point will be inferred from the ephemeris.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
