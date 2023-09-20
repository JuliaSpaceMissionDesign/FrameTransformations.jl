# [Points graphs creation and handling: Spacecraft Geometry Example](@id tutorial_01_points)


```julia
using FrameTransformations
```

`FrameTransformations` provides the possibility to create graphs of points completely defined 
by the user and use unique capabilities of the `FrameSystem` to handle different types of them. 
At the time being, the following points types are allowed:

- **Root point**: it is the root of the graph and has an assigned set of axes. All the other points are then assumed to belong to the same axes.
- **Fixed points**: are those whose positions have a constant `offset` with respect their `parent` points in the given set of `axes`.
- **Updatable points**: differently from all the other classes, the state vector for updatable points (expressed in the set of input `axes`) shall be manually updated before being used for any other computations.
- **Dynamical points**:  for them the state vector for these points depends only on time and is computed through the custom functions provided by the user.
- **Ephemeris points**: is intended for points whose state-vector is read from ephemeris kernels.

This tutorial will cover only the first four types since a dedicated tutorial is available for
ephemeris points within the `Frames` tutorial.

## Graph creation

Let's assume we want to create a points computational graph whose points are assigned
w.r.t. to a dummy frame called `SatelliteFrame`.

There are two things to be defined to create the computational graph: the `order` 
of the graph (i.e. 1 for position, 2 for position and velocity, ...) and the `timescale`
in which time is represented within the graph.


```julia
fs = FrameSystem{2, Float64, BarycentricDynamicalTime}()
```




    FrameSystem{2, Float64, BarycentricDynamicalTime, JSMDUtils.NullEphemerisProvider}(
      eph: JSMDUtils.NullEphemerisProvider(),
      points: EMPTY
      axes: EMPTY
    )




We can see that, within the frame system there are both `points` and `axes` graphs. 
In this case, at the moment, they are completely empty.

## Register the root point

The first step in the construction of the graph is the definition of the root point axes model.
A new axes can be created by means of the [`@axes`](@ref) macro. 


```julia
# Create the new axes
@axes SATF 1 SatelliteFrame

# Register the axes 
add_axes_inertial!(fs, SATF)
```

Once the root point axes are defined the root point can be also defined. 
Generally speaking, there are always two ways of working with points within `FrameTransformations`:
using aliases or not.


```julia
# Define Spacecraft center alias 
@point SC -10000 SpacecraftCenter

# Register the new point as a root 
add_point_root!(fs, SC, SATF) # Using aliases 
# add_point_root!(fs, :SC, -10000, SATF) # Not using aliases
```

Now, the computational graph starts to be populated by the root point and root axes: 


```julia
fs
```




    FrameSystem{2, Float64, BarycentricDynamicalTime, JSMDUtils.NullEphemerisProvider}(
      eph: JSMDUtils.NullEphemerisProvider(),
      points: 	 
    	 SC
    	 
      axes: 	
    	SATF
    	
    )




## Register fixed points

To register fixed points, the method `add_point_fixed!` can be called. For example let's assume
that we want to add the location of an antenna and a some solar solar panels:


```julia
# Create the new points
@point SACL -10101 SolarArrayCenterLeft
@point SACR -10102 SolarArrayCenterRight
@point Antenna -10001

# Define offsets 
sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixed!(fs, SACL, SC, SATF, sa_offset_left)
add_point_fixed!(fs, SACR, SC, SATF, sa_offset_right)
add_point_fixed!(fs, Antenna, SC, SATF, an_offset)
```

As a result the graph is not populated with new points:


```julia
fs
```




    FrameSystem{2, Float64, BarycentricDynamicalTime, JSMDUtils.NullEphemerisProvider}(
      eph: JSMDUtils.NullEphemerisProvider(),
      points: 	 
    	 SC
    	  ├── SACL 
    	  ├── SACR 
    	  ├── Antenna 
    	 
      axes: 	
    	SATF
    	
    )




We can now use the `vectorX` methods to compute the desired quantities out of the graph:


```julia
# We can evaluate the positions of the different points in the graph:
vector3(fs, SC, SACL, SATF, 0.0)
```




    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
     1.0
     0.0
     0.0




```julia
# and we can also combine paths
vector3(fs, SACL, SACR, SATF, 0.0)
```




    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
     -2.0
     -0.0
     -0.0



## Register updatable points

Updatable points are a class of point whose states shall be manually updated before the other
computations. They can be inserted in the computational graphs as follows:


```julia
# Let's create a new point, which will be "updatable"
@point UA -10002 UpdatableAppendage

# Register the new point 
add_point_updatable!(fs, UA, SC, SATF)
```

If we call any `vectorX` method without performing the update of the desired state, we'll
have an error.  Therefore, we shall first call the `update_point!` method and the evaluate the computational graph:


```julia
update_point!(fs, UA, [0.0, -1.0, 0.0], 0.0)
vector3(fs, Antenna, UA, SATF, 0.0)
```




    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
     -0.0
     -1.0
      1.0



Note that in the previous example, only the position has been updated but the current frame
system is of order two! Therefore, in this case, calling `vector6` will give an error since 
the computational graph is order-sentitive.

To correct that behaviour, also the higher order shall be updated:


```julia
update_point!(fs, UA, [1.0, -1.0, 0.0, 0.0, 0.0, 0.0], 0.0)
vector6(fs, Antenna, UA, SATF, 0.0)
```




    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
      1.0
     -1.0
      1.0
     -0.0
     -0.0
     -0.0



## Register dynamical points

It is possible to register also time-dependent points (these are not generally ephemeris ones,
but could be). To do that, first a time-dependent function shall be defined and then 
the `add_point_dynamical!` can be used to insert the point.


```julia
# Create the new point 
@point TimeDependantAppendage -10003

# Define how the point behaves in time 
fun(t::T) where T = [cos(t), sin(t), 0]

# Register the point 
add_point_dynamical!(fs, TimeDependantAppendage, UA, SATF, fun)
```

Note that it is possible to create parent-child relations between any kind of points.

Once registered, the usual `vectorX` methods could be called. Note that, differently from the
updatable points case, if the user-defined function returns a state vector which is smaller 
than the one associated to the order of the computational graph,the other orders are computed via autodiff.


```julia
fs
```




    FrameSystem{2, Float64, BarycentricDynamicalTime, JSMDUtils.NullEphemerisProvider}(
      eph: JSMDUtils.NullEphemerisProvider(),
      points: 	 
    	 SC
    	  ├── SACL 
    	  ├── SACR 
    	  ├── Antenna 
    	  ├── UA 
    	   ├── TimeDependantAppendage 
    	 
      axes: 	
    	SATF
    	
    )




Let's not evaluate the computational graph at a new point and see how it behaves:


```julia
update_point!(fs, UA, [1.0, -1.0, 0.0, 0.0, 0.0, 0.0], π/3)
vector6(fs, SC, TimeDependantAppendage, SATF, π/3)
```




    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
      1.5
     -0.1339745962155614
      0.0
     -0.8660254037844386
      0.5000000000000001
      0.0


