# Points graphs creation and handling


```julia
using Basic
```

`Basic` provides the possibility to create graphs of points completely defined by the user 
and use unique capabilities of the `FrameSystem` to handle different types of them.

At the time being, the following points types are allowed:
- **Root point**: it is the root of the graph and has an assigned set of axes. All the other points are then assumed to belong to the same axes.
- **Ephemeris points**: is intended for points whose state-vector is read from ephemeris kernels.
- **Fixed points**: are those whose positions have a constant `offset` with respect their `parent` points in the given set of `axes`.
- **Updatable points**: differently from all the other classes, the state vector for updatable points (expressed in the set of input `axes`) shall be manually updated before being used for any other computations.
- **Dynamical points**:  for them the state vector for these points depends only on time and is computed through the custom functions provided by the user.

## Graph creation

Then, let's assume we want to create a points computational graph whose points are assigned
w.r.t. to a dummy frame called `SatelliteFrame`.

There is actually two other things to be assigned to create the computational graph: the `order` 
of the graph, i.e. if position, velocity, acceleration,... shall be computed and the timescale
in which time is represented within the graph.


```julia
G = FrameSystem{2, Float64, BarycentricDynamicalTime}()
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider}(
      eph: Basic.Ephemeris.NullEphemerisProvider(),
      points: EMPTY
      axes: EMPTY
    )



We can see that, within the frame system there are both `points` and `axes` graphs. In
this case, at the moment, they are completely empty.

## Register the root point


Let's now register our root node! For that purpose, we can exploit the `add_point_root!` function:


```julia
# Create axes
@axes SATF 1 SatelliteFrame

# Create the root point
@point SC -10000 Spacecraft

# Register the axes 
add_axes_inertial!(G, SATF)

# Register the root point 
add_point_root!(G, SC, SATF)
```

Now, our computational graph starts to be populated by the root point and the root axes:


```julia
G
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider}(
      eph: Basic.Ephemeris.NullEphemerisProvider(),
      points: 	 
    	 SC
    	 
      axes: 	
    	SATF
    	
    )



## Register fixed points

To register fixed points, the method `add_point_fixed!` can be called. For example let's assume
that we want to add the location of an antenna and a solar panel:


```julia
# Create the new points
@point SolarArrayCenter -10001
@point AntennaCenter -10002

# Define offsets 
sa_offset = [1.0, 1.0, 0.0]
an_offset = [-1.0, 1.0, 0.0]

add_point_fixed!(G, SolarArrayCenter, SC, SATF, sa_offset)
add_point_fixed!(G, AntennaCenter, SC, SATF, an_offset)
```

Then, the graph is populated by the new points:


```julia
G
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider}(
      eph: Basic.Ephemeris.NullEphemerisProvider(),
      points: 	 
    	 SC
    	  ├── SolarArrayCenter 
    	  ├── AntennaCenter 
    	 
      axes: 	
    	SATF
    	
    )



We can now use the `vectorX` methods to compute the desired quantities:


```julia
vector3(G, SolarArrayCenter, AntennaCenter, SATF, 0.0)
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
@point UpdatableAppendage -10003

# Register the new point 
add_point_updatable!(G, UpdatableAppendage, SC, SATF)
```

Now if we call any `vectorX` method without performing the update of the desired state, we'll
have an error:


```julia
vector3(G, AntennaCenter, UpdatableAppendage, SATF, 0.0)
```

    UpdatablePoint with NAIFId -10003 has not been updated at time 0.0 for order 1
    Stacktrace:
     [1] _compute_vector3(point::Basic.Frames.FramePointNode{2, Float64, 6}, t::Float64)
       @ Basic.Frames ~/Documents/Gitlab/Astronaut/Basic/src/Frames/transform.jl:273
     [2] _compute_vector3
       @ ~/Documents/Gitlab/Astronaut/Basic/src/Frames/transform.jl:259 [inlined]
     [3] _get_vector3_backwards(frame::FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider, 6}, t::Float64, path::Vector{Int64})
       @ Basic.Frames ~/Documents/Gitlab/Astronaut/Basic/src/Frames/transform.jl:234
     [4] _compute_vector3(frame::FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider, 6}, t::Float64, axesid::Int64, path::Vector{Int64})
       @ Basic.Frames ~/Documents/Gitlab/Astronaut/Basic/src/Frames/transform.jl:195
     [5] vector3(frame::FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider, 6}, from::AntennaCenterPoint, to::UpdatableAppendagePoint, axes::SatelliteFrameAxes, t::Float64)
       @ Basic.Frames ~/Documents/Gitlab/Astronaut/Basic/src/Frames/transform.jl:182
     [6] top-level scope
       @ ~/Documents/Gitlab/Astronaut/Basic/docs/src/Tutorials/t_04_points.ipynb:1


Therefore, we shall first call the `update_point!` method and the evaluate the computational graph:


```julia
update_point!(G, UpdatableAppendage, [1.0, -1.0, 0.0], 0.0)
vector3(G, AntennaCenter, UpdatableAppendage, SATF, 0.0)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
      2.0
     -2.0
     -0.0


Note that in the previous example, only the position has been updated! Therefore calling `vector6`
will give an error since the computational graph is order-sentitive.

To correct that behaviour, also the higher order shall be updated:


```julia
update_point!(G, UpdatableAppendage, [1.0, -1.0, 0.0, 0.0, 0.0, 0.0], 0.0)
vector6(G, AntennaCenter, UpdatableAppendage, SATF, 0.0)
```


    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
      2.0
     -2.0
     -0.0
     -0.0
     -0.0
     -0.0


## Register dynamical points

It is possible to register also time-dependent points (these are not generally ephemeris ones,
but could be). To do so, the following steps could be performed:


```julia
# Create the new point 
@point TimeDependantAppendage -10004

# Define how the point behaves in time 
fun(t::T) where T = [cos(t), sin(t), 0]

# Register the point 
add_point_dynamical!(G, TimeDependantAppendage, UpdatableAppendage, SATF, fun)
```

Note that it is possible to create parent-child relations between any kind of points.
Once registered, the usual `vectorX` methods could be called. Note that if the user-defined 
function returns a state vector which is smaller than the one of the computational graph,
the other orders are computed via autodiff.


```julia
G
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, Basic.Ephemeris.NullEphemerisProvider}(
      eph: Basic.Ephemeris.NullEphemerisProvider(),
      points: 	 
    	 SC
    	  ├── SolarArrayCenter 
    	  ├── AntennaCenter 
    	  ├── UpdatableAppendage 
    	   ├── TimeDependantAppendage 
    	 
      axes: 	
    	SATF
    	
    )



Let's not evaluate the computational graph at a new point and see how it behaves:


```julia
update_point!(G, UpdatableAppendage, [1.0, -1.0, 0.0, 0.0, 0.0, 0.0], π/3)
vector6(G, SC, TimeDependantAppendage, SATF, π/3)
```


    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
      1.5
     -0.1339745962155614
      0.0
     -0.8660254037844386
      0.5000000000000001
      0.0


## Register ephemeris points

Note that the previous graph  not point to any ephemeris provider, therefore _ephemeris points_
cannot be registered there.

To do that, we need to load, through an ephemeris provider, some ephemeris kernels. 
In this case we use the `CalcephProvider`:


```julia
# Load ephemeris to memory
eph = ephem_load(CalcephProvider, ["/home/andrea/Documents/Kernels/spk/de440.bsp"])

# Create the graph
G = FrameSystem{2, Float64}(eph)

# Create points
@point SSB 0 SolarSystemBarycenter
@point EMB 3 EarthMoonBarycenter 
@point Sun 10
@point Earth 399

# Create default axes 
@axes ICRF 1 InternationalCelestialReferenceFrame

# Register the ICRF axes 
add_axes_inertial!(G, ICRF)

# Register the root node 
add_point_root!(G, SSB, ICRF)

# Register the other nodes 
add_point_ephemeris!(G, EMB)
add_point_ephemeris!(G, Earth)
add_point_ephemeris!(G, Sun)
```

Now the points graph is fully operative and we can use `vectorX` methods to get the required quantities:


```julia
vector3(G, Sun, Earth, ICRF, 0.0)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
     -2.6499033677425094e7
      1.3275741733833946e8
      5.755671847053819e7


Note that the parents of the ephemeris points are assigned via `ephem_available_points` method. 
In fact, the comutational graph of points already has some parent-child relations:


```julia
G
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, CalcephProvider}(
      eph: CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x000000000694aae0)),
      points: 	 
    	 SSB
    	  ├── EMB 
    	   ├── Earth 
    	  ├── Sun 
    	 
      axes: 	
    	ICRF
    	
    )


