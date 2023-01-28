# [Reference frames creation and transformations](@id tutorial_06_frames)


```julia
using Basic, ReferenceFrameRotations
```

## Introduction

`Basic` provides the capability to create, expand and efficienty differentiate Frame Systems. 
A frame system allows you to compute the relative position and orientation of a registered object
with respect to any other registered objects at any time. This includes time derivatives of motion
up to order 3 (jerk) and well integrated with `ForwardDiff`, allowing to compute partial 
derivatives of motion with respect to _any_ parameter.

`Frames` is a rather generic framework which can be extended and configured by the user. The Frame 
System, indeed, is the _ultimate goal_ of the `Basic` computational environment. Therefore, 
the aim of this tutorial is to provide an overview of the core functionalities of a Frame System 
together with some relevant informations relative to its structure and conceptual modeling, 
in such a way any user would be capable to understand, use and extend it in the proper way.

## The `FrameSystem`

The entry-point in the `Frames` computational environment, is the `FrameSystem` object. A 
`FrameSystem` manages a collection of points and axes in the form of `FramePointNode`s and 
`FrameAxesNode`s, respectively. These two kind of _nodes_, as you can see from the 
[Points Graphs Tutorial](@ref tutorial_04_points) and the [Axes Graphs Tutorial](@ref tutorial_05_axes),
are organized in a graph form and represent two precise entities:

- **Axes**: defines an orientation in space. These are related each other by means of a `Rotation` 
    transformation which relate one axes to a parent axes in a certain time interval.
- **Points**: defines a location in space. These are related each other by means of a `Translation`
    transformation which relate one point to a parent point in a particular axes in a certain 
    time interval.

Any node can have several transformations defining their orientation or position with 
respect to other `FrameAxesNode` or `FramePointNode`, each applicable during a particular 
time period. Moreover, as you have seen in the dedicated tutorials:

- Nodes can be **created** independently of each other (by means of `@axes` or `@point` macros).
- They shall be **registered** within the `FrameSystem` to be used, connecting with a parent by
    means of a transformation. This step can be usually performed by means of a dedicated method,
    i.e. `add_axes_xxxx!` or `add_point_yyyy!`, where `xxxx` and `yyyy` shall be specified 
    depending on the particular axes or point type.

!!! note
    Don't confuse the concept of frame of reference and axes! In fact, while in general these
    two concept are mixed together, they have a precise meaning within `Basic`:

    - Axes define coordinate systems, which are a mathematical concept;
    - A frame of reference is a physical concept related to the state of motion of something.

    In practice, we use a _coordinate system_ to specify a frame of reference because it is 
    convenient from a computational perspective but to each frame of reference there should 
    in theory correspond a family of an infinite numbers of coordinate systems comoving 
    together. 
    
    To clarify more, the velocity of a point with respect to a physical three-dimensional 
    body is an absolute concept (in Newtonian mechanics) but to express it as a set of 3 
    components we need to choose a coordinate system. A natural choice is to use a 
    coordinate system comoving with that body but it is also possible to use a coordinate 
    system rotating with respect to that body.

    It is normal to express quantities computed in a particular frame of reference in a 
    different coordinate system. For example, one might solve the motion of spacecraft in 
    an inertial frame of reference centered on the Earth, but express its velocity coordinates 
    along axes rotating with the Earth. This is not the same as computing the velocity 
    of the spacecraft in the Earth rotating system.

# Use Case: HiFi

Once the general structure of the `FrameSystem` is understood, we can pass to a use case in which
we want to build and exploit our computational graph collection to perform some computations.

### Creating an empty `FrameSystem`

The initial step is usually to create an empty `FrameSystem` object. Note that, depending
on the specific application, you might require or not the use of an ephemeris provider.
In this example we assume to be working in a high-fidelity environment and to load some ephemeris
via the `Calceph` interface:


```julia
# Load ephemeris to memory
eph = ephem_load(
    CalcephProvider, 
    [
        "/home/andrea/Documents/SpiceKernels/spk/de440.bsp", 
        "/home/andrea/Documents/SpiceKernels/pck/moon_pa_de440_200625.bpc"
    ]
)
```


    CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x00000000061ec350))


At this point we are ready to create the empty instance of our `FrameSystem`:


```julia
const FRAMES = FrameSystem{2, Float64}(eph)
```

    WARNING: redefinition of constant FRAMES. This may fail, cause incorrect answers, or produce other errors.



    FrameSystem{2, Float64, BarycentricDynamicalTime, CalcephProvider}(
      eph: CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x00000000061ec350)),
      points: EMPTY
      axes: EMPTY
    )



### Registering root axes & points

Once the graph has been created, the root axes shall be assigned. This shall be an **inertial** axes. 
In this case we consider the `GCRF`:


```julia
# Create the new axes 
@axes GCRF 1 GeocentricCelestialReferenceFrame

# Register the new axes in the graph as root axes
add_axes_inertial!(FRAMES, GCRF)
```

The axes here called `GCRF` is practically coincident with the `ICRF` axes and with the 
so called `ECI` frame. For that reason, let's create an alias of `GCRF` and call it `ICRF`. 
To do so we can register a zero-offset inertial axes (i.e. using identity matrix). Note that
the same can be achieved using a `fixed-offset` set of axes.


```julia
# Create the new axes 
@axes ICRF 2 InternationalCelestialReferenceFrame 

# Register the new axes in the graph as root axes
add_axes_inertial!(FRAMES, ICRF, parent=GCRF, dcm=DCM(1.0I))
```

Now, we assume we are working within the Cislunar environment, i.e. we want to insert the 
`Earth`, `Moon` and `Sun` bodies in the simulation considering the `Sun` the parent of the `Earth`, 
and the `Earth` the parent of the `Moon`. We want to use the `Earth` as root point as the `GCRF` 
is centered on it.


```julia
# Create the points 
@point Earth 399 
@point Sun 10 
@point Moon 301

# Register the root node 
add_point_root!(FRAMES, Earth, GCRF)
# Register the other nodes 
add_point_ephemeris!(FRAMES, Sun, Earth)
add_point_ephemeris!(FRAMES, Moon, Earth)
```

### Registering Earth and Lunar body axes

At this point, we want to be able to transform back and forth from body-fixed axes of the Earth 
and the Moon to the `GCRF`. The Earth and the Moon cases are actually _special cases_ because
the `Basic` computational environment provides the possibility to create both low and high 
precision frames. Let's start with the low precision one: `IAU_EARTH` and `IAU_MOON`.

To register such frames, we first need to parse a `.tcp` file with the required constants:


```julia
iau = Basic.load(TPC("/home/andrea/Documents/SpiceKernels/pck/pck00010.tpc"));
```

And then create the new axes. Remember, there shall be associated to a `ICRF` frame from 
definition. 

There are actually different axes which can be created starting from IAU Rotational Elements, 
in this case we are interested in Body-Centric Rotating, True-of-Date axes, or BCR-TOD (which 
are the `Basic` way of representing SPICE `IAU_XXX` axes).


```julia
# Create new axes 
@axes IAU_EARTH 3 
@axes IAU_MOON 4

# Register the new axes 
add_axes_bcrtod!(FRAMES, iau, Earth, IAU_EARTH, ICRF);
add_axes_bcrtod!(FRAMES, iau, Moon, IAU_MOON, ICRF);
```

    ┌ Warning: ignoring orient_axes_dd_icrf_to_bcr_tod_Earth, frame system order is less than 3
    └ @ Basic.Frames /home/andrea/.julia/dev/Basic/src/Frames/axes.jl:374
    ┌ Warning: ignoring orient_axes_dd_icrf_to_bcr_tod_Moon, frame system order is less than 3
    └ @ Basic.Frames /home/andrea/.julia/dev/Basic/src/Frames/axes.jl:374


Now let us insert also an inertial axes for the Moon, for convenience. This can be done using 
the IAU-based Body-Centric Inertial J2000 axes available in `Basic`, which a brings from `ICRF` to
the local equatorial plane of the body at J2000 epoch.


```julia
# Create new axes 
@axes LME2000 5

# Register the new axes 
add_axes_bci2000!(FRAMES, iau, Moon, LME2000, ICRF);
```

At this point, to conclude the graph, we want insert the high-precision Earth and Lunar body
axes. For the former, we can exploit `ITRF` axes for that purpose. For the latter, `PA440` ones.


```julia
# Create new axes 
@axes ITRF 6 
@axes MOONPA_DE440 31008 

# Register the new axes
add_axes_itrf!(FRAMES, ITRF, GCRF) # default IAUModel is iau2006b
add_axes_pa440!(FRAMES, MOONPA_DE440, ICRF)
```

    ┌ Warning: ignoring #154, frame system order is less than 3
    └ @ Basic.Frames /home/andrea/.julia/dev/Basic/src/Frames/axes.jl:374
    ┌ Warning: ignoring #155, frame system order is less than 4
    └ @ Basic.Frames /home/andrea/.julia/dev/Basic/src/Frames/axes.jl:374


### Using the `FrameSystem`

Now let's assume we have a spacecraft in orbit about the Moon in the `LME2000` axes, and we want
to compute the `GCRF` position of it, its position relative to the Moon surface and the Earth surface. 
Therefore, to fully exploit the capabilities of the `FrameSystem`, first the S/C shall be 
inserted within the computational graph. To do that, an **UpdatablePoint** can be used:


```julia
# Create the new point
@point SC -1_900_000

# Register the new point
add_point_updatable!(FRAMES, SC, Moon, LME2000)
```

Therefore the final `FrameSystem` results:


```julia
FRAMES
```


    FrameSystem{2, Float64, BarycentricDynamicalTime, CalcephProvider}(
      eph: CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x00000000061ec350)),
      points: 	 
    	 Earth
    	  ├── Sun 
    	  ├── Moon 
    	   ├── SC 
    	 
      axes: 	
    	GCRF
    	 ├── ICRF 
    	  ├── IAU_EARTH 
    	  ├── IAU_MOON 
    	  ├── LME2000 
    	 ├── ITRF 
    	 ├── MOONPA_DE440 
    	
    )



Now to exploit the Frame System, let's assume we have the `SC` state assigned as a point 
on a circular orbit at 500 km altitude, equatorial orbit about the Moon, so that, at epoch 
`e`, the state vector of the `SC` is:


```julia
e = Epoch("2020-01-01T12:45:30.0 TDB");
x = [2274.0, 0.0, 0.0, 0.0, sqrt(4904.87/2274.0), 0.0];

# Update! 
update_point!(FRAMES, SC, x, e)
```

!!! note 
    The timescale used for the `Epoch` **shall** be the same used in the `FrameSystem`.

At this point we are ready to get the desired (transformed) states:


```julia
# To GCRF 
vector6(FRAMES, Earth, SC, GCRF, e)
```


    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
     -284588.45268360287
     -271663.34366381646
      -78291.19908247433
           0.7300629920301382
           0.6818719329922416
           0.3078533603313829



```julia
# To IAU_MOON
vector3(FRAMES, Moon, SC, IAU_MOON, e)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
      1681.6697518585638
     -1530.7066491241687
         0.003318830657668824



```julia
# To MOONPA_DE440
vector3(FRAMES, Moon, SC, MOONPA_DE440, e)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
      1681.1972462850058
     -1531.2254587454092
         0.6430969025562603



```julia
# To IAU_EARTH
vector3(FRAMES, Earth, SC, IAU_EARTH, e)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
       20597.321286475617
     -392896.05374002125
      -78291.2054867196



```julia
# To ITRF
vector3(FRAMES, Earth, SC, ITRF, e)
```


    3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
     -391545.0010915323
       38505.3058372699
      -78300.22405992186


# Use Case: CR3BP
