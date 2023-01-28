# [Axes graphs creation and handling](@id tutorial_05_axes)


```julia
using Basic, ReferenceFrameRotations
```

`Basic` provides the possibility to create graphs of axes completely defined by the user 
and use unique capabilities of the `FrameSystem` to handle different types of them.

At the time being, the following axes types are allowed:
- **Inertial axes**: these are the only ones which can be used as root. Can be also other than root axes.
- **Fixed offset axes**:  have a constant orientation with respect to their `parent` axes, represented by `dcm`, a Direction Cosine Matrix (DCM).
- **Rotating axes**:  the orientation of these axes depends only on time and is computed through the custom functions provided by the user
- **Computable axes**:  differ from rotating axes because they are computed through two vectors that are defined within the frame system itself. Computable axes are the equivalent of SPICE's parameterized two-vector frames.
- **Projected axes**:  the orientation of these axes depends only on time and is computed through the custom functions provided by the user. Projected axes are similar to rotating axis, except that all the positions, velocity, etc ... are rotated by the 0-order rotation (i.e. the derivatives of the rotation matrix are null, despite the rotation depends on time).
- **Ephemeris axes**: these are constructed by looking for euler angles and their derivatives needed by the rotation from the ephemeris kernels loaded in `frames`.

## Graph creation

Then, let's assume we want to create an axes computational graph whose axes are assigned w.r.t.
the `ICRF` frame.

There is actually two other things to be assigned to create the computational graph: the `order` 
of the graph, i.e. if position, velocity, acceleration,... shall be computed and the timescale
in which time is represented within the graph. In this case, the latter is embedded in the 
constructor since we are using an ephemeris provider.


```julia
# Load ephemeris to memory
eph = ephem_load(
    CalcephProvider, 
    [
        "/home/andrea/Documents/Kernels/spk/de440.bsp", 
        "/home/andrea/Documents/Kernels/pck/moon_pa_de440_200625.bpc"
    ]
)

# Create the graph
G = FrameSystem{3, Float64}(eph)
```


    FrameSystem{3, Float64, BarycentricDynamicalTime, CalcephProvider}(
      eph: CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x0000000003ee44c0)),
      points: EMPTY
      axes: EMPTY
    )



We can see that, within the frame system there are both `points` and `axes` graphs. In
this case, at the moment, they are completely empty.

## Register the root axes

Once the computational graph has been created, the root axes shall be registered. These shall
be **inertial** axes, without any offset. In this example, let us consider the `ICRF`.


```julia
# Create the new axes 
@axes ICRF 1 InternationalCelestialReferenceFrame 

# Register the new axes in the graph as root axes
add_axes_inertial!(G, ICRF)
```

Once the root axes are registered any other axis type can be added.

## Register fixed-offset axes

Fixed offset axes are a simple way to represent a constant orientation w.r.t. their parent axes.

!!! note
    While inertial axes do not rotate with respect to the star background, fixed offset axes are only 
    constant with respect to their parent axes, but might be rotating with respect to some other 
    inertial axes.


```julia
# Create new axis 
@axes FO1 2 

# Register 
add_axes_fixedoffset!(G, FO1, ICRF, angle_to_dcm(π/4, :Z))
```

Now we can call any `rotationX` method to get the rotation matrix:


```julia
R = rotation6(G, ICRF, FO1, 125*86400.0)
```


    Rotation{2, Float64}(([0.7071067811865476 0.7071067811865475 0.0; -0.7071067811865475 0.7071067811865476 0.0; 0.0 0.0 1.0], [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]))


Then, to apply the rotation to a vector, we can exploit the `Base.*` overloads available in
`Basic` and simply call:


```julia
R * [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```


    6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
      0.7071067811865476
     -0.7071067811865475
      0.0
      0.0
      0.0
      0.0


## Register rotating axes

Rotating axes are, generic, time-dependant, non-inertial axes. In order to register this axes
of this kind in the graph, first a function (and optionally their derivatives) shall be defined.

The function shall return a  Direction Cosine Matrix (DCM). Let us call it `fun`.

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed via automatic differentiation.
Here `δⁿ` are the functions that return the DCM and the n-th order derivative.


```julia
# Define new axes
@axes RotAx 3

# Define the transformation
fun(t) = angle_to_dcm(t, :Z)'

# Register the new axes 
add_axes_rotating!(G, RotAx, FO1, fun)
```

Now we can call any `rotationX` method to get any desired rotation matrix:


```julia
rotation6(G, ICRF, RotAx, π/4)
```


    Rotation{2, Float64}(([1.0 1.0146536357569526e-17 0.0; 1.0146536357569526e-17 1.0 0.0; 0.0 0.0 1.0], [-1.0146536357569526e-17 -1.0 0.0; 1.0 1.0146536357569526e-17 0.0; 0.0 0.0 0.0]))


## Register computable axes

Computable axes are a kind of _time dependant axes_. In this case, differently from the 
rotating axes, the axes and their derivatives are computed through two time-dependant vectors
that are extracted from the registered ephemeris. These axes are the equivalent of SPICE's
two-vector frames.


```julia
# Create new axes 
@axes SunFrame 4

# Create points to be used
@point SSB 0 SolarySystemBarycenter 
@point Sun 10 

# Register the root point 
add_point_root!(G, SSB, ICRF)

# Register other points 
add_point_ephemeris!(G, Sun, SSB)

# Create principal direction
v1 = ComputableAxesVector(10, 0, 1)

# Create secondary direction
v2 = ComputableAxesVector(10, 0, 2)

# Register the new axes
add_axes_computable!(G, SunFrame, ICRF, v1, v2, :XY)
```


```julia
rotation6(G, ICRF, SunFrame, 0.0)
```


    Rotation{2, Float64}(([-0.930764538489918 -0.34524125434920866 -0.12035717753850385; 0.36466636541206515 -0.8528478114547875 -0.3737232297221213; 0.026378321152114293 -0.39173854391949264 0.91969837304468], [4.875905229140774e-9 -1.1403314091867927e-8 -4.997003351253353e-9; 1.2447345420279773e-8 4.583247734107517e-9 1.686583449794906e-9; -3.0651424534885396e-11 7.168470364139439e-11 3.141268419372045e-11]))


## Register projected axes

Projected axes are a particular type of inertial axes. In this case the rotation is built by
means of a time dependant function `f(t)`. However, all the derivatives of `f(t)` are assumed
to be zero. This axes type is usually used to build True-of-Date (TOD) axes sets.


```julia
# Define new axes
@axes ProjAx 5

# Define the transformation
fun(t) = angle_to_dcm(t, :Z)'

# Register the new axes 
add_axes_projected!(G, ProjAx, ICRF, fun)
```


```julia
rotation6(G, ICRF, ProjAx, 1.0)
```


    Rotation{2, Float64}(([0.5403023058681398 -0.8414709848078965 0.0; 0.8414709848078965 0.5403023058681398 0.0; 0.0 0.0 1.0], [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]))


## Register ephemeris axes

Ephemeris axes a are a type of time-dependent axes which are build by means of Euler angles
contained within an ephemeris file (`.bpc`). These are used to build Lunar frames or other
high-precision axes.


```julia
# Define axes 
# These shall have the ID of the ephemeris kernel
@axes MoonDummy 31008

# Register new axes 
add_axes_ephemeris!(G, MoonDummy, :ZXZ)
```


```julia
rotation6(G, ICRF, MoonDummy, 86400.0)
```


    Rotation{2, Float64}(([0.6219858271328162 0.7077417415402774 0.33501531028811493; -0.7827094708541836 0.5741658942616172 0.24020701509723086; -0.022349834027802825 -0.41162501521201655 0.9110792126761739], [-2.08357886027015e-6 1.5281551322405242e-6 6.400225290645924e-7; -1.6557201080619335e-6 -1.8840092793849879e-6 -8.917888467237793e-7; -3.897629745766924e-10 -4.729986090501684e-10 -2.2326181369018686e-10]))

