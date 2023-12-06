# [Use Case: CR3BP](@id tutorial_03_cr3bp)

The power of the [`FrameSystem`](@ref) is its capability to handle axes transformations
and point translations of both high-accuracy and simplified models. The use-case here presented includes the case of the Circular-Restricted Three-Body Problem (CR3BP) rotating frame transformation handling.

In particular, when dealing with the [CR3BP](https://orbital-mechanics.space/the-n-body-problem/circular-restricted-three-body-problem.html), mission analysis are used to exploit non-dimensional, rotating coordinates to express the equations of motion and perform the computations. 

In this tutorial, we create a [`FrameSystem`](@ref) to handle transformations within the Earth-Moon 
CR3BP, which is characterized by a mass ratio of approximately `μ = 0.012`. We start off by creating a frame system without any ephemeris provider, since we are using a simplified model.

```@setup init 
using FrameTransformations
```

```@repl init
CR3BP = FrameSystem{2, Float64}()
```

As always, the first step requires the definition of the root axes and points. In this case, we use the a generic set of _inertial axes_ and the Earth-Moon Barycenter (EMB).

```@setup rootCR3BP
using FrameTransformations 
CR3BP = FrameSystem{2, Float64}()
```

```@repl rootCR3BP
@axes InertialAx 1 InertialFrame 

add_axes_inertial!(CR3BP, InertialAx)

@point EMBc 1 EarthMoonBarycenterCr3bp

add_point_root!(CR3BP, EMBc, InertialAx)
```

We now proceed to add our synodic axes: in the CR3BP these are uniformly rotating with respect to the `InertialAx` about the Z-axis. Therefore, we leverage the [rotating axes](@ref rot_axes) type: 

```@setup postRootCR3BP
using FrameTransformations 

@axes InertialAx 1 InertialFrame 
@point EMBc 1 EarthMoonBarycenterCr3bp

CR3BP = FrameSystem{2, Float64}()

add_axes_inertial!(CR3BP, InertialAx)
add_point_root!(CR3BP, EMBc, InertialAx)
```

```@repl postRootCR3BP
using ReferenceFrameRotations

f(t) = angle_to_dcm(t, :Z) 

@axes SynodicAx 2 SynodicFrame 

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)
```

Note that there is no need to specify the rotation derivatives, as they'll be computed by 
automatic differentiation via the [ForwardDiff](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package. For performace-critical transformations, however, it is reccomended to manually define these derivatives.

Now, let's assume we have our spacecraft. Most likely, its states will be expressed in the synodic frame.
In this case, we leverage [updatable points](@ref updatable_points), since we desire to manually update its state at each time.

```@setup postSynodicCR3BP
using FrameTransformations 
using ReferenceFrameRotations

@axes InertialAx 1 InertialFrame 
@axes SynodicAx 2 SynodicFrame 
@point EMBc 1 EarthMoonBarycenterCr3bp

CR3BP = FrameSystem{2, Float64}()

add_axes_inertial!(CR3BP, InertialAx)
add_point_root!(CR3BP, EMBc, InertialAx)

f(t) = angle_to_dcm(t, :Z) 

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)
```

```@repl postSynodicCR3BP
@point Spacecraft -1_900_000

add_point_updatable!(CR3BP, Spacecraft, EMBc, SynodicAx)
```

We know assume that at `t = 0.8` our spacecraft is at L4, therefore we update its state accordingly:

```@setup postSpacecraftCR3BP
using FrameTransformations 
using ReferenceFrameRotations

@axes InertialAx 1 InertialFrame 
@axes SynodicAx 2 SynodicFrame 
@point EMBc 1 EarthMoonBarycenterCr3bp
@point Spacecraft -1_900_000
CR3BP = FrameSystem{2, Float64}()

f(t) = angle_to_dcm(t, :Z) 

add_axes_inertial!(CR3BP, InertialAx)
add_point_root!(CR3BP, EMBc, InertialAx)

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)
add_point_updatable!(CR3BP, Spacecraft, EMBc, SynodicAx)
```

```@repl postSpacecraftCR3BP
μ = 0.012
xL4 = [1/2-μ, sqrt(3)/2, 0.0, 0.0, 0.0, 0.0]
t = 0.8

update_point!(CR3BP, Spacecraft, xL4, t)
```

Finally we can retrieve the spacecraft state in both the synodic as well as the inertial axes: 

```@setup postUpdateCR3BP
using FrameTransformations 
using ReferenceFrameRotations

@axes InertialAx 1 InertialFrame 
@axes SynodicAx 2 SynodicFrame 
@point EMBc 1 EarthMoonBarycenterCr3bp
@point Spacecraft -1_900_000
CR3BP = FrameSystem{2, Float64}()

f(t) = angle_to_dcm(t, :Z) 

add_axes_inertial!(CR3BP, InertialAx)
add_point_root!(CR3BP, EMBc, InertialAx)

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)
add_point_updatable!(CR3BP, Spacecraft, EMBc, SynodicAx)

μ = 0.012
xL4 = [1/2-μ, sqrt(3)/2, 0.0, 0.0, 0.0, 0.0]
t = 0.8

update_point!(CR3BP, Spacecraft, xL4, t)

t = 0.8
```

```@repl postUpdateCR3BP
vector6(CR3BP, EMBc, Spacecraft, SynodicAx, t)
vector6(CR3BP, EMBc, Spacecraft, InertialAx, t)
```