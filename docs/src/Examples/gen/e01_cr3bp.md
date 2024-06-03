```@meta
EditURL = "../e01_cr3bp.jl"
```

# [Use Case: CR3BP](@id example_01_cr3bp)
_This example was generated on 2024-05-31T10:00:28.409._

The power of the [`FrameSystem`](@ref) is its capability to handle axes
transformations and point translations of both high-accuracy and simplified
models. The use-case here presented includes the case of the Circular-Restricted
Three-Body Problem (CR3BP) rotating frame transformation handling.

In particular, when dealing with the [CR3BP](https://orbital-mechanics.space/the-n-body-problem/circular-restricted-three-body-problem.html),
mission analysis are used to exploit non-dimensional, rotating coordinates
to express the equations of motion and perform the computations.

In this tutorial, we create a [`FrameSystem`](@ref) to handle transformations
within the Earth-Moon CR3BP, which is characterized by a mass ratio of
approximately `μ = 0.012`. We start off by creating a frame system without any
ephemeris provider, since we are using a simplified model.

````@example e01_cr3bp
using FrameTransformations

CR3BP = FrameSystem{2, Float64}()
````

As always, the first step requires the definition of the root axes and points.
In this case, we use the a generic set of _inertial axes_ and the Earth-Moon
Barycenter (EMB).

````@example e01_cr3bp
add_axes_root!(CR3BP, :InertialAx, 1)
add_point_root!(CR3BP, :EMBc, 1, 1)
````

We now proceed to add our synodic axes: in the CR3BP these are uniformly
rotating with respect to the `InertialAx` about the Z-axis. Therefore, we
leverage rotating axes.

````@example e01_cr3bp
using ReferenceFrameRotations

f(t) = angle_to_dcm(t, :Z)

add_axes_rotating!(CR3BP, :SynodicAx, 2, 1, f)
````

Note that there is no need to specify the rotation derivatives, as they'll
be computed by  automatic differentiation via the
[ForwardDiff](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package.
For performace-critical transformations, however, it is reccomended to manually
define these derivatives.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

