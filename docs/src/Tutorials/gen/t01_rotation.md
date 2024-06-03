```@meta
EditURL = "../t01_rotation.jl"
```

# [Rotations](@id tutorial_01_rotation)
_This example was generated on 2024-05-31T11:33:00.402._

Before diving into the creation of the axes graph, it is worth highlighting that transformations
that express the relative orientation or its time-derivatives between two generic set of
axes are represented by a [`Rotation`](@ref) object, which stores a Direction Cosine Matrix
(DCM) and its derivatives. This package leverages the already available
[ReferenceFrameRotations.jl](https://github.com/JuliaSpace/ReferenceFrameRotations.jl)
to define the DCM objects.

A time-fixed rotation between two axes and its derivative can then be expressed as follows:

````@example t01_rotation
using FrameTransformations
using ReferenceFrameRotations

dcm  = angle_to_dcm(π/3, :Z)
δdcm = DCM(0I)

R = Rotation(dcm, δdcm)
````

````@example t01_rotation
R[1]
````

````@example t01_rotation
R[2]
````

A rotation object is returned by all the rotation functions that are applied to the `FrameSystem`.
It provide overloads to the basic algebraic operations so that multiplication and inversions
can be efficiently computed leveraging the properties of rotation matrixes.

For example, to rotate a generic vector `v`, we can simply do:

````@example t01_rotation
v = [1., -6., 3., 0., 5., 0]
R*v
````

The inverse can instead be taken as:

````@example t01_rotation
inv(R)
````

See the [Rotation API](@ref rotation_api) for more information on this object.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

