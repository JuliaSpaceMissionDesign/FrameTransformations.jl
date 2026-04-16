## [Frame System](@id frames_api)

```@docs 
FrameSystem
order 
FrameTransformations.timescale 
```

### Points

```@docs 
has_point
points_graph
points_alias
```

### Axes

```@docs 
has_axes 
axes_graph
axes_alias
```

### Directions

```@docs 
has_direction
directions
```

## [Rotations](@id rotation_api)

```@docs 
Rotation
Base.inv
Translation
```

## [Transformations](@id transformations_api)

### [Points](@id points_transform_api)

```@docs 
vector3
vector6
vector9 
vector12 
```

### [Axes](@id axes_transform_api)

```@docs 
rotation3
rotation6
rotation9
rotation12
```

### [Directions](@id directions_transform_api) 

```@docs 
direction3 
direction6 
direction9
direction12
```

## [Compiled Fast-Path](@id compiled_api)

The compiled fast-path provides zero-overhead, AD-transparent callables that bypass the 
`FunctionWrapper` type-erasure barrier used internally by `FrameSystem`. Use these when you 
need full inlining, custom AD-backend support (e.g., Mooncake, Zygote), or maximum 
performance in hot loops such as ODE right-hand sides.

### Types

```@docs 
CompiledRotation
CompiledTranslation
CompiledDirection
```

### Constructors

```@docs 
compile_rotation
compile_translation
compile_direction
```
