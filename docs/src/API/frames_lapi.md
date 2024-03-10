# [Low-level API](@id low_frames_api)

## Axes

```@docs
FrameTransformations.AbstractFrameAxes
FrameTransformations.FrameAxesNode
FrameTransformations.axes_name
FrameTransformations.axes_id

FrameTransformations.build_axes

FrameTransformations.ComputableAxesProperties

```

## Points

```@docs
FrameTransformations.AbstractFramePoint
FrameTransformations.FramePointNode

FrameTransformations.point_name
FrameTransformations.point_id 

FrameTransformations.build_point

FrameTransformations._get_comp_axes_vector3
FrameTransformations._get_comp_axes_vector6
FrameTransformations._get_comp_axes_vector9
FrameTransformations._get_comp_axes_vector12

```

## Two Vectors

```@docs
FrameTransformations.twovectors_to_dcm 
FrameTransformations.twovectors_to_δdcm 
FrameTransformations.twovectors_to_δ³dcm 
FrameTransformations.twovectors_to_δ²dcm

FrameTransformations._twovectors_to_dcm 
FrameTransformations._two_vectors_to_rot6
FrameTransformations._two_vectors_to_rot9
FrameTransformations._two_vectors_to_rot12
FrameTransformations._two_vectors_basis

```