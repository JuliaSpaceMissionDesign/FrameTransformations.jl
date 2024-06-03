using FrameTransformations

CR3BP = FrameSystem{2, Float64}()

add_axes_root!(CR3BP, :InertialAx, 1)
add_point_root!(CR3BP, :EMBc, 1, 1)

using ReferenceFrameRotations

f(t) = angle_to_dcm(t, :Z)

add_axes_rotating!(CR3BP, :SynodicAx, 2, 1, f)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
