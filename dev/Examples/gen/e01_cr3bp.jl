using FrameTransformations

CR3BP = FrameSystem{2, Float64}()

@axes InertialAx 1 InertialFrame

add_axes_inertial!(CR3BP, InertialAx)

@point EMBc 1 EarthMoonBarycenterCr3bp

add_point_root!(CR3BP, EMBc, InertialAx)

using ReferenceFrameRotations

f(t) = angle_to_dcm(t, :Z)

@axes SynodicAx 2 SynodicFrame

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)

@point Spacecraft -1_900_000

add_point_updatable!(CR3BP, Spacecraft, EMBc, SynodicAx)

μ = 0.012
xL4 = [1/2-μ, sqrt(3)/2, 0.0, 0.0, 0.0, 0.0]
t = 0.8

update_point!(CR3BP, Spacecraft, xL4, t)

vector6(CR3BP, EMBc, Spacecraft, SynodicAx, t)

vector6(CR3BP, EMBc, Spacecraft, InertialAx, t)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
