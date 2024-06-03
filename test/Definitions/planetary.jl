using FrameTransformations
using ReferenceFrameRotations
using LinearAlgebra
using JSMDInterfaces.Bodies
using Test

function Bodies.body_rotational_elements(::Number, ::Val{-1})
    return -pi/2, pi/2, 0.0
end

frames = FrameSystem{4, Float64}()
@test_throws Exception add_axes_bci2000!(frames, :BCI, -1, -1)
@test_throws Exception add_axes_bcrtod!(frames, :BCR, -2, -1)
add_axes_icrf!(frames)

@test_nowarn add_axes_bci2000!(frames, :BCI, -1, -1)

R = rotation12(frames, :ICRF, :BCI, 12.345)

@test R[1] ≈ DCM(1.0I)
@test R[2] ≈ DCM(0.0I)
@test R[3] ≈ DCM(0.0I)
@test R[4] ≈ DCM(0.0I) 

@test_nowarn add_axes_bcrtod!(frames, :BCR, -2, -1; deriv=false)

R = rotation12(frames, :ICRF, :BCR, 12.345)

@test R[1] ≈ DCM(1.0I)
@test R[2] ≈ DCM(0.0I)
@test R[3] ≈ DCM(0.0I)
@test R[4] ≈ DCM(0.0I) 

@test_throws Exception add_axes_bci2000!(frames, :BCR, -2, -10)
@test_throws Exception add_axes_bcrtod!(frames, :BCR, -2, -10)