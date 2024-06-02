using FrameTransformations
using ReferenceFrameRotations
using LinearAlgebra
using Tempo
using Test

frames = FrameSystem{4, Float64}()
add_axes_icrf!(frames)
add_axes_mod!(frames, :MOD, 399)

e = Epoch("2022-01-01T12:00:00.0000 UTC")
@test_nowarn add_axes_frozen!(frames, :FROZEN, -399, 399, e)

eTDB = convert(TDB, e)
R0 = rotation3(frames, :ICRF, :MOD, eTDB)
R = rotation3(frames, :ICRF, :FROZEN, eTDB)

@test R == R0

R = rotation12(frames, :ICRF, :FROZEN, eTDB)

@test R[2] ≈ DCM(0.0I)
@test R[3] ≈ DCM(0.0I)
@test R[4] ≈ DCM(0.0I) 
