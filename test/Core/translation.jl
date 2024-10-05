using FrameTransformations
using StaticArrays
using Test

# SVector3
@test FrameTransformations.svector3_eltype(FrameTransformations.SVector3{Float64}) == Float64
vft1 = FrameTransformations.SVector3(1.,2.,3.)
vit1 = FrameTransformations.SVector3(4, 5, 6)
@test FrameTransformations.svector3_eltype(vft1) == Float64
@test FrameTransformations.svector3_eltype(vit1) == Int64
@test FrameTransformations.promote_svector3_eltype((vft1, vit1)) == Float64

# Basic ctor
t = Translation((vft1, vit1))
@test t[1] == vft1
@test t[2] == vit1
@test size(t) == (2,)
@test length(t) == 2
@test keys(t) == Base.OneTo(2)
@test order(t) == 2

# Varargs constructor
@test t == Translation(1,2,3,4,5,6.)
@test_throws DimensionMismatch Translation(1,2.)
# Convert to Tuple 
@test Tuple(t) == (1,2,3,4,5,6)
# Convert to SVector
@test SVector(t) == SA[1,2,3,4,5,6]
@test_throws DimensionMismatch Translation{1}(SA[1,2]) 

# Empty constructor
t0 = Translation{3, Float64}()
@test t0[1] == zeros(3)
@test t0[2] == zeros(3)
@test t0[3] == zeros(3)

# Convert to a different order auto-fill of missing SVector3
t2 = Translation{1}(t)
@test length(t2) == 1
@test t2[1] == vft1

t3 = Translation{3}(t2)
@test length(t3) == 3
@test_throws DimensionMismatch t3 == t2
@test t3[2] == zeros(3)
@test t3[3] == zeros(3)

t4 = Translation{1}(SA[1,2,3])
@test length(t4) == 1
@test t4[1] == [1,2,3]

t5 = Translation{2}(SA[1,2,3])
@test length(t5) == 2
@test t5[2] == zeros(3)

t6 = t4 - t4
@test t6[1] == zeros(3)

t7 = t5 - t4
@test length(t7) == 2
@test t7[1] == zeros(3)
@test t7[2] == zeros(3)
@test t4 - t5 == t5 - t4

t8 = t2 + t0
@test length(t8) == 3
@test t8[1] == [1,2,3]
@test t8[2] == zeros(3)
@test t8[3] == zeros(3)

t9 = t0 + t2 
@test t8 == t9

t10 = -t2
@test t10[1] == [-1, -2, -3]
